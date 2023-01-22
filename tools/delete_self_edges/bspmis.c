#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "bspedupack.h"
#include "bspsparse_input.h"

#define MAX_WORD_LENGTH 1024
#define DEADEDGE 0
#define DEADVERTEX 1
long P;
//char *mtxfilepath;
char distrmtxfilepath[MAX_WORD_LENGTH] = "/home/rick/CLionProjects/ParallelAlgorithms/mis/data/testmatrix/test.mtx-P2";

void bsp_process_recvd_msgs(bool *Alive_edge, bool *Alive_vertex, long const *Adj, long const *destproc, long const *v1, long const *Start, long const *v0, long *nalive, long const nedges){
    bsp_nprocs_t nmessages; // total number of messages received
    bsp_size_t nbytes; //size of nmessages
    bsp_qsize(&nmessages,&nbytes);

    for (long i=0 ; i<nmessages; i++){
        bsp_size_t status; // not used
        long tag;
        bsp_get_tag(&status, &tag);
        long e;
        bsp_move(&e, sizeof(long));
        long v = v0[e]; //local vertex of halo edge
        if (tag == DEADVERTEX){
            if (Alive_vertex[v]) {
                Alive_vertex[v] = false;
                (*nalive)--;
            }
            for (long k=Start[v]; k<Start[v+1]; k++){
                long f = Adj[k];
                Alive_edge[f] = false;
                if (f>nedges){
                    long alive_tag = DEADEDGE; //communicate to remote processor that edge is dead
                    bsp_send(destproc[e - nedges], &alive_tag, &(v1[e]), sizeof(long));
                }
            }
        }
        if (tag == DEADEDGE)
            Alive_edge[e] = false;
    }
}

/* TODO: Insert Mondriaan part */
//void write_distributed_mtx(char *inputfilepath){
//    sprintf(distrmtxfilepath, "%s-P%ld", inputfilepath, P);
//    FILE *OutputFile;
//    if ((OutputFile = fopen(distrmtxfilepath,"r")) == NULL) { //check if output file already exists
//        OutputFile = fopen(distrmtxfilepath,"w");
//        FILE *InputFile;
//        InputFile = fopen(inputfilepath, "r");
//        /* This will contain the Mondriaan options. */
//        struct opts Options;
//        /* This structure will contain the input matrix. */
//        struct sparsematrix inputmatrix;
//        /* Set the default options. */
//        SetDefaultOptions(&Options);
//        /* If we are done setting the options, we check and apply them. */
//        if (!ApplyOptions(&Options)) {
//            printf("Invalid options!\n");
//        }
//        /* Read it from the file. */
//        if (!MMReadSparseMatrix(InputFile, &inputmatrix)) {
//            printf("Unable to read matrix!\n");
//            fclose(InputFile);
//        }
//        fclose(InputFile);
//        if (!DistributeMatrixMondriaan(&inputmatrix, P, 0.03, &Options, NULL)) {
//            printf("Unable to distribute matrix!\n");
//        }
//
//        /* Write the distributed matrix to file */
//        MMWriteSparseMatrix(&inputmatrix, OutputFile, NULL, &Options);
//        fclose(OutputFile);
//    }
//}

void bspmis(){

    bsp_begin(P);

    /***** Part 0: prepare input *****/

    long p= bsp_nprocs(); // p=P
    long s= bsp_pid();

    /* Input of sparse matrix into triple storage */
    long n, nz, *ia, *ja;
    double *weight;
    double suma= bspinput2triple(distrmtxfilepath,&n,&nz,&ia,&ja,&weight);/* Sequential part */
    long maxops=0;
    bsp_push_reg(&maxops,sizeof(long));
    bsp_sync();

//    if (s==0){
//        printf("Please enter the maximum number of operations per superstep\n");
//        printf("    (0 if no maximum)\n"); fflush(stdout);
//        scanf("%ld",&maxops);
//        if (maxops>0)
//            printf("Maximum number of operations per superstep = %ld\n\n", maxops);
//        else
//            printf("No maximum number of operations per superstep\n\n");
//        for (long t=0; t<p; t++)
//            bsp_put(t,&maxops,&maxops,0,sizeof(long));
//    }

    /* Convert data structure to incremental compressed row storage (ICRS),
       where ia contains the local column index increments,
       ja the local column indices, and Start the starting points of rows. */

    long nrows, ncols, *rowindex, *colindex, *Start;
    triple2icrs(n,nz,ia,ja,weight,&nrows,&ncols,&rowindex,&colindex,&Start);

    vecfreei(ia); // increments are not needed

    /* Translate to graph language. Here, nz is the number of edges
       including symmetric duplicates. */
    long nvertices= nrows; // number of local vertices (with degree > 0)
    long *v0= vecalloci(nz);
    long *v1= vecalloci(nz);
    long *weight1= vecalloci(nz);
    long *degree= vecalloci(nvertices);

    /* Search for the local row (vertex) corresponding to each global column.
      We use that both rowindex and colindex are ordered by increasing
      global index. */

    long *rowvertex= vecalloci(ncols);
    long r=0; //local row index
    for (long j=0; j<ncols; j++){
        long jglob= colindex[j];
        while (r<nrows && rowindex[r]<jglob)
            r++;
        if (r<nrows && rowindex[r]==jglob)
            rowvertex[j]= r;     // local vertex
        else
            rowvertex[j]= DUMMY; // nonlocal vertex
    }

    /* Initialize v0, v1, weight1, degree */
    long nhalo= 0; // number of halo edges
    double sum_max= 0.0; // sum of the maximum weights of the local rows
    for (long i=0; i<nrows; i++){
        degree[i]= Start[i+1] - Start[i];
        double maxw= 0.0; // maximum weight of row i
        for (long k=Start[i]; k<Start[i+1]; k++){
            v0[k]= i; // local index of row vertex
            v1[k]= rowvertex[ja[k]]; // local index of column vertex
            if (v1[k] == DUMMY){ // halo edge
                nhalo++;
                weight1[k]= rowindex[i] + colindex[ja[k]];
            } else
                weight1[k]= 2*n + rowindex[i] + colindex[ja[k]];
            if (weight[k] > maxw)
                maxw= weight[k];
        }
        sum_max += maxw;
    }
    if ((nz - nhalo)%2==1)
        bsp_abort("Error on input: nz-nhalo is odd\n");
    long nedges= (nz - nhalo)/2;

    /* Determine new edge numbers */
    long *new=vecalloci(nz);
    long count= 0;
    for (long k=0; k<nz; k++){
        if (v1[k]!=DUMMY && v0[k]<v1[k]){ // local edges are registered first,
            new[k]= count;                 // and only once
            count++;
        }
    }
    for (long k=0; k<nz; k++){
        if (v1[k]==DUMMY){ // halo edges second
            new[k]= count;
            count++;
        }
    }

    /* Insert edges into adjacency lists */
    long *Adj= vecalloci(nz);
    long *free= vecalloci(nrows);

    for (long i=0; i<nrows; i++)
        free[i]= Start[i]; // first free position for row i

    for (long k=0; k<nz; k++){
        if (v1[k]==DUMMY){
            Adj[free[v0[k]]]= new[k];
            free[v0[k]]++;
        } else if (v0[k] < v1[k]){ // insert edge in two lists
            Adj[free[v0[k]]]= new[k]; free[v0[k]]++;
            Adj[free[v1[k]]]= new[k]; free[v1[k]]++;
        }
    }
    vecfreei(free);

    /*  Copy the edges according to the new numbering */
    long nedges_tot= nedges + nhalo;
    long *v0new= vecalloci(nedges_tot);
    long *v1new= vecalloci(nedges_tot);
    long *janew= vecalloci(nedges_tot);
    double *weightnew= vecallocd(nedges_tot);
    long *weight1new= vecalloci(nedges_tot);

    for (long k=0; k<nz; k++){
        if (v1[k]==DUMMY || v0[k] < v1[k]){
            v0new[new[k]]= v0[k];
            janew[new[k]]= ja[k]; // keep a copy for output printing
            if (v1[k]==DUMMY){
                v1new[new[k]]= ja[k]; // register the local column index
            } else
                v1new[new[k]]= v1[k];
            weightnew[new[k]]= weight[k];
            weight1new[new[k]]= weight1[k];
        }
    }

    vecfreei(new);
    vecfreei(weight1);
    vecfreed(weight);
    vecfreei(v1);
    vecfreei(v0);
    /* Initialize communication data structure */
    long np= nloc(p,s,n);
    long *tmpproc=vecalloci(np); // temporary array for storing the owners
    // of the vertices
    bsp_push_reg(tmpproc,np*sizeof(long));

    /* Initialize owner to 0 as a default for vertices with degree 0 */
    for (long i=0; i<np; i++)
        tmpproc[i]= 0;

    /* Set tagsize for communication of halo edge numbers */
    bsp_size_t tagsize= sizeof(indextriple);
    bsp_set_tagsize(&tagsize);
    bsp_sync();

    /* Announce my vertices (rows) */
    for (long i=0; i<nrows; i++){
        long iglob= rowindex[i];
        bsp_put(iglob%p,&s,tmpproc,(iglob/p)*sizeof(long),sizeof(long));
    }
    bsp_sync();

    /* Determine the owner proc[j] of vertex colindex[j],
       for each local column j, by reading the announcements */
    long *proc= vecalloci(ncols);
    for (long j=0; j<ncols; j++){
        if (rowvertex[j]==DUMMY){
            long jglob= colindex[j];
            bsp_get(jglob%p,tmpproc,(jglob/p)*sizeof(long),
                    &(proc[j]),sizeof(long));
        } else
            proc[j]= s;
    }
//    vecfreei(rowvertex);
    bsp_sync();
    bsp_pop_reg(tmpproc);

    /* Initialize destproc[e-nedges]= owner of halo edge e
       and send the local edge number e to this owner */
    long *destproc= vecalloci(nhalo);
    indextriple tag;
    for (long e=nedges; e<nedges_tot; e++){
        destproc[e-nedges]= proc[v1new[e]]; // local column index
        tag.iloc = v0new[e]; // local row index
        tag.i= rowindex[v0new[e]]; // global row index
        tag.j= colindex[v1new[e]]; // global column index
        bsp_send(destproc[e-nedges], &tag, &e, sizeof(long));
    }
    bsp_sync();

    vecfreei(proc);
    vecfreei(tmpproc);

    /* Receive triples (i,j,e), where i,j are global indices and
       e is the local edge number on the sending processor */

    bsp_nprocs_t nmsg; // total number of messages received
    bsp_size_t nbytes;    // total size in bytes received
    bsp_qsize(&nmsg,&nbytes);
    if (nmsg != nhalo)
        bsp_abort("Error: number of messages <> nhalo\n");

    long *ilocmsg= vecalloci(nmsg);
    long *imsg= vecalloci(nmsg);
    long *jmsg= vecalloci(nmsg);
    long *emsg= vecalloci(nmsg);
    for (long k=0; k<nmsg; k++){
        bsp_size_t status; // not used
        bsp_get_tag(&status, &tag);
        ilocmsg[k] = tag.iloc;
        imsg[k]= tag.i;
        jmsg[k]= tag.j;
        bsp_move(&(emsg[k]), sizeof(long));
    }

    /* Sort nmsg edges by primary key j and secondary key i,
       using as radix the smallest power of two >= sqrt(n).
       The div and mod operations are cheap for powers of two.
       A radix of about sqrt(n) minimizes memory and time. */

    long radix;
    for (radix=1; radix*radix<n; radix *= 2)
        ;

    /* Duplicate imsg and jmsg for sorting ilocmsg */
    long *imsgcopy = vecalloci(nmsg);
    long *jmsgcopy = vecalloci(nmsg);
    for (long k=0; k<nmsg; k++){
        imsgcopy[k] = imsg[k];
        jmsgcopy[k] = jmsg[k];
    }

    sort(n,nmsg,imsg,jmsg,emsg,radix,MOD); // imsg goes first
    sort(n,nmsg,imsg,jmsg,emsg,radix,DIV);
    sort(n,nmsg,jmsg,imsg,emsg,radix,MOD); // jmsg goes first
    sort(n,nmsg,jmsg,imsg,emsg,radix,DIV);

    sort(n,nmsg,imsgcopy,jmsgcopy,ilocmsg,radix,MOD); // imsg goes first
    sort(n,nmsg,imsgcopy,jmsgcopy,ilocmsg,radix,DIV);
    sort(n,nmsg,jmsgcopy,imsgcopy,ilocmsg,radix,MOD); // jmsg goes first
    sort(n,nmsg,jmsgcopy,imsgcopy,ilocmsg,radix,DIV);

    /* Couple the local halo edges with the remote halo edges.
       The local halo edges e = (i,j) have been sorted
       in the CRS data structure by primary key i and secondary key j,
       with i and j global indices.

       A remote edge (j,i) corresponds to a local edge (i,j).
       For this reason, the received edges have been sorted
       by primary key j and secondary key i. */


    long *destvertex= vecalloci(nhalo);
    for (long e=nedges; e<nedges_tot; e++) {
        v1new[e] = emsg[e - nedges];
        destvertex[e - nedges] = ilocmsg[e-nedges];
    }

    vecfreei(ilocmsg);
    vecfreei(emsg);
    vecfreei(jmsg);
    vecfreei(imsg);

    long nmatch= 0;          // number of matches found
    long *match= vecalloci(nvertices); // matches found
    long nsteps= 0;  // number of (mixed) supersteps taken
    long nops= 0;    // number of elemental operations carried out

    /***** Part 1: run misfinder *****/
    double time0= bsp_time();

    /* Initialize local mis array */
    long *locmis = vecalloci(nvertices*sizeof(long));
    long miscount = 0;

    /* Initialize array for storing alive value of remote vertices */
    long nalive = nvertices;
    bool *Alive_edge = vecallocb(nedges_tot*sizeof(bool));
    for(long e=0; e<nedges_tot; e++){
        Alive_edge[e] = true;
    }
    bool *Alive_vertex = vecallocb(nvertices*sizeof(bool));
    for(long v=0; v<nvertices; v++){
        Alive_vertex[v] = true;
    }

    /* Initialize array for storing remote vertex rowindex */
    long *v0newrem = vecalloci(nhalo * sizeof(long));
    for (long e = 0; e < nhalo; e++) { //initialize array
        v0newrem[e] = -1;
    }

    /* Initialize array for storing random values of vertexes */
    long *randval_v0 = vecalloci(nvertices);
    long *randval_v1 = vecalloci(ncols);
    bsp_push_reg(v0new, nedges_tot * sizeof(long)); // TODO: set nedges_tot to max(nedges_tot) over all processors (read from mondriaan MTX file?)
    bsp_push_reg(randval_v0, nvertices * sizeof(long));

    long *Done= vecalloci(p);
    bsp_push_reg(Done,p*sizeof(long));

    bsp_sync();

    bool alldone = false;
    srand(time(NULL) * (s + 1));
    while (!alldone) {
        /* Initialize all processors to not done yet */
        for (long t=0; t<p; t++)
            Done[t]= false;
        if (nalive == 0){
            long done = true;
            for (long t=0; t<p; t++)
                bsp_put(t,&done,Done,s*sizeof(long),sizeof(long));
        }
        if (!Done[s]) {
            bsp_process_recvd_msgs(Alive_edge, Alive_vertex, Adj, destproc, v1new, Start, v0new, &nalive, nedges);
            bsp_sync();
            /* Create array of random values of (local) row vertices */
            for (long v = 0; v < nvertices; v++) {
                if (Alive_vertex[v])
                    randval_v0[v] = rand();
            }
            /* Align random values of row vertices with column vertices */
            for (long k = 0; k < ncols; k++) {
                if (rowvertex[k] != DUMMY)
                    randval_v1[k] = randval_v0[rowvertex[k]];
                else
                    randval_v1[k] = -1;
            }

//            /* Retrieve local rowindex of remote vertex */
//            for (long e = nedges; e < nedges_tot; e++) { // TODO: reduce into a single loop and a single get?
//                if (Alive_edge[e]) {
//                    bsp_get(destproc[e - nedges], v0new, v1new[e] * sizeof(long), &(v0newrem[e - nedges]),
//                            sizeof(long));
//                }
//            }
            bsp_sync();

            /* Using local rowindex of remote vertex, retrieve corresponding random value */
            for (long e = nedges; e < nedges_tot; e++) {
                if (Alive_edge[e]) {
                    bsp_get(destproc[e - nedges], randval_v0, destvertex[e - nedges] * sizeof(long),
                            &(randval_v1[janew[e]]),
                            sizeof(long));
                }
            }
            bsp_sync();

            for (long v = 0; v < nvertices ; v++) {
                if (Alive_vertex[v]) {
                    bool ismax = true;
                    for (long j = Start[v]; j < Start[v + 1]; j++) {
                        printf("%ld: v = %ld, v0 = %ld, v1 = %ld, colindex[janew[Adj[j]] = %ld\n",s, rowindex[v], rowindex[v0new[Adj[j]]], v1new[Adj[j]], colindex[janew[Adj[j]]]);
                        if (Alive_edge[Adj[j]] == true && randval_v0[v] <= randval_v1[janew[Adj[j]]]) {
                            ismax = false;
                        }
                    }
                    if (ismax) {
                        locmis[miscount] = rowindex[v];
                        miscount++;
                        Alive_vertex[v] = false; // set alive value of mis vertex to false
                        nalive--;
                        for (long k = Start[v];
                             k < Start[v + 1]; k++) { // set alive value of adjacent edges and vertices to false
                            long e = Adj[k];
                            Alive_edge[e] = false;
                            if (e < nedges) { // check if edge is a local edge
                                long w;
                                if (v0new[e] == v)
                                    w = v1new[e];
                                else
                                    w = v0new[e];
                                if (Alive_vertex[w]) {
                                    Alive_vertex[w] = false;
                                    nalive--;
                                }
                                for (long l = Start[w]; l < Start[w + 1]; l++) {
                                    e = Adj[l];
                                    if (e < nedges)
                                        Alive_edge[e] = false;
                                    else {
                                        long alive_tag = DEADEDGE; //communicate to remote processor that edge is dead
                                        bsp_send(destproc[e - nedges], &alive_tag, &(v1new[e]), sizeof(long));
                                    }
                                }
                            } else {
                                long alive_tag = DEADVERTEX; // communicate to remote processor that vertex is dead
                                bsp_send(destproc[e - nedges], &alive_tag, &(v1new[e]), sizeof(long));
                            }
                        }
                    }
                }
            }
            alldone = true;
            for (long t = 0; t < p; t++) {
                if (Done[t] == false) {
                    alldone = false;
                    break;
                }
            }
        }
        bsp_sync();
    }
    printf("Local MIS at processor %ld :", s);
    for (long k=0; k<miscount; k++){
        printf("%ld, ", locmis[k]+1);
    }
    fflush(stdout);

    bsp_pop_reg(v0new);
    bsp_pop_reg(randval_v0);
    bsp_sync();
    double time1 = bsp_time();
    if (s==0){
        printf("\nFound a MIS in only %.6lf seconds.\n\n",time1-time0);
    }
    bsp_end();
}

int main(int argc, char **argv){
    bsp_init(bspmis, argc, argv);
    /* Sequential part */
    P = 2;
//    mtxfilepath = "/home/rick/CLionProjects/ParallelAlgorithms/mis/data/ca-GrQc/ca-GrQc.mtx";
//    write_distributed_mtx(mtxfilepath);
    /* SPMD part */
    bspmis();
    return 0;
}


