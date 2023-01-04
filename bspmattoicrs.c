#include <stdio.h>
#include <stdlib.h>
#include "math.h"
#include "bspedupack.h"
#include "bspsparse_input.h"
//#include <Mondriaan.h>

long P;

struct edge{
    long v0;
    long v1;
};

//struct edge* load_graphdata(int s){
//    FILE *graphfile = fopen("/home/rick/CLionProjects/ParallelAlgorithms/mis/CA-GrQc.txt", "r");
//    if (graphfile == NULL){
//        printf("File not found\n");
//    }
//    long vertices, edges;
//    fscanf(graphfile, "# Nodes: %ld Edges: %ld", &vertices, &edges);
//    fscanf(graphfile, "%*[^\n]"); //skip next line
//    printf("Nodes: %ld Edges: %ld", vertices, edges);
//    long blocksize = ceil((double)vertices/P);
//    struct edge* graphdata = malloc(edges*sizeof(*graphdata));
//    long edgecount = 0;
//    for (long i=0; i<edges; i++){
//        long v0, v1;
//        fscanf(graphfile, "%ld %ld", &v0, &v1);
//        if (v0 >= s*blocksize && v0 < (s+1)*blocksize){
//            if (v1 > v0 || v1 <s*blocksize){
//                graphdata[edgecount] = (struct edge){v0, v1};
//                edgecount++;
//            }
//
//        }
//    }
//    graphdata = realloc(graphdata, edgecount*sizeof(*graphdata));
//    return graphdata;
//}
void matrix_to_icrs(){
    bsp_begin(P);
    long p = bsp_nprocs();
    long s = bsp_pid();


    long n, nz, *ia, *ja;
    double *weight;
    printf("This is processor %ld\n", s);
    double suma= bspinput2triple(&n,&nz,&ia,&ja,&weight);

    long nrows, ncols, *rowindex, *colindex, *Start;
    triple2icrs(n,nz,ia,ja,weight,&nrows,&ncols,&rowindex,&colindex,&Start);

    vecfreei(ia); // increments are not needed

    printf("%ld :", s);
    for(int i=0; i<nrows; i++){
        printf("%ld, ", rowindex[i]);
    }
    bsp_end();
}

    /* Translate to graph language. Here, nz is the number of edges
       including symmetric duplicates. */
//    long nvertices= nrows; // number of local vertices (with degree > 0)
//    long *v0= vecalloci(nz);
//    long *v1= vecalloci(nz);
//    long *weight1= vecalloci(nz);
//    long *degree= vecalloci(nvertices);

int main(int argc, char **argv){
    bsp_init(matrix_to_icrs, argc, argv);
    /* Sequential part */
    printf("How many processors do you want to use?\n");
    fflush(stdout);

    scanf("%ld",&P);

    matrix_to_icrs();
    return 0;
}


