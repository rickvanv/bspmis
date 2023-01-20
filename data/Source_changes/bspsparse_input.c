#include "bspedupack.h"
#include "bspsparse_input.h"

double bspinput2triple(char *filepath, long *pnA, long *pnz,
                       long **pia, long **pja, double **pa){
  
    /* This function reads a sparse matrix in distributed
       Matrix Market format from the input file and
       distributes matrix triples to the processors.
       The matrix must be square and contain numerical values
       (real or integer) or it must be a square pattern matrix.

       The input consists of a number of header lines starting with '%'
       followed by one line
           m n nz p  (number of rows, columns, nonzeros, processors)
       followed by p+1 lines with the starting numbers
       of the processor parts
           Pstart[0]
           Pstart[1]
           ...
           Pstart[p]
       which means that processor q will get all nonzeros
       numbered Pstart[q]..Pstart[q+1]-1.
       This is followed by nz lines in the format
           i j a     (row index, column index, numerical value),
       or in case of a pattern matrix
           i j       (row index, column index).
       Nonzeros from pattern matrices will be given the default value 1.0.

       The input indices are assumed by Matrix Market to start
       counting at one, but they are converted to start from zero.
       The triples are stored on the responsible processor
       into three local arrays ia, ja, a, in arbitrary order.
       
       Output:
       nA is the global matrix size.
       nz is the number of local nonzeros.
       a[k] is the numerical value of the k'th local nonzero,
            0 <= k < nz.
       ia[k] is the global row index of the k'th local nonzero.
       ja[k] is the global column index.

       The function called by P(0) returns the sum of the
       matrix values for checksum purposes. The others return 0.
    */

    long p= bsp_nprocs(); // p = number of processors obtained
    long s= bsp_pid(); // processor number
    double sum= 0.0;
    bool pattern= true; // pattern matrix (true) or real/integer matrix (false)

    /* Initialize file pointer to NULL */
    FILE *fp = NULL;
    
    /* Initialize data and register global variables */
    long *Pstart= vecalloci(p+1);
    long nzA; //number of nonzeros of A
    bsp_push_reg(pnA,sizeof(long));
    bsp_push_reg(pnz,sizeof(long));
    bsp_push_reg(&nzA,sizeof(long));

    bsp_size_t tagsize= sizeof(indexpair);
    bsp_set_tagsize(&tagsize);
    bsp_sync();

    if (s==0){
        /* Open the matrix file and read the header */
//        char filepath[] = "/home/rick/CLionProjects/ParallelAlgorithms/mis/data/ca-GrQc/ca-GrQc.mtx-P4";
        fp=fopen(filepath,"r");

        /* A is an mA by nA matrix with nzA nonzeros
           distributed over pA processors. */

        /* Read header */
        char line[MAX_LINE_LENGTH];
        char banner[STRLEN], object[STRLEN], format[STRLEN],
             field[STRLEN], symmetry[STRLEN];
        long mA, pA;

        fgets(line, MAX_LINE_LENGTH, fp);

        /* Read banner line if it exists */
        if (sscanf(line, "%s %s %s %s %s", banner, object,
                   format, field, symmetry) != 5 ||
            strcmp(banner, "%%MatrixMarket")!= 0 ){
                /* Not a banner line */
                pattern= false; // default is a real matrix
        } else {
            /* A banner line */
            if (strcmp(field, "pattern")==0)
                pattern= true;
            else
                pattern= false;
            fgets(line, MAX_LINE_LENGTH, fp);
        }


        /* Skip lines starting with '%' */
        while (line[0] == '%')
            fgets(line, MAX_LINE_LENGTH, fp);

        /* Read matrix size, number of nonzeros, number of processors */
        sscanf(line,"%ld %ld %ld %ld\n", &mA, pnA, &nzA, &pA);

        if(mA != *pnA)
            bsp_abort("Error: matrix is not square");
        if(pA != p)
            bsp_abort("Error: p not equal to p(A)\n"); 

        /* Read distribution sizes */
        for (long q=0; q<=p; q++)
            fscanf(fp,"%ld\n", &Pstart[q]);
        for (long q=0; q<p; q++){
            bsp_put(q,pnA,pnA,0,sizeof(long));
            bsp_put(q,&nzA,&nzA,0,sizeof(long));
            long nzq= Pstart[q+1]-Pstart[q];
            bsp_put(q,&nzq,pnz,0,sizeof(long));
        }
    }
    bsp_sync();

    /* Determine the number of supersteps for input of the nonzeros */
    long nsteps= (nzA%MAXSEND==0 ? nzA/MAXSEND : nzA/MAXSEND +1);

    /* Handle the nonzeros in batches of size MAXSEND, which saves buffer
       memory, at the expense of extra syncs. The buffer
       memory needed for communication of the input is at most
       2*MAXSEND integers and MAXSEND reals. */

    long length= *pnz+1; // number of local nonzeros 
                         // + 1 extra for an ICRS dummy nonzero
    *pa= vecallocd(length);
    *pia= vecalloci(length);  
    *pja= vecalloci(length);

    long q=0;     // destination processor
    long count=0; // number of nonzeros stored locally so far

    for (long step=0; step<nsteps; step++){      
        indexpair t;
        if (s==0){
            /* Read the nonzeros from the matrix file and
               send them to their destination */
            for (long k=step*MAXSEND; k<(step+1)*MAXSEND && k<nzA; k++){
                double value=1.0; // default for pattern matrices
                if (pattern)
                    fscanf(fp,"%ld %ld\n", &t.i, &t.j);
                else
                    fscanf(fp,"%ld %ld %lf\n", &t.i, &t.j, &value);
                /* Convert indices to range 0..n-1, 
                   assuming they were in 1..n */
                t.i--;
		t.j--;
                /* Find processor of nonzero k */
                while(k==Pstart[q+1])
                    q++;
                /* Send a triple to P(q), with Pstart[q] <= k < Pstart[q+1]. 
                   The tag is an index pair (i,j).
                   The payload is a numerical value */
                bsp_send(q,&t,&value,sizeof(double));
                sum += value;
            }
        }
        bsp_sync();
        
        /* Store the received nonzeros */
        bsp_size_t nbytes;
        bsp_nprocs_t nmessages;
        bsp_qsize(&nmessages,&nbytes);
        for (long k=count; k<count+nmessages; k++){
            bsp_size_t status;
            bsp_get_tag(&status,&t);
            (*pia)[k]= t.i;
            (*pja)[k]= t.j;
            bsp_move(*pa+k,sizeof(double));
        }
        count += nmessages;
    }


    if (s==0)
        fclose(fp);
    bsp_pop_reg(&nzA);
    bsp_pop_reg(pnz);
    bsp_pop_reg(pnA);
    vecfreei(Pstart);
    bsp_sync();

    return sum;
    
} /* end bspinput2triple */


long key(long i, long radix, long keytype){
   /* This function computes the key of an index i
      according to the keytype */
   
       if (keytype==DIV)
           return i/radix;
       else /* keytype=MOD */
           return i%radix;
           
} /* end key */


void sort(long n, long nz, long *ia, long *ja, long *a,
          long radix, long keytype){

   /* This function sorts the nonzero elements of an n by n
      sparse matrix A stored in triple format in arrays ia, ja, a.
      The sort is by counting. 

      If keytype=DIV, the triples are sorted by increasing value of
      ia[k] div radix.
      if keytype=MOD, the triples are sorted by increasing value of
      ia[k] mod radix.

      The sorting is stable: ties are decided so that the original
      precedences are maintained. For a complete sort by increasing
      index ia[k], this function should be called twice:
      first with keytype=MOD, then with keytype=DIV.
      
      Input: 
      n is the size of the matrix.
      nz is the number of nonzeros.
      a[k] is the numerical value of the k'th nonzero of the
           sparse matrix A, 0 <= k < nz. Note: this value is assumed
           to be an integer. If numerical values are doubles,
           the array a can be used for registering the permutation of
           the nonzeros and performing a copy of the doubles afterwards. 
      ia[k] is the row index of the k'th nonzero, 0 <= ia[k] < n.
      ja[k] is the column index of the k'th nonzero.
      The indices can be local or global.
      radix >= 1.
      
      Output: ia, ja, a in sorted order.
   */
   
   if (nz <= 0)
       return;

   long *ia1= vecalloci(nz);
   long *ja1= vecalloci(nz); 
   long *a1 = vecalloci(nz);
   
   /* Allocate bins */
   long nbins;
   if (keytype==DIV)
       nbins= (n%radix==0 ? n/radix : n/radix+1);
   else {
       // keytype=MOD
       nbins= radix;
   }
   long *startbin= vecalloci(nbins);
   long *lengthbin= vecalloci(nbins);
       
   /* Count the elements in each bin */
   for (long r=0; r<nbins; r++)
       lengthbin[r]= 0;
   for (long k=0; k<nz; k++){
       long r= key(ia[k],radix,keytype);
       lengthbin[r]++;
   }
    
   /* Compute the starting position of each bin */
   startbin[0]= 0;
   for (long r=1; r<nbins; r++)
       startbin[r]= startbin[r-1] + lengthbin[r-1];
       
   /* Enter the elements into the bins in temporary arrays (ia1,ja1,a1) */
   for (long k=0; k<nz; k++){
       long r= key(ia[k],radix,keytype);
       long k1= startbin[r];
       ia1[k1]= ia[k];
       ja1[k1]= ja[k];
       a1[k1] = a[k];
       startbin[r]++;
   }
  
   /* Copy the elements back to the orginal arrays */
   for (long k=0; k<nz; k++){
       ia[k]= ia1[k];
       ja[k]= ja1[k];
       a[k] = a1[k];
   }
   
   vecfreei(lengthbin);
   vecfreei(startbin);
   vecfreei(a1);
   vecfreei(ja1);
   vecfreei(ia1);
   
} /* end sort */


void triple2icrs(long n, long nz, long *ia,  long *ja, double *a,
                 long *pnrows, long *pncols,
                 long **prowindex, long **pcolindex, long **pstart){

    /* This function converts a sparse matrix A given in triple
       format with global indices into a sparse matrix in
       incremental compressed row storage (ICRS) format with 
       local indices. It also outputs a start array for the
       regular CRS format.

       The conversion needs time and memory O(nz + sqrt(n))
       on each processor, which is O(nz(A)/p + n/p + p).
       
       Input:
       n is the global size of the matrix.
       nz is the local number of nonzeros.
       a[k] is the numerical value of the k'th nonzero
            of the sparse matrix A, 0 <= k < nz.
            Note: this value is assumed to be a double.
       ia[k] is the global row index of the k'th nonzero.
       ja[k] is the global column index of the k'th nonzero.
  
       Output:
       nrows is the number of local nonempty rows
       ncols is the number of local nonempty columns
       rowindex[i] is the global row index of the i'th
                   local row, 0 <= i < nrows.
       colindex[j] is the global column index of the j'th
                   local column, 0 <= j < ncols.
       start[i] is the start of the i'th local row in arrays ia, ja, a.
       a[k] is the numerical value of the k'th local nonzero of the
            sparse matrix A, 0 <= k < nz. The array is sorted by
            increasing row index, ties being decided by 
            increasing column index.
       ia[k]= inc[k] is the increment in the local column index
           of the k'th local nonzero, compared to the column
           index of the (k-1)th nonzero if this nonzero is
           in the same row; otherwise, ncols is added
           to the difference.
           Here, the column index of the -1'th nonzero is taken as 0.
       ja[k] is the local column index.
   */
    
   /* radix is the smallest power of two >= sqrt(n).
      The div and mod operations are cheap for powers of two.
      A radix of about sqrt(n) minimizes memory and time. */

   long radix;
   for (radix=1; radix*radix<n; radix *= 2)
       ;
   
   /* Initialize bookkeeping permutation orig to identity permutation */
   long *orig= vecalloci(nz);
   for (long k=0; k<nz; k++)
       orig[k]= k; // original index of k'th nonzero

   /* Sort nonzeros by column index using radix sort */
   sort(n,nz,ja,ia,orig,radix,MOD); // ja goes first
   sort(n,nz,ja,ia,orig,radix,DIV);
   
   /* Count the number of local columns */
   long ncols= 0;
   long jglob_prev= DUMMY;
   for (long k=0; k<nz; k++){
       long jglob= ja[k];
       if(jglob!=jglob_prev)
           /* new column index */
           ncols++;
       jglob_prev= jglob;
   }
   long *colindex= vecalloci(ncols);
   
   /* Convert global column indices to local ones.
      Initialize colindex */
   long j= 0; // local index
   jglob_prev= DUMMY;
   for (long k=0; k<nz; k++){
       long jglob= ja[k];
       if(jglob!=jglob_prev){
           colindex[j]= jglob;
           j++; 
           // j= first free position
       }
       ja[k]= j-1; // last entered local index
       jglob_prev= jglob;
   }
   
   /* Sort nonzeros by row index using radix sort */
   sort(n,nz,ia,ja,orig,radix,MOD); // ia goes first
   sort(n,nz,ia,ja,orig,radix,DIV);

   /* Permute the numerical values */
   double *a1= vecallocd(nz);

   for (long k=0; k<nz; k++)
       a1[k]= a[orig[k]];

   for (long k=0; k<nz; k++)
       a[k]= a1[k];

   vecfreed(a1);
   vecfreei(orig);

   /* Count the number of local rows */
   long nrows= 0;
   long iglob_prev= DUMMY;
   for (long k=0; k<nz; k++){
       long iglob= ia[k];
       if(iglob!=iglob_prev)
           /* new row index */
           nrows++;
       iglob_prev= iglob;
   }
   long *rowindex= vecalloci(nrows);
   long *start= vecalloci(nrows+1);
                              
   /* Convert global row indices to local ones.
      Initialize rowindex and inc */
   long i= 0; // local index
   iglob_prev= DUMMY;
   for (long k=0; k<nz; k++){
       long inck;
       if (k==0)
           inck= ja[k];
       else
           inck= ja[k] - ja[k-1];

       long iglob= ia[k]; 
       if(iglob!=iglob_prev){
           if(k>0)
               inck += ncols;
           rowindex[i]= iglob;
           start[i]= k;
           i++;
       } 
       ia[k]= inck; // ia is used to store inc
       iglob_prev= iglob;
   }
   start[nrows]= nz; // add dummy row for CRS

   /* Add dummy triple for ICRS */
   if (nz==0)
       ia[nz]= ncols;
   else 
       ia[nz]= ncols - ja[nz-1];
   ja[nz]= 0;                                                  
   a[nz]= 0.0;     
   
   *pnrows= nrows;
   *pncols= ncols;
   *prowindex= rowindex;
   *pcolindex= colindex;
   *pstart= start;
   
} /* end triple2icrs */
