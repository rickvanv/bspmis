
double bspinput2triple(char *filepath, long *pnA, long *pnz,
                     long **pia, long **pja, double **pa);
  
void sort(long n, long nz, long *ia, long *ja, long *a,
          long radix, long keytype);

void triple2icrs(long n, long nz, long *ia,  long *ja, double *a,
                 long *pnrows, long *pncols,
                 long **prowindex, long **pcolindex, long **pstart);

#define DIV 0
#define MOD 1
#define STRLEN 100
#define MAX_LINE_LENGTH 500
#define MAXSEND 10000 // Maximum number of nonzeros to send
                      // in one superstep, to save buffer memory

typedef struct {long i,j;} indexpair;
typedef struct {long iloc,i,j;} indextriple;
