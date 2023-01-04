#include <stdlib.h>
#include <stdio.h>
/* Make sure to include the Mondriaan headers. */
#include <Mondriaan.h>

int main(int argc, char **argv)
{
    /* This will contain the Mondriaan options. */
    struct opts Options;
    /* This file pointer will be the opened Matrix Market file. */
    FILE *File;
    /* This structure will contain the arc130 matrix. */
    struct sparsematrix Arc130;
    /* Variables used for calculating the communication volume. */
    long ComVolumeRow, ComVolumeCol, Dummy;
    /* Variables used to calculate the imbalance. */
    long MaxNrNz;
    double Imbalance;

    /* Set the default options. */
    SetDefaultOptions(&Options);

    /* We can also read options from disk. */
    /* if (!SetOptionsFromFile(&Options, "options.txt")) printf("Unable to set options from disk!\n"); */

    /* If we are done setting the options, we check and apply them. */
    if (!ApplyOptions(&Options))
    {
        printf("Invalid options!\n");
        return EXIT_FAILURE;
    }

    /* Open the arc130 matrix file. */
    if (!(File = fopen("/home/rick/Documents/mondriaan-master/tests/arc130.mtx", "r")))
    {
        printf("Unable to open arc130!\n");
        return EXIT_FAILURE;
    }

    /* Read it from the file. */
    if (!MMReadSparseMatrix(File, &Arc130))
    {
        printf("Unable to read arc130!\n");
        fclose(File);
        return EXIT_FAILURE;
    }

    fclose(File);

    /* Distribute the matrix over two processors with an allowed imbalance of 3% and the options provided above. */
    if (!DistributeMatrixMondriaan(&Arc130, 2, 0.03, &Options, NULL))
    {
        printf("Unable to distribute arc130!\n");
        return EXIT_FAILURE;
    }

    /* Calculate the communication volume. */
    CalcCom(&Arc130, NULL, ROW, &ComVolumeRow, &Dummy, &Dummy, &Dummy, &Dummy);
    CalcCom(&Arc130, NULL, COL, &ComVolumeCol, &Dummy, &Dummy, &Dummy, &Dummy);

    /* Calculate the imbalance, making use of the fact that we only distributed the nonzeros over two processors. */
    MaxNrNz = MAX(Arc130.Pstart[2] - Arc130.Pstart[1], Arc130.Pstart[1] - Arc130.Pstart[0]);
    Imbalance = (double)(2*MaxNrNz - Arc130.NrNzElts)/(double)Arc130.NrNzElts;

    if (Imbalance > 0.03)
    {
        printf("Imbalance is too large!\n");
        return EXIT_FAILURE;
    }

    /* Display information about this partitioning. */
    printf("Succesfully distributed %ld nonzeros over two processors: %ld are assigned to processor 0 and %ld to processor 1.\n", Arc130.NrNzElts, Arc130.Pstart[1] - Arc130.Pstart[0], Arc130.Pstart[2] - Arc130.Pstart[1]);
    printf("This distribution has a total communication volume equal to %ld and imbalance equal to %.1f%%.\n", ComVolumeRow + ComVolumeCol, 100.0*Imbalance);

    /* Free matrix data. */
    MMDeleteSparseMatrix(&Arc130);

    return EXIT_SUCCESS;
}
