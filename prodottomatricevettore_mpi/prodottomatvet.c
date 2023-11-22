#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "ex2.h"
#include "ex3.h"

void prodottomatvet(int *globalptr, const int myrow, const int mycol, const int myRank, const int nproc, 
                      const int blocks[2], const int globalsizes[2], const int localsizes[2],
                      int **localdata, MPI_Comm rowComm, MPI_Comm colComm, int*xvec, int** result, int strat);

int main (int argc, char* argv[]) {

    int myRank, numProcesses;
    int **matrix, **localMatrix;
    int totalRows; // Numero totale di righe nella matrice
    int totalCols; // Numero totale di colonne nella matrice
    int localRows, localCols;
    MPI_Comm grid, gridRow, gridCol;
    int coordinates[2];
    int i,j;
    int p,q;
    int strat;
    int* resultsvec;
    int* xvec;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    if (argc != 6) {
        perror("Numero di argomenti inseriti errato\n");
        MPI_Finalize();
        return 1;
    }

    p = atoi(argv[1]);
    q = atoi(argv[2]);
    totalRows = atoi(argv[3]);
    totalCols = atoi(argv[4]);
    strat = atoi(argv[5]);

    if (numProcesses != p * q) {
        if (myRank == 0) {
            printf("Il numero di processi deve formare una griglia pxq.\n");
        }
        MPI_Finalize();
        return 1;
    }

    // Crea la griglia 2D
    createGrid(&grid, &gridRow, &gridCol, myRank, numProcesses, p, q, coordinates);

    // Calcola le dimensioni locali della matrice per ogni processo
    localRows = totalRows / p;
    if (coordinates[0] < totalRows % p) {
        localRows++;
    }

    localCols = totalCols / q;
    if (coordinates[1] < totalCols % q) {
        localCols++;
    }

    // Alloca spazio per la matrice globale e la matrice locale
    if (myRank == 0) {
        matrix = allocint2darray(totalRows, totalCols);
    }
    localMatrix = allocint2darray(localRows, localCols);

    xvec = malloc(totalCols*sizeof(int));
    resultsvec = malloc(totalRows*sizeof(int));

    // Inizializza la matrice sul processo radice (rank 0)
    if (myRank == 0) {
        for (i = 0; i < totalRows; i++) {
            for (j = 0; j < totalCols; j++) {
                matrix[i][j] = i * totalCols + j;
            }
        }
        // Stampa la matrice totale
        printLocalMatrix(matrix, totalRows, totalCols, myRank);
    }

    int blocks[2] = {p, q};
    int globalsizes[2] = {totalRows, totalCols};
    int localsizes[2] = {localRows, localCols};
    
    // Distribuisci la matrice tra i processi
    if (myRank == 0)
        prodottomatvet(&(matrix[0][0]), coordinates[0], coordinates[1], myRank, numProcesses, blocks, globalsizes, localsizes, localMatrix, gridRow, gridCol, xvec, &resultsvec, strat);
    else 
        prodottomatvet(NULL, coordinates[0], coordinates[1], myRank, numProcesses, blocks, globalsizes, localsizes, localMatrix, gridRow, gridCol, xvec, &resultsvec, strat);


    printvec(resultsvec, totalRows, myRank);

    // Deallocazione della memoria
    if (myRank == 0) {
        freeint2darray(matrix);
    }
    freeint2darray(localMatrix);
    MPI_Comm_free(&grid);
    MPI_Comm_free(&gridRow);
    MPI_Comm_free(&gridCol);

    MPI_Finalize();
    return 0;

}

void prodottomatvet(int *globalptr, const int myrow, const int mycol, const int myRank, const int nproc, 
                      const int blocks[2], const int globalsizes[2], const int localsizes[2],
                      int **localdata, MPI_Comm rowComm, MPI_Comm colComm, int* xvec, int** result, int strat) {

    

}

void printvec(int* vec, int size, int myRank) {
    int i,j;
    char buf [1000];
    char addbuf [100];

    snprintf(buf, sizeof(buf), "Processo %d: Vettore risultante:\n", myRank);
    for (i = 0; i < localRows; i++) {
        snprintf(addbuf, sizeof(addbuf), "%4d ", vec[i]);
        strcat(buf,addbuf);
    }
    snprintf(addbuf, sizeof(addbuf), "\n");
    strcat(buf,addbuf);
    printf("%s", buf);
}