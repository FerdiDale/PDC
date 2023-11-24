#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "ex2.h"
#include "ex3.h"

void prodottoMatVet(int *globalptr, const int myrow, const int mycol, const int myRank, const int nproc, 
                      const int blocks[2], const int globalsizes[2], const int localsizes[2],
                      int **localdata, MPI_Comm rowComm, MPI_Comm colComm, int*xVec, int*xVecLoc, int* yVec, int* yVecLoc, int strat);

void distributeGeneralVector(const int myrow, const int myRank,
                      const int blocks[2], const int globalsizes[2], const int localsizes[2],  MPI_Comm colComm, int*xVec, int* xVecLoc);

void printvec(int* vec, int size, int myRank);

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
    int* yVec;
    int* xVec;
    int* xVecLoc;
    int* yVecLoc;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    // p, q, num righe, num colonne, strategia
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

    switch(strat) {
        case 1: 
            createGrid(&grid, &gridRow, &gridCol, myRank, numProcesses, 1, numProcesses, coordinates);
            break;
        case 2: 
            createGrid(&grid, &gridRow, &gridCol, myRank, numProcesses, numProcesses, 1, coordinates);
            break;
        case 3: 
            createGrid(&grid, &gridRow, &gridCol, myRank, numProcesses, p, q, coordinates);
            break;
        default:
            perror("Numero di strategia errato");
            MPI_Finalize();
            exit(1);
    }

    // Crea la griglia 2D

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
        xVec = malloc(totalCols*sizeof(int));
        yVec = calloc(totalRows,sizeof(int));
    }
    localMatrix = allocint2darray(localRows, localCols);
    xVecLoc = malloc(localCols*sizeof(int));
    yVecLoc = calloc(localRows,sizeof(int));    

    // Inizializza la matrice sul processo radice (rank 0)
    if (myRank == 0) {
        for (i = 0; i < totalRows; i++) {
            for (j = 0; j < totalCols; j++) {
                matrix[i][j] = i * totalCols + j;
            }
        }
        for (i = 0; i < totalCols; i++) {
            xVec[i] = 1;
        }
        // Stampa la matrice totale
        printLocalMatrix(matrix, totalRows, totalCols, myRank);
    }

    int blocks[2] = {p, q};
    int globalsizes[2] = {totalRows, totalCols};
    int localsizes[2] = {localRows, localCols};
    
    // Distribuisci la matrice tra i processi
    if (myRank == 0)
        prodottoMatVet
    (&(matrix[0][0]), coordinates[0], coordinates[1], myRank, numProcesses, blocks, globalsizes, localsizes, localMatrix, gridRow, gridCol, xVec, xVecLoc, yVec, yVecLoc, strat);
    else 
        prodottoMatVet
    (NULL, coordinates[0], coordinates[1], myRank, numProcesses, blocks, globalsizes, localsizes, localMatrix, gridRow, gridCol,  NULL, xVecLoc, NULL, yVecLoc, strat);

    printvec(yVec, totalRows, myRank);

    // Deallocazione della memoria
    if (myRank == 0) {
        freeint2darray(matrix);
        free(xVec);
        free(yVec);
    }
    freeint2darray(localMatrix);
    free(xVecLoc);
    free(yVecLoc);
    MPI_Comm_free(&grid);
    MPI_Comm_free(&gridRow);
    MPI_Comm_free(&gridCol);

    MPI_Finalize();
    return 0;

}

void prodottoMatVet(int *globalptr, const int myrow, const int mycol, const int myRank, const int nproc, 
                      const int blocks[2], const int globalsizes[2], const int localsizes[2],
                      int **localdata, MPI_Comm rowComm, MPI_Comm colComm, int* xVec, int* xVecLoc, int* yVec, int* yVecLoc, int strat) {

    int i,j;

    int* sumVecLoc = calloc(localsizes[0], sizeof(int));

    distributeMatrix(globalptr, myrow, mycol, myRank, nproc, 
                      blocks, globalsizes, localsizes,
                      localdata, rowComm, colComm);

    distributeGeneralVector(myrow, myRank, blocks, globalsizes, localsizes, colComm, xVec, xVecLoc);

    for (i=0; i<localsizes[0]; i++){
        for (j=0; j<localsizes[1]; j++){
            yVecLoc[i] += localdata[i][j]* xVecLoc[j];
        }
    }

    MPI_Allreduce(yVecLoc, sumVecLoc, localsizes[0], MPI_INT, MPI_SUM, rowComm); 

    MPI_Allgather(sumVecLoc, localsizes[0], MPI_INT, yVec, localsizes[0], MPI_INT, colComm);

    free(sumVecLoc);

}

void distributeGeneralVector(const int myrow, const int myRank,
                      const int blocks[2], const int globalsizes[2], const int localsizes[2],  MPI_Comm colComm, int*xVec, int* xVecLoc) {
    if (myrow == 0) {
        distributeVector(xVec, globalsizes[1], xVecLoc, localsizes[1], myRank, blocks[0]);
    }
    MPI_Bcast(xVecLoc, localsizes[1], MPI_INT, 0, colComm);
}

void printvec(int* vec, int size, int myRank) {
    int i,j;
    char buf [1000];
    char addbuf [100];

    snprintf(buf, sizeof(buf), "Processo %d: Vettore risultante:\n", myRank);
    for (i = 0; i < size; i++) {
        snprintf(addbuf, sizeof(addbuf), "%4d ", vec[i]);
        strcat(buf,addbuf);
    }
    snprintf(addbuf, sizeof(addbuf), "\n");
    strcat(buf,addbuf);
    printf("%s", buf);
}
