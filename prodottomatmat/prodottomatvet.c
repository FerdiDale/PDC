#include <mpi.h>
#include "prodottomatvet.h"
void prodottoMatVet(int* globalptr, const int myrow, const int mycol, const int myRank, const int nproc,
    const int blocks[2], const int globalsizes[2], const int localsizes[2],
    int** localdata, MPI_Comm rowComm, MPI_Comm colComm, int* xVec, int* xVecLoc, int* yVec, int* yVecLoc);

void distributeGeneralVector(const int myrow, const int mycol, const int myRank,
    const int blocks[2], const int globalsizes[2], const int localsizes[2], MPI_Comm rowComm, MPI_Comm colComm, int* xVec, int* xVecLoc);

void printvec(int* vec, int size, int myRank);

void gatherVector(int* locVec, int* totalVec, MPI_Comm communicator, int nsend, int nProc, int N);

/*
int main(int argc, char* argv[]) {

    int myRank, numProcesses;
    int** matrix = NULL;
    int** localMatrix = NULL;
    int totalRows; // Numero totale di righe nella matrice
    int totalCols; // Numero totale di colonne nella matrice
    int localRows, localCols;
    MPI_Comm grid, gridRow, gridCol;
    int coordinates[2];
    int i, j;
    int p, q;
    int strat;
    int* yVec = NULL;
    int* xVec = NULL;
    int* xVecLoc = NULL;
    int* yVecLoc = NULL;

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

    switch (strat) {
    case 1:
        p=numProcesses;
        q=1;
        break;
    case 2:
        p=1;
        q=numProcesses;
        break;
    case 3:
        break;
    default:
        perror("Numero di strategia errato");
        MPI_Finalize();
        exit(1);
    }

    createGrid(&grid, &gridRow, &gridCol, myRank, numProcesses, p, q, coordinates);

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
        xVec = calloc(totalCols, sizeof(int));
    }
    yVec = calloc(totalRows, sizeof(int));
    localMatrix = allocint2darray(localRows, localCols);
    xVecLoc = calloc(localCols, sizeof(int));
    yVecLoc = calloc(localRows, sizeof(int));

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
        printvec(xVec, totalCols, myRank);
    }

    int blocks[2] = { p, q };
    int globalsizes[2] = { totalRows, totalCols };
    int localsizes[2] = { localRows, localCols };

    // Distribuisci la matrice tra i processi
    if (myRank == 0)
        prodottoMatVet
        (&(matrix[0][0]), coordinates[0], coordinates[1], myRank, numProcesses, blocks, globalsizes, localsizes, localMatrix, gridRow, gridCol, xVec, xVecLoc, yVec, yVecLoc);
    else
        prodottoMatVet
        (NULL, coordinates[0], coordinates[1], myRank, numProcesses, blocks, globalsizes, localsizes, localMatrix, gridRow, gridCol, xVec, xVecLoc, yVec, yVecLoc);

    printvec(yVec, totalRows, myRank);

    // Deallocazione della memoria
    if (myRank == 0) {
        freeint2darray(matrix);
        free(xVec);
    }
    free(yVec);
    freeint2darray(localMatrix);
    free(xVecLoc);
    free(yVecLoc);
    MPI_Comm_free(&grid);
    MPI_Comm_free(&gridRow);
    MPI_Comm_free(&gridCol);

    MPI_Finalize();
    return 0;

}
*/

void prodottoMatVet(int* globalptr, const int myrow, const int mycol, const int myRank, const int nproc,
    const int blocks[2], const int globalsizes[2], const int localsizes[2],
    int** localdata, MPI_Comm rowComm, MPI_Comm colComm, int* xVec, int* xVecLoc, int* yVec, int* yVecLoc) {

    int i, j;

    int* sumVecLoc = calloc(localsizes[0], sizeof(int));

    // printf("\nProcesso=%d riga %d colonna %d : Arrivo a prima di distribuzione matrice\n", myRank, myrow, mycol);

    distributeMatrix(globalptr, myrow, mycol, myRank, nproc,
        blocks, globalsizes, localsizes,
        localdata, rowComm, colComm);

    // printf("\nProcesso=%d riga %d colonna %d dim vec %d: Arrivo a prima di distribuzione vettore\n", myRank, myrow, mycol, localsizes[1]);

    distributeGeneralVector(myrow, mycol, myRank, blocks, globalsizes, localsizes, rowComm, colComm, xVec, xVecLoc);

    for (i = 0; i < localsizes[0]; i++) {
        for (j = 0; j < localsizes[1]; j++) {
            yVecLoc[i] += localdata[i][j] * xVecLoc[j];
        }
    }

    MPI_Allreduce(yVecLoc, sumVecLoc, localsizes[0], MPI_INT, MPI_SUM, rowComm);

    gatherVector(sumVecLoc, yVec, colComm, localsizes[0], blocks[0], globalsizes[0]);

    free(sumVecLoc);

}

void gatherVector(int* locVec, int* totalVec, MPI_Comm communicator, int nsend, int nProc, int N) {

    int *recvcounts, *displs;
    int i;

    // Allocazione di array per determinare il numero di elementi da inviare a ciascun processo e gli spostamenti
    recvcounts = (int *)malloc(nProc * sizeof(int));
    displs = (int *)malloc(nProc * sizeof(int));

    // Calcola il numero medio di elementi per processo e gli elementi rimanenti
    int avg_elements = N / nProc;
    int remaining_elements = N % nProc;

    // Calcola il numero di elementi da inviare a ciascun processo e gli spostamenti
    for (i = 0; i < nProc; i++) {
        recvcounts[i] = avg_elements;
        if (i < remaining_elements) {
            recvcounts[i]++;
        }

        displs[i] = (i > 0) ? (displs[i - 1] + recvcounts[i - 1]) : 0;
    }

    MPI_Allgatherv(locVec, nsend, MPI_INT, totalVec, recvcounts, displs, MPI_INT, communicator);

    // Libera la memoria allocata per gli array di invio e spostamento
    free(recvcounts);
    free(displs);

}

void distributeGeneralVector(const int myrow, const int mycol, const int myRank,
    const int blocks[2], const int globalsizes[2], const int localsizes[2], MPI_Comm rowComm, MPI_Comm colComm, int* xVec, int* xVecLoc) {
    if (myrow == 0) {
        distributeVector(xVec, globalsizes[1], xVecLoc, localsizes[1], mycol, blocks[1], rowComm);
        // printf("\nProcesso=%d : Arrivo a dopo di distribuzione vettore\n", myRank);
        printvec(xVecLoc, localsizes[1], myRank);
    }
    MPI_Bcast(xVecLoc, localsizes[1], MPI_INT, 0, colComm);
}

void printvec(int* vec, int size, int myRank) {
    int i;
    char buf[1000];
    char addbuf[100];

    snprintf(buf, sizeof(buf), "Processo %d: Vettore risultante:\n", myRank);
    for (i = 0; i < size; i++) {
        snprintf(addbuf, sizeof(addbuf), "%4d ", vec[i]);
        strcat(buf, addbuf);
    }
    snprintf(addbuf, sizeof(addbuf), "\n");
    strcat(buf, addbuf);
    printf("%s", buf);
}
