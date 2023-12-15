#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "prodottomatvet.h"
#include "ex2.h"
#include "ex3.h"

int main(int argc, char* argv[]) {

    int myRank, numProcesses;
    int** A = NULL;
    int** localA = NULL;
    int** broadcastA = NULL;
    int** B = NULL;
    int** localB = NULL;
    int** localC= NULL;
    int totalN; // Numero totale di righe nella matrice
    int localN;
    MPI_Comm grid, gridRow, gridCol;
    int coordinates[2];
    int i, j;
    int p;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    //N taglia della matrice
    if (argc != 2) {
        if (myRank == 0)
            fprintf(stderr, "%s", "Numero di argomenti inseriti errato\n");
        MPI_Finalize();
        return 1;
    }

    totalN = atoi(argv[1]);

    if ((numProcesses & (numProcesses-1)) != 0) { //numero di processi non potenza di 2, sfruttiamo il codice binario particolare delle potenze di 2, composto da tutti 0 eccetto per un 1
            //Se nroc è potenza di 2, nproc-1 sarà quindi composto da tutti 1 fino ad uno 0 laddove nproc aveva il suo unico 1. L'and bit a bit quindi ci restituirà 0.
        if (myRank == 0)
            fprintf(stderr, "%s", "Il numero di processi deve essere un quadrato perfetto\n");
        MPI_Finalize();
        return 1;
    }

    p = sqrt(numProcesses)

    if (totalN % p != 0) {
        if (myRank == 0)
            fprintf(stderr, "%s", "La taglia della matrice quadrata deve essere multiplo della taglia della griglia di processi\n");
        MPI_Finalize();
        return 1;
    }

    createGrid(&grid, &gridRow, &gridCol, myRank, numProcesses, p, p, coordinates);

    // Crea la griglia 2D

    // Calcola le dimensioni locali delle matrici (N/p) per ogni processo
    localN = totalN / p;

    // Alloca spazio per le matrici globali e locali
    if (myRank == 0) {
        A = allocint2darray(totalN, totalN);
        B = allocint2darray(totalN, totalN);
    }
    localA = allocint2darray(localN, localN);
    localB = allocint2darray(localN, localN);
    localC = allocint2darray(localN, localN);
    broadcastA = allocint2darray(localN, localN);

    // Inizializza la matrice sul processo radice (rank 0)
    if (myRank == 0) {
        for (i = 0; i < totalN; i++) {
            for (j = 0; j < totalN; j++) {
                A[i][j] = i * totalN + j;
                B[i][j] = 1;
            }
        }
        // Stampa la matrice totale
        printLocalMatrix(A, totalN, totalN, myRank);
        printLocalMatrix(B, totalN, totalN, myRank);
    }

    for (i = 0; i < localN; i++) {
            for (j = 0; j < localN; j++) {
                localC[i][j] = 0;
            }
        }

    prodottoMatMat(A, B, localA, broadcastA, localB, localC, totalN, localN, numProcesses, p, myRank, coordinates, gridRow, gridCol);

    // Distribuisci la matrice tra i processi
    // if (myRank == 0)
    //     prodottoMatVet
    //     (&(matrix[0][0]), coordinates[0], coordinates[1], myRank, numProcesses, blocks, globalsizes, localsizes, localMatrix, gridRow, gridCol, xVec, xVecLoc, yVec, yVecLoc);
    // else
    //     prodottoMatVet
    //     (NULL, coordinates[0], coordinates[1], myRank, numProcesses, blocks, globalsizes, localsizes, localMatrix, gridRow, gridCol, xVec, xVecLoc, yVec, yVecLoc);

    printLocalMatrix(localC, localN, localN, myRank);

    // Deallocazione della memoria
    if (myRank == 0) {
        freeint2darray(A);
        freeint2darray(B);
    }
    freeint2darray(localA);
    freeint2darray(localB);
    freeint2darray(localC);
    freeint2darray(broadcastA);
    MPI_Comm_free(&grid);
    MPI_Comm_free(&gridRow);
    MPI_Comm_free(&gridCol);

    MPI_Finalize();
    return 0;

}

void prodottoMatMat(int** A, int** B, int** localA, int** broadcastA, int** localB, int** localC, int totalN, int localN, int numProcesses, int p, int myRank, int* coordinates, MPI_Comm gridRow, MPI_Comm gridCol) {

    int k;

    if (myRank == 0) {
        distributeMatrix(&A[0][0], coordinates[0], coordinates[1], myRank, numProcesses, {p, p}, {totalN, totalN}, {localN, localN}, localA, gridRow, gridCol); //Distribuzione della matrice A
        distributeMatrix(&B[0][0], coordinates[0], coordinates[1], myRank, numProcesses, {p, p}, {totalN, totalN}, {localN, localN}, localB, gridRow, gridCol); //Distribuzione della matrice B
    }
    distributeMatrix(NULL, coordinates[0], coordinates[1], myRank, numProcesses, {p, p}, {totalN, totalN}, {localN, localN}, localA, gridRow, gridCol); //Distribuzione della matrice A
    distributeMatrix(NULL, coordinates[0], coordinates[1], myRank, numProcesses, {p, p}, {totalN, totalN}, {localN, localN}, localB, gridRow, gridCol); //Distribuzione della matrice B

    for (k = 0; k < p; k++) { //Algoritmo principale che si ripete sulla k-esima diagonale

        //FASE DI BROADCAST

        if (coordinates[1] - coordinates[0] == k) { //I processi sulla k-esima diagonale principale dovranno inviare per primi il proprio blocco locale, quindi broadcastA coincide con localA
                                                    //I processi sulla diagonale k-esima avranno la differenza tra gli indici j ed i pari a k
            broadcastA = localA;
        }
        MPI_Bcast(broadcastA[0], localN*localN, MPI_INT, coordinates[0] + k, gridRow); 
        //A e' allocata in modo che A[0] punti ad un blocco di memoria contiguo che contiene l'intera matrice
        //Ad ogni riga, il processo sulla diagonale k-esima (con posizione row+k) inviera' agli altri processi

        //FASE DI MULTIPLY
        multiply(broadcastA, localB, localC, localN);

        //FASE DI ROLL
        roll(localB, localN, myRank, numProcesses, p, coordinates, gridCol);

    }

}

void multiply(int** broadcastA, int** localB, int** localC, int localN) {
    int i, j, k;
    for (i = 0; i < localN; i++) {
        for (j = 0; j < localN; j++) {
            for (k = 0; k < localN; k++) {
                localC[i][j] += broadcastA[i][k]*localB[k][j];
            }
        }
    }
}

void roll(int** localB, int localN, int myRank, int numProcesses, int p, int* coordinates, MPI_Comm gridCol) {
    int** recvB = allocint2darray(localN, localN);
    //INVIO ASINCRONO, RICEZIONE ASINCRONA, ATTESA SULLA RICEZIONE, MODIFICO localB
}