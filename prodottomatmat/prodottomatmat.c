#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "prodottomatvet.h"
#include "ex2.h"
#include "ex3.h"

void prodottoMatMat(int** A, int** B, int** localA, int** broadcastA, int** localB, int** localC, int totalN, int localN, int numProcesses, int p, int myRank, int* coordinates, MPI_Comm gridRow, MPI_Comm gridCol);

int isPerfectSquare(int num);

void multiply(int** broadcastA, int** localB, int** localC, int localN);

void roll(int** localB, int localN, int myRank, int numProcesses, int p, int* coordinates, MPI_Comm gridCol);

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
            fprintf(stderr, "%s", "Errore: numero di argomenti inseriti errato\n");
        MPI_Finalize();
        return 1;
    }

    totalN = atoi(argv[1]);

    if (totalN < 1) { 
        if (myRank == 0)
            fprintf(stderr, "%s", "Errore: la dimensione della matrice è non positiva\n");
        MPI_Finalize();
        return 1;
    }

    if (!isPerfectSquare(numProcesses)) { //numero di processi non potenza di 2, sfruttiamo il codice binario particolare delle potenze di 2, composto da tutti 0 eccetto per un 1
            //Se nroc è potenza di 2, nproc-1 sarà quindi composto da tutti 1 fino ad uno 0 laddove nproc aveva il suo unico 1. L'and bit a bit quindi ci restituirà 0.
        if (myRank == 0)
            fprintf(stderr, "%s", "Errore: il numero di processi deve essere un quadrato perfetto\n");
        MPI_Finalize();
        return 1;
    }

    p = sqrt(numProcesses);

    if (totalN % p != 0) {
        if (myRank == 0)
            fprintf(stderr, "%s", "Errore: la taglia della matrice quadrata deve essere multiplo della taglia della griglia di processi\n");
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

        // Verifica se le allocazioni sono riuscite
        if (A == NULL || B == NULL) {
            fprintf(stderr, "Errore: Allocazione di memoria per A o B non riuscita.\n");
            MPI_Finalize();
            return 1;
        }
    }

    localA = allocint2darray(localN, localN);
    localB = allocint2darray(localN, localN);
    localC = allocint2darray(localN, localN);
    broadcastA = allocint2darray(localN, localN);

    // Verifica se le allocazioni sono riuscite
    if (localA == NULL || localB == NULL || localC == NULL || broadcastA == NULL) {
        fprintf(stderr, "Errore: Allocazione di memoria per matrici locali non riuscita.\n");
        MPI_Finalize();
        return 1;
    }

    // Inizializza la matrice sul processo radice (rank 0)
    if (myRank == 0) {
        for (i = 0; i < totalN; i++) {
            for (j = 0; j < totalN; j++) {
                A[i][j] = i * totalN + j;
                B[i][j] = 1;
            }
        }
        // Stampa la matrice totale
        // printLocalMatrix(A, totalN, totalN, myRank);
        // printLocalMatrix(B, totalN, totalN, myRank);
    }

    for (i = 0; i < localN; i++) {
        for (j = 0; j < localN; j++) {
            localC[i][j] = 0;
        }
    }

    prodottoMatMat(A, B, localA, broadcastA, localB, localC, totalN, localN, numProcesses, p, myRank, coordinates, gridRow, gridCol);
    
    if (totalN <= 20) {
        for(i=0;i<numProcesses;i++){
            if(myRank==i){
                printLocalMatrix(localC, localN, localN, myRank);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

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

    int k, iter;
    int gridDims [2] = {p, p};
    int totalMatrixDims [2] = {totalN, totalN};
    int localMatrixDims [2] = {localN, localN};
    double start_time, end_time, elapsed_time, total_time;

    if (myRank == 0) {
        distributeMatrix(&A[0][0], coordinates[0], coordinates[1], myRank, numProcesses, gridDims, totalMatrixDims, localMatrixDims, localA, gridRow, gridCol); //Distribuzione della matrice A
        distributeMatrix(&B[0][0], coordinates[0], coordinates[1], myRank, numProcesses, gridDims, totalMatrixDims, localMatrixDims, localB, gridRow, gridCol); //Distribuzione della matrice B
    } else {
        distributeMatrix(NULL, coordinates[0], coordinates[1], myRank, numProcesses, gridDims, totalMatrixDims, localMatrixDims, localA, gridRow, gridCol); //Distribuzione della matrice A
        distributeMatrix(NULL, coordinates[0], coordinates[1], myRank, numProcesses, gridDims, totalMatrixDims, localMatrixDims, localB, gridRow, gridCol); //Distribuzione della matrice B
    }

    for (iter = 0; iter < 10; iter++) { //Ripetiamo l'intero processo 10 volte per avere una misura del tempo di esecuzione piu' accurata

        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime(); //Prendiamo il tempo di inizio, dopo aver sincronizzato i vari processi con la precedente Barrier

        for (k = 0; k < p; k++) { //Algoritmo principale che si ripete sulla k-esima diagonale

            //FASE DI BROADCAST

            if ((coordinates[0]+k)%p == coordinates[1]) { //I processi sulla k-esima diagonale principale dovranno inviare per primi il proprio blocco locale, quindi broadcastA coincide con localA
                                                        //I processi sulla diagonale k-esima avranno la differenza tra gli indici j ed i pari a k
                broadcastA = localA;
            }
            MPI_Bcast(broadcastA[0], localN*localN, MPI_INT, (coordinates[0] + k)%p, gridRow); 
            //A e' allocata in modo che A[0] punti ad un blocco di memoria contiguo che contiene l'intera matrice
            //Ad ogni riga, il processo sulla diagonale k-esima (con posizione row+k mod p) inviera' agli altri processi

            //FASE DI MULTIPLY
            multiply(broadcastA, localB, localC, localN);

            if (k == p-1)
                break; //L'ultimo roll non è necessario, risparmiamo tempo

            //FASE DI ROLL
            roll(localB, localN, myRank, numProcesses, p, coordinates, gridCol);

        }


        end_time = MPI_Wtime();//Otteniamo il tempo di fine

        elapsed_time = end_time - start_time;
        MPI_Reduce(&elapsed_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); //Prendiamo il massimo dei tempi calcolati dai vari processi

        if (myRank == 0)
            printf("Dimensione matrice: %dx%d,\nNumero di processi: %d,\nTempo impiegato: %e\n\n\n\n", totalN, totalN, numProcesses, total_time);

    }

}

int isPerfectSquare(int num) {
    // Calcola la radice quadrata
    double squareRoot = sqrt((double)num);

    // Verifica se la radice quadrata è un numero intero
    return squareRoot == floor(squareRoot);
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
    MPI_Request sendRequest;
    int i, colRank;
    int** recvB = allocint2darray(localN, localN);

    MPI_Comm_rank(gridCol, &colRank);

    // Invio la matrice locale alla riga inferiore
    MPI_Isend(localB[0], localN * localN, MPI_INT, (colRank + p - 1) % p, 30 + colRank, gridCol, &sendRequest);
    
    // Ricevo la matrice locale dalla riga inferiore
    MPI_Recv(recvB[0], localN * localN, MPI_INT, (colRank + 1) % p, 30 + ((colRank + 1) % p), gridCol, MPI_STATUS_IGNORE);
    
    // Attendo la fine dell'invio
    MPI_Wait(&sendRequest, MPI_STATUS_IGNORE);

    // Copio i dati dalla matrice ricevuta a localB
    memcpy(localB[0], recvB[0], localN * sizeof(int));

    // Liberare la memoria allocata per recvB
    free(recvB[0]);
    free(recvB);
}

