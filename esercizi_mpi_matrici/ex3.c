#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Funzione per creare una griglia di comunicazione MPI con righe e colonne
void createGrid(MPI_Comm *grid, MPI_Comm *gridRow, MPI_Comm *gridCol, int myRank, int numProcesses, int numRows, int numCols, int *coords);

// Funzione per distribuire una matrice tra i processi
void distributeMatrix(int *globalptr, const int myrow, const int mycol, const int rank, const int size, 
                      const int blocks[2], const int globalsizes[2], const int localsizes[2],
                      int **localdata, MPI_Comm rowComm, MPI_Comm colComm);

// Funzione per stampare la matrice locale di ogni processo
void printLocalMatrix(int **localMatrix, int localRows, int localCols, int myRank);

// Allocazione di una matrice bidimensionale di interi
int **allocint2darray(int n, int m);

// Liberazione della memoria allocata per una matrice bidimensionale di interi
void freeint2darray(int **a);

int main(int argc, char **argv) {
    int myRank, numProcesses;
    int **matrix, **localMatrix;
    int totalRows = 9; // Numero totale di righe nella matrice
    int totalCols = 7; // Numero totale di colonne nella matrice
    int localRows, localCols;
    MPI_Comm grid, gridRow, gridCol;
    int coordinates[2];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    // Verifica che il numero di processi formi una griglia pxq
    int p = 2; // Numero di processori lungo ogni riga della griglia
    int q = numProcesses / p; // Numero di processori lungo ogni colonna della griglia
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

    // Inizializza la matrice sul processo radice (rank 0)
    if (myRank == 0) {
        for (int i = 0; i < totalRows; i++) {
            for (int j = 0; j < totalCols; j++) {
                matrix[i][j] = i * totalCols + j;
            }
        }
        // Stampa la matrice locale su ogni processo
        printLocalMatrix(matrix, totalRows, totalCols, myRank);
    }

    int blocks[2] = {p, q};
    int globalsizes[2] = {totalRows, totalCols};
    int localsizes[2] = {localRows, localCols};
    
    // Distribuisci la matrice tra i processi
    distributeMatrix(&(matrix[0][0]), coordinates[0], coordinates[1], myRank, numProcesses, blocks, globalsizes, localsizes, localMatrix, gridRow, gridCol);

    // Stampa la matrice locale su ogni processo
    printLocalMatrix(localMatrix, localRows, localCols, myRank);

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

void createGrid(MPI_Comm *grid, MPI_Comm *gridRow, MPI_Comm *gridCol, int myRank, int numProcesses, int numRows, int numCols, int *coords) {
    int dim = 2, dims[2], reorder, periods[2], subCoords[2];

    // Imposta le dimensioni della griglia
    dims[0] = numRows;
    dims[1] = numCols;

    // Imposta i periodi per rendere la griglia non periodica su nessuna dimensione
    periods[0] = periods[1] = 0;

    // Non è necessario riordinare i ranghi dei processi
    reorder = 0;

    // Crea la griglia cartesiana
    MPI_Cart_create(MPI_COMM_WORLD, dim, dims, periods, reorder, grid);

    // Ottiene le coordinate del processo corrente nella griglia
    MPI_Cart_coords(*grid, myRank, 2, coords);

    // Imposta le coordinate per creare la griglia delle righe
    subCoords[0] = 0;
    subCoords[1] = 1;
    MPI_Cart_sub(*grid, subCoords, gridRow);

    // Imposta le coordinate per creare la griglia delle colonne
    subCoords[0] = 1;
    subCoords[1] = 0;
    MPI_Cart_sub(*grid, subCoords, gridCol);
}

void distributeMatrix(int *globalptr, const int myrow, const int mycol, const int rank, const int size,
                      const int blocks[2], const int globalsizes[2], const int localsizes[2],
                      int **localdata, MPI_Comm rowComm, MPI_Comm colComm) {
    int **rowdata = NULL;
    int rowBlocksize = globalsizes[0] / blocks[0];
    int colBlocksize = globalsizes[1] / blocks[1];
    int remainingRows = globalsizes[0] % blocks[0];
    int remainingCols = globalsizes[1] % blocks[1];

    /* Prima di tutto, distribuire l'array per righe, con il processore nella colonna 0 corrispondente a ogni riga
     * ricevendo i dati */
    if (mycol == 0) {
        int *sendcounts = (int *)malloc(blocks[0] * sizeof(int));
        int *senddispls = (int *)malloc(blocks[0] * sizeof(int));
        senddispls[0] = 0;

        for (int row = 0; row < blocks[0]; row++) {
            /* ogni processore ottiene blocksize righe, ognuna di dimensione globalsizes[1]... */
            sendcounts[row] = rowBlocksize;
            if (row < remainingRows) {
                sendcounts[row]++;
            }
            sendcounts[row] *= globalsizes[1];
            if (row > 0)
                senddispls[row] = senddispls[row - 1] + sendcounts[row - 1];
        }

        /* allocare i dati della mia riga */
        rowdata = allocint2darray(sendcounts[myrow], globalsizes[1]);

        /* eseguire la distribuzione delle righe */
        MPI_Scatterv(globalptr, sendcounts, senddispls, MPI_INT,
                     &(rowdata[0][0]), sendcounts[myrow], MPI_INT, 0, colComm);

        free(sendcounts);
        free(senddispls);
    }

    /* Ora, all'interno di ciascuna riga di processori, possiamo distribuire le colonne.
     * Possiamo farlo come nell'esempio precedente; creare un vettore
     * (e un vettore locale) e distribuire di conseguenza */
    int locnrows = rowBlocksize;
    if (myrow < remainingRows)
        locnrows++;

    MPI_Datatype vec, localvec;
    MPI_Type_vector(locnrows, 1, globalsizes[1], MPI_INT, &vec);
    MPI_Type_create_resized(vec, 0, sizeof(int), &vec);
    MPI_Type_commit(&vec);

    MPI_Type_vector(locnrows, 1, localsizes[1], MPI_INT, &localvec);
    MPI_Type_create_resized(localvec, 0, sizeof(int), &localvec);
    MPI_Type_commit(&localvec);

    int *sendcounts = (int *)malloc(blocks[1] * sizeof(int));
    int *senddispls = (int *)malloc(blocks[1] * sizeof(int));
    senddispls[0] = 0;
    for (int col = 0; col < blocks[1]; col++) {
        sendcounts[col] = colBlocksize;
        if (col < remainingCols) {
            sendcounts[col]++;
        }
        if (col > 0)
            senddispls[col] = senddispls[col - 1] + sendcounts[col - 1];
    }
    int *rowptr = (mycol == 0) ? &(rowdata[0][0]) : NULL;

    MPI_Scatterv(rowptr, sendcounts, senddispls, vec,
                 &(localdata[0][0]), sendcounts[mycol], localvec, 0, rowComm);

    MPI_Type_free(&localvec);
    MPI_Type_free(&vec);

    if (mycol == 0)
        freeint2darray(rowdata);

    free(sendcounts);
    free(senddispls);
}



void printLocalMatrix(int **localMatrix, int localRows, int localCols, int myRank) {
    printf("Process %d: Local Matrix:\n", myRank);
    for (int i = 0; i < localRows; i++) {
        for (int j = 0; j < localCols; j++) {
            printf("%4d", localMatrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int **allocint2darray(int n, int m) {
    int **ptrs = calloc(n,sizeof(int *));
    ptrs[0] = calloc(n*m,sizeof(int));
    for (int i=1; i<n; i++) 
        ptrs[i] = ptrs[i-1] + m;
    return ptrs;
}

void freeint2darray(int **a) {
    free(a[0]);
    free(a);
}

