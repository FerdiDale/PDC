#ifndef EX3_H
#define EX3_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

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

#endif
