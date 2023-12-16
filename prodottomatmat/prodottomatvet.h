#ifndef MATVET_H
#define MATVET_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "ex2.h"
#include "ex3.h"

void prodottoMatVet(int* globalptr, const int myrow, const int mycol, const int myRank, const int nproc,
    const int blocks[2], const int globalsizes[2], const int localsizes[2],
    int** localdata, MPI_Comm rowComm, MPI_Comm colComm, int* xVec, int* xVecLoc, int* yVec, int* yVecLoc);

void distributeGeneralVector(const int myrow, const int mycol, const int myRank,
    const int blocks[2], const int globalsizes[2], const int localsizes[2], MPI_Comm rowComm, MPI_Comm colComm, int* xVec, int* xVecLoc);

void printvec(int* vec, int size, int myRank);

void gatherVector(int* locVec, int* totalVec, MPI_Comm communicator, int nsend, int nProc, int N);

#endif
