#ifndef EX2_H
#define EX2_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void distributeVector(int *V, int N, int *localV, int localN, int rank, int size);

#endif
