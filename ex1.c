#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

int main(int argc, char **argv) {
    int menum, nproc, P, p;
    int px, q;
    int dim, *ndim, reorder, *period, *coordinate;
    MPI_Comm comm_grid;
    int menum_grid;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &menum);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    if (argc != 3) {
        if (menum == 0) {
            printf("Usage: %s P p\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    P = atoi(argv[1]);
    p = atoi(argv[2]);

    if (P != nproc || P % p != 0) {
        if (menum == 0) {
            printf("Il numero di processi specificato (%d) non è corretto per P (%d) o p (%d)\n", nproc, P, p);
        }
        MPI_Finalize();
        return 1;
    }

    q = P / p;

    dim = 2;
    coordinate = (int *)calloc(dim, sizeof(int));

    ndim = (int *)calloc(dim, sizeof(int));
    ndim[0] = p;
    ndim[1] = q;

    period = (int *)calloc(dim, sizeof(int));
    period[0] = period[1] = 0;
    reorder = 0;

    MPI_Cart_create(MPI_COMM_WORLD, dim, ndim, period, reorder, &comm_grid);
    MPI_Comm_rank(comm_grid, &menum_grid);

    MPI_Cart_coords(comm_grid, menum_grid, dim, coordinate);

    printf("Processo %d: Coordinata nella griglia: (%d, %d)\n", menum, coordinate[0], coordinate[1]);

    MPI_Finalize();
    return 0;
}
