#include <stdio.h>
#include “mpi.h”
main(int argc, char *argv[]){
    int menum, nproc,... ;
    int n, nloc, tag, i, strat, start_pos, sum, parz_sum;
    int *x, *xloc;
    MPI_Status status;

    MPI_Init(&argv, &argc);

    MPI_Comm_rank(MPI_COMM_WORLD, &menum);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    //Acquisizione dell'input
    if (menum==0){
        if (argc < 3) {
            perror("Numero di argomenti insufficiente");
            return 1;
        }
        n = atoi(argv[1]);   
        x* = (int*) malloc(n*(sizeof(int)));
        if (n <= 20) {
            for (i = 0 ; i < n; i++)
                x[i] = atoi(argv[2+i]);
        } 
        else {
            for (i = 0; i < n; i++)
                x[i] = 1;
        }

        strat = atoi(argv[n+2]); 
    }

    //Comunicazione dei dati
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&strat, 1, MPI_INT, 0, MPI_COMM_WORLD);

    nloc=n/nproc
    rest=n%nproc
    if (menum<rest) nloc=nloc+1

    if (menum == 0) {
        xloc = x;
        tmp = nloc;
        start_pos = 0;
        for (i = 1; i < nproc ; i++) {
            start_pos += tmp;
            if (i == rest) 
                tmp--;
            tag = 22+i;
            MPI_Send(&x[start_pos], tmp, MPI_INT, i, tag, MPI_COMM_WORLD);
        }
    }
    else {
        xloc = (int*)malloc(nloc*sizeof(int));
        tag = 22 + menum;
        MPI_Recv(xloc, nloc, MPI_INT, 0, tag, MPI_COMM_WORLD);
    }

    //Calcolo locale
    sum = xloc[0];
    for (i = 1; i < nloc; i++) 
        sum += xloc[i];


    //Comunicazione
    if (strat == 1) { //I Strategia
        if (menum == 0) {
            for (i = 1; i < nproc; i++) {
                tag = 80;
                MPI_Recv(&parz_sum, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
                sum+=parz_sum; 
            }
        }
        else {
            tag = 80+menum;
            MPI_Send(&sum, 1, MPI_INT, 0, tag, MPI_COMM_WORLD)
        }
    }

}