#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

main(int argc, char *argv[]){
    int menum, nproc;
    int n, nloc, tag, i, strat, start_pos, sum, parz_sum, rest, tmp;
    int *x, *xloc;
    double log_nproc;
    double start_time, end_time, elapsed_time, total_time;
    double time_mean = 0.0;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &menum);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    for (int iter = 0; iter++; iter < 10) { //Eseguiamo tutto 10 volte per avere una media dei tempi 

        //Acquisizione dell'input
        if (menum==0){
            if (argc < 3) {
                perror("Numero di argomenti insufficiente");
                return 1;
            }
            n = atoi(argv[1]);   
            x = (int*) malloc(n*(sizeof(int)));
            if (n <= 20) {
                for (i = 0 ; i < n; i++)
                    x[i] = atoi(argv[2+i]);


                strat = atoi(argv[n+2]); 
            } 
            else {
                for (i = 0; i < n; i++)
                    x[i] = 1;
                
                strat = atoi(argv[2]);
            }

            if (n <= 0) {
                perror("Numero di input non valido!");
                return 1;
            }

            if (n < nproc) {
                perror("Numero di input minore del numero di processori utilizzati");
                return 1;
            }
        }

        //Comunicazione dei dati
        MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&strat, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //Eventuale cambio di strategia
        if (strat == 2 || strat == 3) { 
            if ((nproc & (nproc-1)) != 0) { //numero di processi non potenza di 2
                printf("E' stata richiesta la strategia %d ma il numero di processi richiesta non e' potenza di due, si prosegue con la strategia 1\n", strat);
                strat = 1;
            } 
        }

        //Comunicazione dei dati

        nloc=n/nproc;
        rest=n%nproc;
        if (menum<rest) nloc=nloc+1;

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
            MPI_Recv(xloc, nloc, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
        }

        //Calcolo di potenze e logaritmi
        int pow [nproc];
        int curr_pow = 1;

        for (i = 0; i < nproc; i++) {
            pow[i] = curr_pow;
            curr_pow*=2;
        }

        log_nproc = log(nproc)/log(2);

        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime();

        //Calcolo locale
        sum = xloc[0];
        for (i = 1; i < nloc; i++) 
            sum += xloc[i];

        //Comunicazione
        switch(strat) { 
            case 1: //I strategia
                if (menum == 0) {
                    for (i = 1; i < nproc; i++) {
                        tag = 80+i;
                        MPI_Recv(&parz_sum, 1, MPI_INT, i, tag, MPI_COMM_WORLD, &status);
                        sum+=parz_sum; 
                    }
                    end_time = MPI_Wtime();
                    printf("La somma totale calcolata con la prima strategia e' %d\n", sum);
                }
                else {
                    tag = 80+menum;
                    MPI_Send(&sum, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
                    end_time = MPI_Wtime();
                }
                break;
            
            case 2: //II strategia
                for (i = 0; i < log_nproc; i++) {
                    if (menum % pow[i] == 0) {
                        if (menum % pow[i+1] == 0) {
                            tag = 120 + menum + pow[i];
                            MPI_Recv(&parz_sum, 1, MPI_INT, menum+pow[i], tag, MPI_COMM_WORLD, &status);
                            sum+=parz_sum; 
                        }
                        else {
                            tag = 120 + menum;
                            MPI_Send(&sum, 1, MPI_INT, menum-pow[i], tag, MPI_COMM_WORLD);
                            break;
                        }
                    }
                }
                end_time = MPI_Wtime();
                if (menum == 0)
                    printf("La somma totale calcolata con la seconda strategia e' %d\n", sum);
                break;
            
            case 3: //III strategia
                for (i = 0; i < log_nproc; i++) {
                    if (menum % pow[i+1] < pow[i]) {
                        tag = 200 + menum + pow[i];
                        MPI_Recv(&parz_sum, 1, MPI_INT, menum+pow[i], tag, MPI_COMM_WORLD, &status);
                        tag = 200 + menum;
                        MPI_Send(&sum, 1, MPI_INT, menum+pow[i], tag, MPI_COMM_WORLD);
                        sum+=parz_sum; 
                    }
                    else {
                        tag = 200 + menum;
                        MPI_Send(&sum, 1, MPI_INT, menum-pow[i], tag, MPI_COMM_WORLD);
                        tag = 200 + menum - pow[i];
                        MPI_Recv(&parz_sum, 1, MPI_INT, menum-pow[i], tag, MPI_COMM_WORLD, &status);
                        sum+=parz_sum; 
                    }
                }
                end_time = MPI_Wtime();
                printf("Sono il processo %d, la somma totale e' %d\n", menum, sum);
                break;
            
            default:
                perror("Input della strategia errato");
                return 1;

        }

        elapsed_time = end_time - start_time;
        MPI_Reduce(elapsed_time, total_time, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

        if (menum == 0)
            time_mean+=total_time;
        
    }
        
    if (menum == 0) {  
        time_mean/=10;
        printf("Tempo impiegato: %e", time_mean);
    }

    MPI_Finalize();

    return 0;

}
