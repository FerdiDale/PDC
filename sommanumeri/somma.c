#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"

int random_int();

int main(int argc, char *argv[]) {
    int menum, nproc, tag, i, start_pos, rest, strat;
    long n, nloc, sum, parz_sum, tmp;
    int *x, *xloc;
    int iter;
    double log_nproc;
    double start_time, end_time, elapsed_time, total_time;
    double time_mean = 0.0;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &menum);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    for (iter = 0; iter < 10; iter++) { //Eseguiamo tutto 10 volte per avere una media dei tempi 

        //Acquisizione dell'input
        if (menum==0){
            if (argc < 3) {
                perror("Numero di argomenti insufficiente");
                return 1;
            }
            n = atoi(argv[1]);   
            x = (int*) malloc(n*(sizeof(int)));
            if (n <= 20) {
                if (argc != n+3) {
                    perror("Numero di argomenti errato!");
                    return 1;
                }
                for (i = 0 ; i < n; i++)
                    x[i] = atoi(argv[2+i]);


                strat = atoi(argv[n+2]); 
            } 
            else {
                for (i = 0; i < n; i++)
                    x[i] = random_int();
                
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
        MPI_Bcast(&n, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&strat, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //Eventuale cambio di strategia
        if (strat == 2 || strat == 3) { 
            if ((nproc & (nproc-1)) != 0) { //numero di processi non potenza di 2, sfruttiamo il codice binario particolare delle potenze di 2, composto da tutti 0 eccetto per un 1
            //Se nroc è potenza di 2, nproc-1 sarà quindi composto da tutti 1 fino ad uno 0 laddove nproc aveva il suo unico 1. L'and bit a bit quindi ci restituirà 0.
                if (menum == 0) 
                    printf("E' stata richiesta la strategia %d ma il numero di processi richiesta non e' potenza di due, si prosegue con la strategia 1\n", strat);
                strat = 1;
            } 
        }

        //Comunicazione dei dati

        nloc=n/nproc; //Dividiamo il numero di input per il numero di processori per distribuirli equamente, dobbiamo però stare attenti al resto, come vediamo in seguito
        rest=n%nproc;
        if (menum<rest) nloc=nloc+1; //Fino al rest_esimo processore dobbiamo dare un numero in input in più

        if (menum == 0) {
            xloc = x;
            tmp = nloc;
            start_pos = 0;
            for (i = 1; i < nproc ; i++) {
                start_pos += tmp;
                if (i == rest) //Arrivati al rest_esimo processore rimuoviamo il numero extra
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

        //Calcolo di potenze e logaritmi, per efficienza negli algoritmi delle strategie seguenti
        int pow [nproc];
        int curr_pow = 1;

        for (i = 0; i < nproc; i++) {
            pow[i] = curr_pow;
            curr_pow*=2;
        }

        log_nproc = log(nproc)/log(2);

        MPI_Barrier(MPI_COMM_WORLD);
        start_time = MPI_Wtime(); //Prendiamo il tempo di inizio, dopo aver sincronizzato i vari processi con la precedente Barrier

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
                        MPI_Recv(&parz_sum, 1, MPI_LONG, i, tag, MPI_COMM_WORLD, &status);
                        sum+=parz_sum; 
                    }
                    end_time = MPI_Wtime();
                    printf("La somma totale calcolata con la prima strategia e' %ld\n", sum);
                }
                else {
                    tag = 80+menum;
                    MPI_Send(&sum, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                    end_time = MPI_Wtime(); //Otteniamo il tempo di fine
                }
                break;
            
            case 2: //II strategia
                for (i = 0; i < log_nproc; i++) { //Avremo log2(nproc) passi, dato che ad ogni passo accumuliamo la somma di 2^i processi
                    if (menum % pow[i] == 0) { //All'i-esimo passo parteciperanno solo i processi il cui rank è divisibile per 2^i
                        if (menum % pow[i+1] == 0) {//I processi che riceveranno i dati sono quelli che parteciperanno alla prossima iterazione, cioè con rank divisibile per 2^(i+1)
                            tag = 120 + menum + pow[i]; //menum + 2^i ci darà il rank del processo accoppiato che invierà i dati
                            MPI_Recv(&parz_sum, 1, MPI_LONG, menum+pow[i], tag, MPI_COMM_WORLD, &status);
                            sum+=parz_sum; 
                        }
                        else {
                            tag = 120 + menum;
                            MPI_Send(&sum, 1, MPI_LONG, menum-pow[i], tag, MPI_COMM_WORLD);
                            break;
                        }
                    }
                }
                end_time = MPI_Wtime();//Otteniamo il tempo di fine
                if (menum == 0)
                    printf("La somma totale calcolata con la seconda strategia e' %ld\n", sum);
                break;
            
            case 3: //III strategia
                for (i = 0; i < log_nproc; i++) {
                    if (menum % pow[i+1] < pow[i]) { //Distinguiamo i processi per formare le coppie, ogni processo comunicherà con quello distante 2^i
                        tag = 200 + menum + pow[i]; //Invia prima il secondo processo, con rank menum+2^i per il ricevente
                        MPI_Recv(&parz_sum, 1, MPI_LONG, menum+pow[i], tag, MPI_COMM_WORLD, &status);
                        tag = 200 + menum;
                        MPI_Send(&sum, 1, MPI_LONG, menum+pow[i], tag, MPI_COMM_WORLD);
                        sum+=parz_sum; 
                    }
                    else {
                        tag = 200 + menum;
                        MPI_Send(&sum, 1, MPI_LONG, menum-pow[i], tag, MPI_COMM_WORLD);
                        tag = 200 + menum - pow[i];//Invia adesso il secondo processo, con rank menum-2^i per il ricevente
                        MPI_Recv(&parz_sum, 1, MPI_LONG, menum-pow[i], tag, MPI_COMM_WORLD, &status);
                        sum+=parz_sum; 
                    }
                }
                end_time = MPI_Wtime(); //Otteniamo il tempo di fine
                printf("Sono il processo %d, la somma totale calcolata con la terza strategia e' %ld\n", menum, sum);
                break;
            
            default:
                if (menum == 0) 
                    perror("Input della strategia errato");
                return 1;

        }

        if (menum != 0)
            free(xloc);
        else 
            free(x);

        elapsed_time = end_time - start_time;
        MPI_Reduce(&elapsed_time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); //Prendiamo il massimo dei tempi calcolati dai vari processi

        if (menum == 0) {
            time_mean+=total_time;
        }
        
    }
        
    if (menum == 0) {  
        time_mean/=10; //Facciamo la media dei risultati delle 10 iterazioni del processo
        printf("Strategia usata: %d,\nNumero di input: %ld,\nNumero di processori:%d,\nTempo impiegato: %e\n\n\n\n", strat, n, nproc, time_mean);
    }

    MPI_Finalize();

    return 0;

}

int random_int() {
    srand(time(NULL));
    return rand()%100;
}
