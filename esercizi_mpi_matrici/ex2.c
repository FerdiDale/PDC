#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

// Funzione per distribuire equamente il vettore tra i processi
void distributeVector(int *V, int N, int *localV, int localN, int rank, int size) {
    int *sendcounts, *displs;

    // Allocazione di array per determinare il numero di elementi da inviare a ciascun processo e gli spostamenti
    sendcounts = (int *)malloc(size * sizeof(int));
    displs = (int *)malloc(size * sizeof(int));

    // Calcola il numero medio di elementi per processo e gli elementi rimanenti
    int avg_elements = N / size;
    int remaining_elements = N % size;

    // Calcola il numero di elementi da inviare a ciascun processo e gli spostamenti
    for (int i = 0; i < size; i++) {
        sendcounts[i] = avg_elements;
        if (i < remaining_elements) {
            sendcounts[i]++;
        }

        displs[i] = (i > 0) ? (displs[i - 1] + sendcounts[i - 1]) : 0;
    }

    // Distribuisci i dati tra i processi
    /* MPI_Scatterv:
    * V è il vettore globale che si desidera distribuire tra i processi.
    * sendcounts specifica il numero di elementi inviati a ciascun processo.
    * displs specifica gli spostamenti degli elementi da V per ciascun processo.
    * MPI_INT è il tipo di dato degli elementi in V e localV.
    * localV è il vettore locale che riceverà i dati.
    * localN è il numero di elementi ricevuti da ciascun processo.
    * 0 è il rango del processo sorgente (in questo caso è il processo root).
    * MPI_COMM_WORLD è il communicator dei processi.
    */
    MPI_Scatterv(V, sendcounts, displs, MPI_INT, localV, localN, MPI_INT, 0, MPI_COMM_WORLD);

    // Libera la memoria allocata per gli array di invio e spostamento
    free(sendcounts);
    free(displs);
}

int main(int argc, char **argv) {
    int rank, size;
    int N, P; // Dimensione del vettore
    int *V, *localV;
    int localN;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Verifica che siano forniti gli argomenti corretti
    if (argc != 3) {
        if (rank == 0) {
            printf("Errore: usage: %s N P\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    // Ottieni N e P dagli argomenti di input
    P = atoi(argv[1]);
    N = atoi(argv[2]);

    // Verifica che il numero di processi P sia valido
    if (P != size) {
        if (rank == 0) {
            printf("Errore: il numero di processi non coincide con la dimensione del communicator\n");
        }
        MPI_Finalize();
        return 1;
    }

    // Calcola la dimensione del vettore locale per ciascun processo
    localN = N / P;
    if (rank < N % P) {
        localN++;
    }

    // Allocazione di spazio per il vettore globale e il vettore locale
    V = (int *)malloc(N * sizeof(int));
    localV = (int *)malloc(localN * sizeof(int));

    // Inizializza il vettore V sul processo root (rank 0)
    if (rank == 0) {
        for (int i = 0; i < N; i++) {
            V[i] = i;
        }
    }

    // Distribuisci il vettore V tra i processi
    distributeVector(V, N, localV, localN, rank, size);

    // Stampa il vettore locale su ogni processo
    printf("Processo %d: Vettore locale: ", rank);
    printf("[");
    for (int i = 0; i < localN; i++) {
        printf("%d ", localV[i]);
    }
    printf("]");
    printf("\n");

    // Deallocazione della memoria
    free(V);
    free(localV);

    MPI_Finalize();
    return 0;
}
