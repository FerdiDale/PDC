#include "ex2.h"

// Funzione per distribuire equamente il vettore tra i processi
void distributeVector(int *V, int N, int *localV, int localN, int rank, int nProc, MPI_Comm communicator) {
    int *sendcounts, *displs;
    int i;

    // Allocazione di array per determinare il numero di elementi da inviare a ciascun processo e gli spostamenti
    sendcounts = (int *)malloc(nProc * sizeof(int));
    displs = (int *)malloc(nProc * sizeof(int));

    // Calcola il numero medio di elementi per processo e gli elementi rimanenti
    int avg_elements = N / nProc;
    int remaining_elements = N % nProc;

    // Calcola il numero di elementi da inviare a ciascun processo e gli spostamenti
    for (i = 0; i < nProc; i++) {
        sendcounts[i] = avg_elements;
        if (i < remaining_elements) {
            sendcounts[i]++;
        }

        displs[i] = (i > 0) ? (displs[i - 1] + sendcounts[i - 1]) : 0;
    }

    // printf("\nProcesso=%d : Prima dello scatter di vettore\n", rank);
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
    MPI_Scatterv(V, sendcounts, displs, MPI_INT, &localV[0], localN, MPI_INT, 0, communicator);

    // Libera la memoria allocata per gli array di invio e spostamento
    free(sendcounts);
    free(displs);
}
