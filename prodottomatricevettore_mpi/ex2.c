#include "ex2.h"

// Funzione per distribuire equamente il vettore tra i processi
void distributeVector(int *V, int N, int *localV, int localN, int rank, int size) {
    int *sendcounts, *displs;
    int i;

    // Allocazione di array per determinare il numero di elementi da inviare a ciascun processo e gli spostamenti
    sendcounts = (int *)malloc(size * sizeof(int));
    displs = (int *)malloc(size * sizeof(int));

    // Calcola il numero medio di elementi per processo e gli elementi rimanenti
    int avg_elements = N / size;
    int remaining_elements = N % size;

    // Calcola il numero di elementi da inviare a ciascun processo e gli spostamenti
    for (i = 0; i < size; i++) {
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
