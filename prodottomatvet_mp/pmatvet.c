#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

// Se il compilatore non supporta la specifica del linguaggio C99 o successive
// allora definisci la macro 'restrict' come un commento vuoto.
#if __STDC_VERSION__ < 199901L
#define restrict /* nothing */
#endif


// Funzione che restituisce un vettore risultante dalla moltiplicazione matrice-vettore
// Parametri:
//   m, n: dimensioni della matrice A (n x m)
//   x: vettore
//   A: matrice
//   elapsedTime: parametro di output per il tempo impiegato
double *matxvet(int m, int n, double * restrict x, double * restrict A, double *elapsedTime);

// Funzione di utilità per stampare un vettore
// Parametri:
//   vec: il vettore da stampare
//   size: la dimensione del vettore
void printVec(double* vec, int size);

// Funzione per generare un numero casuale tra 0 e 1 (non utilizzata attualmente)
double randomDouble();

int main(int argc, char* argv[]) {
    srand(time(NULL));

    int n, m, i, j;
    double *A;
    double *x, *b;
    int iter;
    double deltat, avgt;

    // Controllo del numero corretto di argomenti da linea di comando
    if (argc != 3) {
        fprintf(stderr, "%s", "Numero di argomenti inseriti errato\n");
        return 1;
    }

    // Recupero delle dimensioni della matrice dalla linea di comando
    n = atoi(argv[1]);
    m = atoi(argv[2]);

    // Controllo sulle dimensioni della matrice
    if (n < 1 || m < 1) {
        fprintf(stderr, "%s", "Dimensione richiesta per la matrice non positiva\n");
        return 1;
    }

    // Controllo sul numero di thread disponibili rispetto alle righe della matrice
    if (n < omp_get_max_threads()) {
        fprintf(stderr, "%s", "Sono presenti meno righe nella matrice che thread a disposizione\n");
        return 1;
    }

    // Allocazione dinamica della memoria per la matrice A e il vettore x
    A = malloc(n * m * sizeof(double));
    x = malloc(m * sizeof(double));

    // Inizializzazione casuale della matrice A e del vettore x
    for (j = 0; j < m; j++) {
        for (i = 0; i < n; i++) {
            A[i * m + j] = randomDouble();
        }
        x[j] = randomDouble();
    }

    avgt = 0;

    // Esecuzione ripetuta per ottenere una media dei tempi di esecuzione
    for (iter = 0; iter < 10; iter++) {
        // Chiamata alla funzione matxvet per calcolare il prodotto matrice-vettore
        b = matxvet(m, n, x, A, &deltat);

        // Aggiornamento del tempo medio
        avgt += deltat;

        // printVec(b, n);

        // Liberazione della memoria allocata per il vettore risultante
        free(b);
    }

    // Calcolo della media dei tempi di esecuzione
    avgt /= 10.0;

    // Stampa dei risultati
    printf("Test con %d thread, taglia della matrice %dx%d, tempo impiegato medio %e\n", omp_get_max_threads(), n, m, avgt);

    // Liberazione della memoria allocata per il vettore x e la matrice A
    free(x);
    free(A);

    return 0;
}

// Implementazione della funzione matxvet
double *matxvet(int m, int n, double * restrict x, double * restrict A, double *elapsedTime) {
    int i, j;
    double *b;
    double t0, t1;
    struct timeval time;

    // Allocazione dinamica della memoria per il vettore risultante con inizializzazione a 0 dei valori iniziali
    b = calloc(n, sizeof(double));

    // Misurazione del tempo di inizio
    gettimeofday(&time, NULL);
    t0 = time.tv_sec + (time.tv_usec / 1000000.0);

    // Calcolo parallelo del prodotto matrice-vettore con OpenMP
    #pragma omp parallel for default(none) shared(m,n,A,x,b) private(i,j)
    for (i = 0; i < n; i++) { //La direttiva pragma omp for riguarda solo questo ciclo esterno, non il ciclo interno: i thread si divideranno i possibili valori di i
        for (j = 0; j < m; j++) //E poi cicleranno mediante j sull'intera i-esima riga della tabella: i thread si suddividono le intere righe
            b[i] += A[i * m + j] * x[j];
    }

    // Misurazione del tempo di fine
    gettimeofday(&time, NULL);
    t1 = time.tv_sec + (time.tv_usec / 1000000.0);

    // Calcolo del tempo trascorso
    *elapsedTime = t1 - t0;

    // Restituzione del vettore risultante
    return b;
}

// Implementazione della funzione printVec
void printVec(double* vec, int size) {
    int i;

    printf("Vettore: ");
    for (i = 0; i < size; i++) {
        printf("%.2f ", vec[i]);
    }
    printf("\n");
}

// Implementazione della funzione randomDouble (non utilizzata attualmente)
double randomDouble() {
    // double div = RAND_MAX / 100;
    // return (rand() / div);
    return 1;
}
