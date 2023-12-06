#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#if __STDC_VERSION__ < 199901L
#define restrict /* nothing */
#endif

double *matxvet(int m, int n, double * restrict x, double * restrict A, double *elapsedTime);

void printVec(double* vec, int size);

double randomDouble();

int main(int argc, char* argv[]){

    srand(time(NULL));

    int n, m, i, j;
    double *A;
    double *x, *b;
    int iter;
    double deltat, avgt;

    if (argc != 3) { //Nome file, n, m
        fprintf(stderr, "%s", "Numero di argomenti inseriti errato\n");
        return 1;
    }

    n = atoi(argv[1]);
    m = atoi(argv[2]);

    if (n < 1 || m < 1) {
        fprintf(stderr, "%s", "Dimensione richiesta per la matrice non positiva\n");
        return 1;
    }

    if (n < omp_get_max_threads()) {
        fprintf(stderr, "%s", "Sono presenti meno righe nella matrice che thread a disposizione\n");
        return 1;
    }

    A = malloc(n*m*sizeof(double));
    x = malloc(m*sizeof(double));

    for (j = 0; j < m; j++) {
        for (i = 0; i < n; i++) {
            A[i*m+j] = randomDouble();
        }
        x[j] = randomDouble();
    }

    avgt = 0;

    for (iter = 0; iter < 10; iter++) { //Ripetiamo il processo 10 volte facendo la media dei tempi per avere un'accuratezza dei risultati maggiore

        b = matxvet(m, n, x, A, &deltat);

        avgt+=deltat;

    }

    printVec(b, n);

    avgt/=10.0;

    printf("Test con %d thread, taglia della matrice %dx%d, tempo impiegato medio %e\n", omp_get_max_threads(), n, m, avgt); 
    
    free(x);
    free(b);
    free(A);

    return 0;
}

double *matxvet(int m, int n, double * restrict x, double * restrict A, double *elapsedTime){
    int i,j;
    double *b;
    double t0,t1;
    struct timeval time;
    
    b = calloc(n,sizeof(double));

    gettimeofday(&time, NULL);
    t0=time.tv_sec+(time.tv_usec/1000000.0);

    #pragma omp parallel for default(none) shared(m,n,A,x,b) private(i,j)
    for (i=0; i<n; i++){ //La direttiva pragma omp for riguarda solo questo ciclo esterno, non il ciclo interno: i thread si divideranno i possibili valori di i
        for (j=0; j<m; j++) //E poi cicleranno mediante j sull'intera i-esima riga della tabella: i thread si suddividono le intere righe
            b[i] += A[i*m+j]*x[j];
    }

    gettimeofday(&time, NULL);
    t1=time.tv_sec+(time.tv_usec/1000000.0);

    *elapsedTime = t1-t0;

    return b;

}

void printVec(double* vec, int size) {
    int i;
    
    printf("Vettore: ");
    for (i = 0; i < size; i++) {
        printf("%.2f ", vec[i]);
    }
    printf("\n");
}

double randomDouble() {
    double div = RAND_MAX / 100;
    return (rand() / div);
    // return 1;
}
