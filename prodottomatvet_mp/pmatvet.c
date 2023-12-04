#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#if __STDC_VERSION__ < 199901L
#define restrict /* nothing */
#endif

double *matxvet(int m, int n, double * restrict x, double * restrict A, double *elapsedTime);

double **allocint2darray(int n, int m);

void freeint2darray(double **a);

void printVec(double* vec, int size);

int get_num_threads();

int main(int argc, char* argv[]){

    int n, m, i, j;
    double *A;
    double *x, *b;
    int iter;
    double deltat, avgt;

    if (argc != 3) { //Nome file, n, m
        perror("Numero di argomenti inseriti errato\n");
        return 1;
    }

    n = atoi(argv[1]);
    m = atoi(argv[2]);

    A = malloc(n*m*sizeof(double));
    x = malloc(m*sizeof(double));

    for (j = 0; j < m; j++) {
        for (i = 0; i < n; i++) {
            A[i*m+j] = 1;
        }
        x[j] = 1;
    }

    avgt = 0;

    for (iter = 0; iter < 10; iter++) { //Ripetiamo il processo 10 volte facendo la media dei tempi per avere un'accuratezza dei risultati maggiore

        b = matxvet(m, n, x, A, &deltat);

        avgt+=deltat;

        // printVec(b, n);

    }

    avgt/=10.0;

    printf("Test con %d thread, taglia della matrice %dx%d, tempo impiegato medio %e\n", get_num_threads(), n, m, avgt); 
    
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
    for (i=0; i<n; i++){

        for (j=0; j<m; j++)
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

int get_num_threads() {
    int num_threads;
    #pragma omp parallel
    {
        #pragma omp single
        num_threads = omp_get_num_threads();
    }
    return num_threads;
}
