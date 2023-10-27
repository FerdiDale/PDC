#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]){

    if(argc != 2){
        perror("Numero di parametri sbagliati!");
        exit(1);
    }

    int n = atoi(argv[1]);

    if(n <= 0){
        perror("Numero di parametri sbagliati!");
        exit(1);
    }

    int i = 0;
    int sommatot = 0;
    int* numeri = malloc(n * sizeof(int*));

    for(i=0; i<n; i++){
        numeri[i] = 1;
    }

    #pragma omp parallel 
    {
        #pragma omp for reduction (+:sommatot), private(i)
        for (i = 0; i < n; i++){
            sommatot += numeri[i];
        }

        // Stampo informazioni thread
        printf("Thread %d, somma parziale: %d\n", omp_get_thread_num(), sommatot);

    }


    printf("Il totale e' : %d\n", sommatot);
        

    return 0;
}
