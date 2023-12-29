#!/bin/bash

# Specifica dei valori per il numero di processori (NCPU), le dimensioni della matrice (DIM) e il numero di cicli (NUM_CYCLES)
NCPU_VALUES="1 4"
DIM_VALUES="500 1000 1500 2000"
NUM_CYCLES=10

# Itera su i valori di NCPU
for NCPU in $NCPU_VALUES; do
    # Itera su i valori di DIM
    for DIM in $DIM_VALUES; do
        # Ciclo interno
        for ((CICLO=1; CICLO<=NUM_CYCLES; CICLO++)); do
            # Chiamare lo script PBS con i valori correnti di NCPU, DIM e CICLO
            qsub -v NCPU=$NCPU,DIM=$DIM,CICLO=$CICLO prodottomatmat_in_script.pbs
        done
    done
done

