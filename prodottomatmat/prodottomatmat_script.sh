#!/bin/bash

# Specifica i valori di NCPU e DIM
NCPU_VALUES="1 4"
DIM_VALUES="500 1000 1500 2000"
NUM_CYCLES=10

# Ciclo su valori di NCPU
for NCPU in $NCPU_VALUES; do
    # Ciclo su valori di DIM
    for DIM in $DIM_VALUES; do
        # Ciclo interno
        for ((CICLO=1; CICLO<=NUM_CYCLES; CICLO++)); do
            # Chiamare lo script PBS con i valori correnti di NCPU e DIM
            qsub -v NCPU=$NCPU,DIM=$DIM,CICLO=$CICLO prodottomatmat_in_script.pbs
        done
    done
done
