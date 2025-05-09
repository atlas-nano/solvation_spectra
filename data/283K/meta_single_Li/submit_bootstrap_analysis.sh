#!/bin/bash

DESTINATION_BASE="/global/cfs/cdirs/m4248/xiaoxusr/solvation_283K/meta_single_Li/"
conda activate ele_machine_clone
# LI_INDEXES=(2946 3012 3000 3012 2148 1500 1260 1104 1176 789 708 636) # id-1
# CONCENTRATIONS=(0.5 1 2 5 7 10 12 15 17 19 21 23)
LI_INDEXES=(3012 3000 3012 2148 1500 1260 1104 1176 789 708 636) # id-1
CONCENTRATIONS=(1 2 5 7 10 12 15 17 19 21 23)
NUMBER_OF_CVS=(3)

for i in "${!CONCENTRATIONS[@]}"; do
    CONC=${CONCENTRATIONS[$i]}
    for j in "${!NUMBER_OF_CVS[@]}"; do  
        cd ${DESTINATION_BASE}/${CONC}M/meta
        echo "We are in directory ${DESTINATION_BASE}/${CONC}M/meta"
        sbatch analysis_bootstrap2.sh
        echo "Analysis for ${CONC}M meta submitted!"
    done
done