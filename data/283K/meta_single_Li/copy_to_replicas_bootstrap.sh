#!/bin/bash
LI_INDEXES=(2946 3012 3000 3012 2148 1500 1260 1104 1176 789 708 636) # id-1
CONCENTRATIONS=(0.5 1 2 5 7 10 12 15 17 19 21 23)
NUMBER_OF_CVS=(3)
TEMP=283


ANALYSIS_SOURCE="/global/cfs/cdirs/m4248/xiaoxusr/solvation_${TEMP}K/meta_single_Li/analysis"
DESTINATION_BASE="/global/cfs/cdirs/m4248/xiaoxusr/solvation_${TEMP}K/meta_single_Li"


for i in "${!CONCENTRATIONS[@]}"; do
    CONC=${CONCENTRATIONS[$i]}
    LI_INDEX=${LI_INDEXES[$i]} # (id-1)!!!
    for j in "${!NUMBER_OF_CVS[@]}"; do  
        NUMBER_OF_CV=${NUMBER_OF_CVS[$j]}     
        cp -r ${ANALYSIS_SOURCE}/*.py ${DESTINATION_BASE}/${CONC}M/meta/
        cat > "${DESTINATION_BASE}/${CONC}M/meta/analysis_bootstrap.sh" <<EOL
#!/bin/bash
#SBATCH -C cpu
#SBATCH -t 12:00:00
#SBATCH -J free_energy_analysis
#SBATCH -o free_energy_analysis.o%j
#SBATCH -A m4248
#SBATCH -N 1
#SBATCH --ntasks-per-node=128
#SBATCH -q regular
# variables
LI_INDEX=${LI_INDEX} # (id-1)!!!
NUMBER_OF_CV=${NUMBER_OF_CV}  # Customize as needed
CONC=${CONC}
BASE_PATH="${DESTINATION_BASE}/${CONC}M/meta"
NSTRIDES=10
O_RADII=2.65
H_RADII=2.95
CL_RADII=3.05
TEMP=${TEMP}

# Run Python script
conda activate ele_machine_clone
python colvars_analyzer_script2.py --base_dir \$BASE_PATH --number_of_cv \$NUMBER_OF_CV &> colvars_analysis.log
echo "Colvars analysis completed!"
RESULTS_DIR="${DESTINATION_BASE}/results_bootstrap/${CONC}M/"
mkdir -p \$RESULTS_DIR
mv \$BASE_PATH/*.csv \$RESULTS_DIR/ 
mv \$BASE_PATH/*.pkl \$RESULTS_DIR/ 
mv \$BASE_PATH/*.png \$RESULTS_DIR/ 
mv \$BASE_PATH/*.log \$RESULTS_DIR/ 
echo "Results for Colvars analysis moved to \$RESULTS_DIR"

# Loop over different SKIP_FRAMES values
for n_skip in {0..9}; do
    SKIP_FRAMES=\$n_skip
    
    # Free Energy Analysis
    python free_energy_analysis2.py --base_path \$BASE_PATH --skip_frames \$SKIP_FRAMES --nstrides \$NSTRIDES --O_radii \$O_RADII --Cl_radii \$CL_RADII --Li_index \$LI_INDEX --T \$TEMP &> "free_energy_analysis_${SKIP_FRAMES}.log"
    echo "Free energy analysis with SKIP_FRAMES=\$SKIP_FRAMES completed!"

    # Free Energy Correction
    python free_energy_correction3.py --base_path \$BASE_PATH --nstrides \$NSTRIDES --O_radii \$O_RADII --H_radii \$H_RADII --Cl_radii \$CL_RADII --Li_index \$LI_INDEX --T \$TEMP --conc \$CONC &> "energy_correction_analysis_${SKIP_FRAMES}.log"
    echo "Free energy correction with SKIP_FRAMES=\$SKIP_FRAMES completed!"

    # Move results to results directory
    RESULTS_DIR="${DESTINATION_BASE}/results_bootstrap/${CONC}M/SKIP_\${SKIP_FRAMES}"
    mkdir -p \$RESULTS_DIR
    mv \$BASE_PATH/*.csv \$RESULTS_DIR/ 
    mv \$BASE_PATH/*.pkl \$RESULTS_DIR/ 
    mv \$BASE_PATH/*.png \$RESULTS_DIR/ 
    mv \$BASE_PATH/*.log \$RESULTS_DIR/ 
    echo "Results for SKIP_FRAMES=\$SKIP_FRAMES moved to \$RESULTS_DIR"
done
EOL

    # Make the generated script executable
    chmod +x "${DESTINATION_BASE}/${CONC}M/meta/analysis_bootstrap.sh"
    echo "Generated script: ${DESTINATION_BASE}/${CONC}M/meta/analysis_bootstrap.sh"

done
    
    echo "Copied Analysis scripts to ${CONC}M directory, meta"

done