# solvation_spectra

### Here are the inputs and scripts for paper- Practical Corrections for Finite Concentrations Molecular Dynamics Simulations
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15377288.svg)](https://doi.org/10.5281/zenodo.15377288)

## Directory Structure
```
solvation spectra
├── README.md
├── data
└── analysis
```
## File Descriptions
```data```
- ```meta_single_Li```: Conduct Metadynamics (MetaD) simulations by applying bias potential on only one single Li+ ion.
For 3 temperatures (283K, 298K, 313K), we firstly conducted classical MD (```classical```) and then use the equilibrium trajectory to start MetaD simulations (```meta```).
  - In each concentration's ```meta```, there are 2 files to generate the inputs for the 10-replica multiple walker metaD. 
    - ```mk-folders-cp.sh```: Copy the base_inputs directory to 10 replica directories
    - ```submit.sh```: Replace the place holders in the scripts with the replica ids and submit jobs.
  - There is one more file for the post-analysis.
    - ```lammpstrj-cp.sh```: Copy and rename the lammps data files and lammps output trajectories for the following free energy analysis.
    - ```analysis-bootstrap.sh```

```analysis```
- Analyze Modules
  - ```python```
    - ```colvars_analysis_tool.py```: Metadynamics collective variable (CV) analyzer (equilibrium check)
    - ```free_energy_tool.py```: analyzer for free energy calculation and correction
- Scripts for each simulation condition (eg. 298K 0.5M)
  - ```colvars_analyzer_script.py```: Check the equilibration of Metadynamics by plotting the histgram and trajectory of CVs sampled over time, and plotting the potential of mean force (PMF)
  - ```free_energy_analysis.py```: Reweight the probability of configurations sampled by Metadynamics to get the free energy surface.
  - ```free_energy_correction.py```: Chemical potential correction.
- Scripts for collective analysis 
  - ```agglomerates_plot_script.py```: Plot agglomerates relative free energy at different concentrations and temperatures.
  - ```free_energy_error_bar_for_each_cluster.py```: Plot the free energies with error bars for specific cluster from the bootstrap analysis results.
  - ```free_energy_error_bar_for_each_conc.py```: Plot the free energies with error bars at specific concentration from the bootstrap analysis results.

## Example Workflow for 298K 
Noted that this workflow is currently only for LiCl water solutions. It can be easily adapted for other systems but careful edits to the scripts will be needed.
- Run Classical MD simulations:
  - ```cd ${project_dir}/298K/meta_single_Li/${CONC}}M/classical```
  - ```sbatch lammps.lammps.slurm```
- After finishing the Classical MD, we can create inputs and run MetaD simulations:
  - ```cd ${project_dir}/298K/meta_single_Li/${CONC}}M/meta```
  - ```sh mk-folders-cp.sh```
  - ```sbatch submit.sh```
- After finishing the MetaD, we can do data preparation.
  - ```sh lammpstrj-cp.sh```
- Generate bootstrap analysis scripts for each concentrations
  - ```cd ${project_dir}/298K/meta_single_Li```
  - ```sh copy_to_replicas_bootstrap.sh```
- Generate bootstrap analysis scripts for each concentrations
  - ```sh analysis_bootstrap.sh```
- Plotting with all the 298K results
  - ```python free_energy_bar_for_each_cluster.py```: Plot free energy of a specific cluster across different concentrations.
  - ```python free_energy_bar_for_each_conc.py```: Plot free energy of clusters at a specific concentration.
When you finish all the calculations
  - ```python agglomerates_plot_script.py```: Compare free energies of agglomerates across 3 temperatures.
