#!/bin/bash
#Job name:
#SBATCH --job-name=lammps
#
#Partition:
#SBATCH --partition=shared
#
#Account
#SBATCH --account=csd799
#
#Nodes
#SBATCH --nodes=1
#
#Processors:
#SBATCH --ntasks-per-node=8
#
#Wall clock limit:
#SBATCH --time=48:00:00

# SCRATCH=/expanse/lustre/scratch/xiruan/temp_project/$SLURM_JOB_ID
# CURR_DIR=`pwd`
# echo $CURR_DIR >> whereWasI_$SLURM_JOB_ID
# mkdir $SCRATCH
# cp * $SCRATCH
# cd $SCRATCH

name=LiCl
rtemp=298
press=1
rep_id=__REPID__
restart_id=__RESTARTID__
i=__i__
conc=18
lmp_equil_file=in_colvars.lammps
# lmp_data_file=lammps.298K.6000000.restart # 0.5M 

module reset # Do not purge, 'module purge' clears all modules, while 'module reset' restores the environment to a predefined set of default modules.
module load cpu/0.15.4  gcc/9.2.0  gsl/2.5  openmpi/3.1.6 fftw/3.3.8

nprocs=$(( $SLURM_NTASKS_PER_NODE * $SLURM_NNODES / 2 ))
PARALLEL="mpirun -n $nprocs -mca btl vader,self"
LMP="/home/xiruan/repos/lammps-2024-04-17/bins/lmp_mpi_colvars_new -screen none -var rtemp $rtemp -var press $press -var i $i -var rep_id $rep_id -var restart_id $restart_id"

echo "LAMMPS dynamics of ${conc}M ${name} lammps at ${rtemp}K"
# echo "running in $SCRATCH"
$PARALLEL $LMP -in ${lmp_equil_file} -log lammps.LiCl.${rtemp}K.equil.lammps.log

# cp -r * ${CURR_DIR}/results
# rm -r $SCRATCH




