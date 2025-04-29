#!/bin/bash
#Job name:
#SBATCH --job-name=lammps
#
#Partition:
#SBATCH --partition=shared
#
#Account
#SBATCH --account=csd450
#
#Nodes
#SBATCH --nodes=1
#
#Processors:
#SBATCH --ntasks-per-node=32
#
#Wall clock limit:
#SBATCH --time=48:00:00

name=LiCl
rtemp=298
press=1
# rep_id=__REPID__
# restart_id=__RESTARTID__
# i=__i__
conc=5
lmp_equil_file=in.lammps
# lmp_data_file=lammps.298K.6000000.restart # 0.5M 

module purge
# module load cpu/0.15.4 slurm gcc/10.2.0 openmpi cmake gsl intel-mkl amdfftw
module load cpu/0.17.3b  gcc/10.2.0/npcyll4 cmake/3.21.4/teqow32 openmpi/4.1.3/oq3qvsv fftw/3.3.10/bmvqsbr python/3.8.12/7zdjza7 intel-mkl/2020.4.304/ghfk3mu gsl/2.7/wtlsmyy gnuplot/5.4.2/mfinpvw sdsc slurm

nprocs=$(( $SLURM_NTASKS_PER_NODE * $SLURM_NNODES / 2 ))
PARALLEL="mpirun -n $nprocs -mca btl vader,self"

# LMP="/expanse/lustre/projects/csd626/tpascal/codes/bin/lmp_expanse -screen none -var rtemp $rtemp -var press $press -var i $i -var rep_id $rep_id -var restart_id $restart_id"
# mkdir ${CURR_DIR}/results
LMP="/expanse/lustre/projects/csd626/tpascal/codes/bin/lmp_expanse -screen none -var rtemp $rtemp -var press $press"

echo "LAMMPS dynamics of ${conc}M ${name} lammps at ${rtemp}K"
# echo "running in $SCRATCH"
$PARALLEL $LMP -in ${lmp_equil_file} -log lammps.1LiCl.${rtemp}K.equil.lammps.log

# cp -r * ${CURR_DIR}/results
# rm -r $SCRATCH




