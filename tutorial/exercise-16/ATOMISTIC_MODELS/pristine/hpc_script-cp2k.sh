#!/bin/bash
## Submission script for scarf
#SBATCH --comment=ALC_EQCM    # Project name
#SBATCH --job-name=Al2O3    # job name
#SBATCH -o %J.out
#SBATCH -e %J.err
#SBATCH --time=1-12:00:00        # days-hh:mm:ss
#
#SBATCH --partition=scarf    # queue (partition)
#SBATCH --ntasks=24          
#SBATCH --nodes=2           
#SBATCH --ntasks-per-node=12          
#SBATCH --cpus-per-task=2           
 
export OMP_NUM_THREADS=2           
export MKL_NUM_THREADS=2           
 
## Load required modules
module load module1
module load module2
 
## Define executable
exec="/home/vol08/scarf628/codes/cp2k.popt"
 
## Execute job
mpirun -np 24 $exec -i input.cp2k -o output.cp2k
