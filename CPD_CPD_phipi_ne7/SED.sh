#!/bin/bash
#SBATCH --job-name=SED
#SBATCH -p COMPUTE
#SBATCH --time=24:00:00
#SBATCH --job-name=SED
#SBATCH -N1
#SBATCH -c 4
#SBATCH --mem-per-cpu=6G 
#SBATCH --mail-user=zcapvwe@ucl.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --chdir=/home/xzcapvwe/MSci_project/models/CPD_sig2127 


# Load required modules
module load OpenMPI
# module load anaconda
echo "Checking the number of processors available..."
nproc
lscpu

# Ensure Python packages and radmc3d are accessible
export PATH=$PATH:~/bin

# Generate SED
srun radmc3d sed incl 0 phi 0 zoomau 36.2 38.2 -1 +1

# Print completion message
echo "Job finished successfully at $(date)"
