#!/bin/bash
#SBATCH --job-name=CPD_dust_simulation
#SBATCH -p COMPUTE
#SBATCH --time=20:00:00
#SBATCH --job-name=test123
#SBATCH -N1
#SBATCH -n 1
#SBATCH --mail-user=zcapvwe@ucl.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --chdir=/home/xzcapvwe/MSci_project/models/CPD_dust_dias_ne9 


# Load required modules
module load OPENMpi
#module load anaconda

# Check number of processors available
echo "Checking the number of processors available..."
nproc
lscpu

# Ensure radmc3d binary is accessible
export PATH=$PATH:~/bin

# Run Python setup script
srun python3 mockmodel.py

# Run RADMC-3D Monte Carlo simulation
srun radmc3d mctherm

# Generate image
srun radmc3d image incl 0 phi 0 loadlambda

# Print completion message
echo "Job finished successfully at $(date)"
