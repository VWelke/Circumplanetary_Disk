#!/bin/bash
#SBATCH --job-name=CPD_dust_simulation
#SBATCH -p COMPUTE
#SBATCH --time=12:00:00
#SBATCH --job-name=test123
#SBATCH -N1
#SBATCH -c 16
#SBATCH --mail-user=zcapvwe@ucl.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --chdir=/home/xzcapvwe/MSci_project/models/CPD_dust_dias_ne9 


# Load required modules
module load OpenMPI
# module load anaconda
echo "Checking the number of processors available..."
nproc
lscpu

# Ensure Python packages and radmc3d are accessible
export PATH=$PATH:~/bin

# Run Python setup script
srun python3 mockmodel_grid_refined.py

# Run RADMC-3D Monte Carlo simulation
srun radmc3d mctherm setthreads 16

# Generate image
srun radmc3d image incl 0 phi 0 loadlambda

# Print completion message
echo "Job finished successfully at $(date)"

