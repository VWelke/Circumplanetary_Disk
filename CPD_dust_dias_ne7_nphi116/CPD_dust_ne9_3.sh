#!/bin/bash
#SBATCH --job-name=CPD_dust_simulation
#SBATCH -p COMPUTE
#SBATCH --time=10:00:00
#SBATCH --job-name=test123
#SBATCH -N1
#SBATCH -n16
#SBATCH --mail-user=zcapvwe@ucl.ac.uk
#SBATCH --mail-type=ALL

# Run Python setup script
srun python3 mockmodel.py

# Run RADMC-3D Monte Carlo simulation
srun radmc3d
srun radmc3d mctherm setthreads 16

# Generate image
srun radmc3d image incl 0 phi 0 loadlambda
