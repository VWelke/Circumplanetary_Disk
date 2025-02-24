#!/bin/bash
#SBATCH --job-name=CPD_dust_simulation
#SBATCH --output=slurm-%j.out
#SBATCH --output=slurm-%j.err
#SBATCH -p COMPUTE
#SBATCH --time=10:00:00
#SBATCH --job-name=test123
#SBATCH -N1
#SBATCH -n12
#SBATCH --mail-user=zcapvwe@ucl.ac.uk
#SBATCH --mail-type=ALL

# Run Python setup script
#srun python3 mockmodel.py
echo "Job done"


#SBATCH --output=arrayJob_%A_%a.out
#SBATCH --error=arrayJob_%A_%a.err
