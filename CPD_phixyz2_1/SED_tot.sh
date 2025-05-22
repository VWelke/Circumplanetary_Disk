#!/bin/bash
#SBATCH --job-name=CPD_dust_simulation
#SBATCH -p COMPUTE
#SBATCH --time=24:00:00
#SBATCH --job-name=nphi_116_1e8
#SBATCH -N1
#SBATCH -c 30
#SBATCH --mem-per-cpu=6G 
#SBATCH --mail-user=zcapvwe@ucl.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --chdir=/home/xzcapvwe/MSci_project/models/CPD_phixyz2_1_Image


# Load required modules
module load OpenMPI
# module load anaconda
echo "Checking the number of processors available..."
nproc
lscpu

# Ensure Python packages and radmc3d are accessible
export PATH=$PATH:~/bin



srun radmc3d sed incl 0 phi 0 setthreads 30

# Print completion message
echo "Job finished successfully at $(date)"

