#!/bin/bash
#SBATCH -p COMPUTE       # Request CPU partition
#SBATCH --time=10:00:00  # Set run time to 10 hours
#SBATCH -N1              # Use 1 node
#SBATCH -n16             # Use 16 CPU cores
#SBATCH --mail-user=zcapvwe@ucl.ac.uk 
#SBATCH --mail-type=ALL  # Get email notifications


# Change to the working directory that I've created and pasted my model.py 
#cd /home/xzcapvwe/MSci_project/models/CPD_dust_dias_ne9  

# Run RADMC-3D
srun python3 mockmodel.py
srun radmc3d mctherm setthreads 16  

# Make image 
srun radmc3d image incl 0 phi 0 loadlambda

