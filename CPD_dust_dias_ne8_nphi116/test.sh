#!/bin/bash
#SBATCH -p COMPUTE                      # Request CPU partition
#SBATCH -N1                             # Use 1 node
#SBATCH -n12                           # Use 16 CPU cores
#SBATCH --mail-user=e.edmondson@ucl.ac.uk 
#SBATCH --mail-type=ALL                   # Get email notifications
echo test

