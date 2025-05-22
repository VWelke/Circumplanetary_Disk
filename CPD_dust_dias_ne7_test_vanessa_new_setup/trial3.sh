#! /bin/bash
#SBATCH --job-name=trial3
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
echo “Hello World”
date
hostname
sleep 60
echo “Bye”
