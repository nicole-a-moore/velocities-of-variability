#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

module load r/4.0.2 

Rscript 19_forcing-pops-marine.R 

