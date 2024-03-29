#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=64GB
#SBATCH --time=130:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

module load r/4.0.2 

Rscript 16_forcing-pops-cluster.R 

