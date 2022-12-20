#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=32GB
#SBATCH --time=216:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load gcc/9.3.0 r/4.0.2 netcdf/4.7.4

chosen_one=$1

Rscript 04_calculating-spectral-exponent_temp.R $chosen_one
