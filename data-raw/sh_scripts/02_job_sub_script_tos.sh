#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load gcc/9.3.0 r/4.0.2 netcdf/4.7.4

chosen_one=$1

Rscript 02_resampling-GCMs_tos.R $chosen_one
