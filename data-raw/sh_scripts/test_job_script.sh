#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=2000M
#SBATCH --time=00:05:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load gcc/9.3.0 r/4.0.2

chosen_one=$1

Rscript test_script.R $chosen_one
