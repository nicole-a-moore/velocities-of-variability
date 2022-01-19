#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

sh wget_04_MPI-ESM-LR_historical.sh
sh wget_04_MPI-ESM-LR_rcp85.sh
sh wget_05_MPI-ESM-MR_historical.sh
sh wget_05_MPI-ESM-MR_rcp85.sh
sh wget_06_GFDL-ESM2G_rcp85.sh
sh wget_07_GFDL-CM3_historical.sh
sh wget_07_GFDL-CM3_rcp85.sh
sh wget_08_GFDL-ESM2M_rcp85.sh
sh wget_09_HadGEM2-CC_historical.sh
sh wget_09_HadGEM2-CC_rcp85.sh
sh wget_10_HadGEM2-ES_historical.sh
sh wget_10_HadGEM2-ES_rcp85.sh