#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load java/13.0.2
  
wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://vesg.ipsl.upmc.fr/thredds/fileServer/cmip5/output1/IPSL/IPSL-CM5A-LR/rcp85/day/atmos/day/r1i1p1/v20111103/tas/tas_day_IPSL-CM5A-LR_rcp85_r1i1p1_20060101-22051231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://vesg.ipsl.upmc.fr/thredds/fileServer/cmip5/output1/IPSL/IPSL-CM5A-LR/rcp85/day/atmos/day/r1i1p1/v20111103/tas/tas_day_IPSL-CM5A-LR_rcp85_r1i1p1_22060101-23001231.nc'

