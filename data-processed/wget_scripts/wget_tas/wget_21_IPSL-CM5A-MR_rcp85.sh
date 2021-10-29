#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load java/13.0.2
  
wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://vesg.ipsl.upmc.fr/thredds/fileServer/cmip5/output1/IPSL/IPSL-CM5A-MR/rcp85/day/atmos/day/r1i1p1/v20111119/tas/tas_day_IPSL-CM5A-MR_rcp85_r1i1p1_20060101-20551231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://vesg.ipsl.upmc.fr/thredds/fileServer/cmip5/output1/IPSL/IPSL-CM5A-MR/rcp85/day/atmos/day/r1i1p1/v20111119/tas/tas_day_IPSL-CM5A-MR_rcp85_r1i1p1_20560101-21001231.nc'

