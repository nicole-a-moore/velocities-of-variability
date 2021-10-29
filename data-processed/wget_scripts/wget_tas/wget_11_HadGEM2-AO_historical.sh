#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load java/13.0.2
  
wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://aims3.llnl.gov/thredds/fileServer/cmip5_css02_data/cmip5/output1/NIMR-KMA/HadGEM2-AO/historical/day/atmos/day/r1i1p1/tas/1/tas_day_HadGEM2-AO_historical_r1i1p1_18600101-20051230.nc'

