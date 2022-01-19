#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load java/13.0.2
  
wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM/rcp85/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM_rcp85_r1i1p1_20060101-20251231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM/rcp85/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM_rcp85_r1i1p1_20260101-20451231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM/rcp85/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM_rcp85_r1i1p1_20460101-20651231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM/rcp85/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM_rcp85_r1i1p1_20660101-20851231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM/rcp85/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM_rcp85_r1i1p1_20860101-21001231.nc'

