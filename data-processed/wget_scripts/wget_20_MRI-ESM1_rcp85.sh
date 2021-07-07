#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load java/13.0.2
  
wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20060101-20151231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20160101-20251231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20260101-20351231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20360101-20451231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20460101-20551231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20560101-20651231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20660101-20751231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20760101-20851231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20860101-20951231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MRI/MRI-ESM1/rcp85/day/atmos/day/r1i1p1/v20140210/tas/tas_day_MRI-ESM1_rcp85_r1i1p1_20960101-21001231.nc'

