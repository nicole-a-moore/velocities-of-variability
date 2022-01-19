#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load java/13.0.2
  
wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20060101-20091231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20100101-20191231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20200101-20291231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20300101-20391231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20400101-20491231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20500101-20591231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20600101-20691231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20700101-20791231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20800101-20891231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password=SundayLab101! 'http://esgf-data1.ceda.ac.uk/thredds/fileServer/esg_dataroot/cmip5/output1/MPI-M/MPI-ESM-LR/rcp85/day/ocean/day/r2i1p1/v20111014/tos/tos_day_MPI-ESM-LR_rcp85_r2i1p1_20900101-21001231.nc'
