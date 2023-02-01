#!/bin/bash
#SBATCH --account=def-jsunday
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:20:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user='nicole.moore@mail.mcgill.ca'

module load java/13.0.2
  
wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password= 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM-CHEM/historical/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM-CHEM_historical_r1i1p1_18700101-18891231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password= 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM-CHEM/historical/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM-CHEM_historical_r1i1p1_18900101-19091231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password= 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM-CHEM/historical/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19100101-19291231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password= 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM-CHEM/historical/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19300101-19491231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password= 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM-CHEM/historical/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19500101-19691231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password= 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM-CHEM/historical/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19700101-19891231.nc'

wget --user='https://esgf-node.llnl.gov/esgf-idp/openid/nicole_a_moore' --password= 'http://esgf-data1.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM-CHEM/historical/day/ocean/day/r1i1p1/v20111129/tos/tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19900101-20051231.nc'

