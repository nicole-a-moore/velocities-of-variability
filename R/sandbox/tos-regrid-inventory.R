## checking output of regridded curvilinear tos
library(raster)

gcm_models <- c("01_CMCC-CESM", "02_CMCC-CM", '03_CMCC-CMS', '04_MPI-ESM-LR', '05_MPI-ESM-MR',
                "06_GFDL-ESM2G", '07_GFDL-CM3', '08_GFDL-ESM2M', '09_HadGEM2-CC', '10_HadGEM2-ES',
                "11_HadGEM2-AO", '12_IPSL-CM5A-LR', '13_IPSL-CM5B-LR', '14_MIROC5', '15_MIROC5-ESM-CHEM',
                '16_MIROC5-ESM', "17_inmcm4", '18_CNRM-CM5', "19_MRI-CGCM3", '20_MRI-ESM1',
                '21_IPSL-CM5A-MR')

## 01:
temps <- stack("data-raw/01_CMCC-CESM/regridded_tos_day_CMCC-CESM_historical_r1i1p1_18680101-18731231.nc")
plot(temps[[1]])

## 03:
temps <- stack("data-raw/03_CMCC-CMS/regridded_tos_day_CMCC-CMS_historical_r1i1p1_18700101-18791231.nc")
plot(temps[[1]])

## 04:
temps <- stack("data-raw/04_MPI-ESM-LR/regridded_tos_day_MPI-ESM-LR_rcp85_r2i1p1_20900101-21001231.nc")
plot(temps[[1]])

## 05:
temps <- stack("data-raw/05_MPI-ESM-MR/regridded_tos_day_MPI-ESM-MR_historical_r1i1p1_18700101-18701231.nc")
plot(temps[[1]])

## 06:
temps <- stack("data-raw/06_GFDL-ESM2G/regridded_tos_day_GFDL-ESM2G_historical_r1i1p1_18660101-18701231.nc")
plot(temps[[1]])

## 07:
temps <- stack("data-raw/07_GFDL-CM3/regridded_tos_day_GFDL-CM3_historical_r1i1p1_18700101-18741231.nc")
plot(temps[[1]])

##08:
## not yet done 

## 09:
temps <- stack("data-raw/09_HadGEM2-CC/")
plot(temps[[1]])

## 10:
temps <- stack("data-raw/10_HadGEM2-ES/")
plot(temps[[1]])