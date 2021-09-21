## visualizing the output:
library(raster)
library(tidyverse)
library(ncdf4)

## set up file names and paths
path = "/Volumes/SundayLab/CMIP5-GCMs/"
gcm_models <- c("01_CMCC-CESM", "02_CMCC-CM", '03_CMCC-CMS', '04_MPI-ESM-LR', '05_MPI-ESM-MR',
                  "06_GFDL-ESM2G", '07_GFDL-CM3', '08_GFDL-ESM2M', '09_HadGEM2-CC', '10_HadGEM2-ES',
                  "11_HadGEM2-AO", '12_IPSL-CM5A-LR', '13_IPSL-CM5B-LR', '14_MIROC5', '15_MIROC5-ESM-CHEM',
                  '16_MIROC5-ESM', "17_inmcm4", '18_CNRM-CM5', "19_MRI-CGCM3", '20_MRI-ESM1',
                  '21_IPSL-CM5A-MR')

folders <- paste(path, gcm_models, "/", sep = "")

gcm = 1
while (gcm < 10) {
  filepath = paste(folders[gcm], "sp_files.rds", sep = "")
  names = readRDS(filepath)
  cut = str_split_fixed(names, "/spatial_temps", n = 2)
  names = paste(folders[gcm], "spatial_temps", cut[,2], sep = "")
  
  l_filenames <- str_replace_all(names, "spatial_temps", 'l-detrended')
  s_filenames <-  str_replace_all(names, "spatial_temps", 's-detrended')
  
  lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
  lon <-seq(from = 0.5, to = 359.5, length.out = 360) 
  
  dates <- readRDS(paste(folders[gcm], "date_new.rds", sep = ""))
  
  date = "05.28"
  lon_index = 0
  lat_index = 90
  count = 1
  while(count < length(l_filenames)+1) {
    
    ## get lat and lon bounds to extract in between:
    lon_bound1 <- lon_index 
    lon_bound2 <- lon_index + 60
    lat_bound1 <- lat_index
    lat_bound2 <- lat_index - 60
    
    num = which(paste("1871", date, sep = '.') == dates)
    seq = seq(from = num, to = 83950, by = 365)
    
    ## retrieve spatial chunk from nc file 
    l_open = nc_open(l_filenames[count])
    l_detrended_tas = ncvar_get(l_open, "var1_1")[,,seq]
    nc_close(l_open)
    
    s_open = nc_open(s_filenames[count])
    s_detrended_tas = ncvar_get(s_open, "var1_1")[,,seq]
    nc_close(s_open)
    
    r_ldtas =  brick(l_detrended_tas,
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    extent(r_ldtas) <- c(lon_bound1, lon_bound2, lat_bound2, lat_bound1)
    
    r_sdtas =  brick(s_detrended_tas, 
                     crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
    extent(r_sdtas) <- c(lon_bound1, lon_bound2, lat_bound2, lat_bound1)
    
    ## merge together
    if (count == 1) {
      l_mosaic  <- r_ldtas
      s_mosaic  <- r_sdtas
    }
    else {
      l_mosaic <- mosaic(l_mosaic, r_ldtas, fun = mean)
      s_mosaic <- mosaic(s_mosaic, r_sdtas, fun = mean)
    }
    
    if (count %in% c(6,12)) {
      lat_index <- lat_index - 60
      lon_index <- 0
    }
    else {
      lon_index <- lon_index + 60
    }
    
    print(paste("count: ", count, sep = ""))
    count = count + 1
  }
  
  saveRDS(s_mosaic, paste("vov-shiny/s_mosaic-", gcm_models[gcm], ".rds", sep = ""))
  saveRDS(l_mosaic, paste("vov-shiny/l_mosaic-", gcm_models[gcm], ".rds", sep = ""))
  
  gcm = gcm + 1
}






