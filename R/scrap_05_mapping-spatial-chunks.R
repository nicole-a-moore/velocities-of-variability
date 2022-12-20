## visualizing the output:
library(raster)
library(tidyverse)
library(ncdf4)

## set up file names and paths
path = "/Volumes/NIKKI/CMIP5-GCMs/"

## create vector of file folders to put data into
gcm_models <- c("01_CMCC-CMS_tas",
                "02_GFDL-CM3_tas",
                "03_GFDL-ESM2G_tas",
                "04_HadGEM2-ES_tas",
                "05_inmcm4_tas",
                "06_IPSL-CM5A-MR_tas",
                "07_MIROC-ESM-CHEM_tas",
                "08_MIROC5_tas",
                "09_MPI-ESM-LR_tas",
                "10_MPI-ESM-MR_tas",
                "11_MRI-CGCM3_tas")

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
