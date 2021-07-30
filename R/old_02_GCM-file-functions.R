## functions for collating data from a GCM file into one clean dataset and detrending data at each location
library(tidyverse)
library(ncdf4)
library(abind)
library(akima)
library(raster)
library(imputeTS)
select <- dplyr::select

# filenames <- paste(path, c("tas_lon-1-120_lat-1-60.rds",
#                            "tas_lon-121-240_lat-1-60.rds",
#                            "tas_lon-241-360_lat-1-60.rds",
#                            "tas_lon-361-480_lat-1-60.rds",
#                            "tas_lon-1-120_lat-61-120.rds",
#                            "tas_lon-121-240_lat-61-120.rds",
#                            "tas_lon-241-360_lat-61-120.rds",
#                            "tas_lon-361-480_lat-61-120.rds",
#                            "tas_lon-1-120_lat-121-180.rds",
#                            "tas_lon-121-240_lat-121-180.rds",
#                            "tas_lon-241-360_lat-121-180.rds",
#                            "tas_lon-361-480_lat-121-180.rds",
#                            "tas_lon-1-120_lat-181-240.rds",
#                            "tas_lon-121-240_lat-181-240.rds",
#                            "tas_lon-241-360_lat-181-240.rds",
#                            "tas_lon-361-480_lat-181-240.rds"), sep = "")
# 
# filenames = paste(path, c("tas_lon-1-48_lat-1-24.rds",
#                           "tas_lon-49-96_lat-1-24.rds",
#                           "tas_lon-97-144_lat-1-24.rds",
#                           "tas_lon-145-192_lat-1-24.rds",
#                           "tas_lon-1-48_lat-25-48.rds",
#                           "tas_lon-49-96_lat-25-48.rds",
#                           "tas_lon-97-144_lat-25-48.rds",
#                           "tas_lon-145-192_lat-25-48.rds",
#                           "tas_lon-1-48_lat-49-72.rds",
#                           "tas_lon-49-96_lat-49-72.rds",
#                           "tas_lon-97-144_lat-49-72.rds",
#                           "tas_lon-145-192_lat-49-72.rds",
#                           "tas_lon-1-48_lat-73-96.rds",
#                           "tas_lon-49-96_lat-73-96.rds",
#                           "tas_lon-97-144_lat-73-96.rds",
#                           "tas_lon-145-192_lat-73-96.rds"), sep = "")

##############################################
#####              FUNCTIONS:           ######
##############################################

extract_and_organize <- function(historical_filenames, rcp85_filenames, path) {
  
  #####      EXTRACTING DATA FROM FILES      #####
  i = 1
  while (i < length(historical_filenames)+length(rcp85_filenames)+1) {
    
    if (i < length(historical_filenames)+1) {
      filename <- paste(path, historical_filenames[i], sep = "")
    }
    else {
      filename <- paste(path, rcp85_filenames[i-length(historical_filenames)], sep = "")
    }
    ncfile <- nc_open(filename)
    
    ## for first file from GCM, extract latitude, lonitude, air temp and time 
    if (i == 1) {
      lat <- ncvar_get(ncfile, "lat_bnds")
      lon <- ncvar_get(ncfile, "lon_bnds") 
      tas <- ncvar_get(ncfile, "tas") ## air surface temperature in an array of [lon, lat, time]
      time <- ncvar_get(ncfile, "time_bnds")
    }
    ## for remaining files, append air temp and time onto existing variables 
    else {
      tas <- abind(tas, ncvar_get(ncfile, "tas"), along = 3) 
      time_modified <- ncvar_get(ncfile, "time_bnds")
      ## modify time variable to link time series together 
      time_modified[1,] <-  time_modified[1,] + time[1,ncol(time)]+1
      time_modified[2,] <-  time_modified[2,] + time[1,ncol(time)]+1
      time <- cbind(time, time_modified)
    }
    
    ## close the file
    nc_close(ncfile)
    
    ## move to next file
    i = i + 1
  }
  
  ## save
  saveRDS(tas, paste(path, "tas.rds", sep = ""))
  
  ## check that time series have 84,371 days
  if (length(tas[1,1,]) != 84371) {
    print("uh oh! wrong time series length")
  }
  
  ## make lat and lon coords represent centre of grid cells:
  lat <- (lat[1,] + lat[2,]) / 2
  lon <- (lon[1,] + lon[2,]) / 2 
  
  ## reorganize data so North America is left of Europe when plotted:
  lon[which(lon >= 180)] = lon[which(lon >= 180)] - 360
  
  #####     REMOVING OUTLIERS/INTERPOLATING    #####
  ## loop through cells
  x = 1
  while (x < nrow(tas) + 1) {
    y = 1
    while (y < ncol(tas) + 1) {
      ## get local time series 
      local_ts <- tas[x,y,] 
      
      ## temperature remove values > 60C (333.15K)
      local_ts[which(local_ts > 333.15)] <- NA
      
      ## interpolate if temps are missing
      if(length(is.na(local_ts)) != 0 & length(is.na(local_ts)) != length(local_ts)) {
        local_ts <- na_kalman(local_ts, smooth = TRUE, model = "StructTS")
        tas[x,y,] <- local_ts
      }
      
      y = y + 1
    }
    x = x + 1
  }
  
  
  #####      STANDARDIZING SPATIAL EXTENT      #####
  ## must resample temperatures so data is on a 1 degree x 1 degree grid 
  ## doing it all in one object exhausts memory limit, so break into 10 time chunks, each spanning entire spatial extent:
  chunk = 1
  while (chunk < 11) {
    
    ## create empty array to store new spatially-standardized temps
    if(chunk == 10) {
      standardized_temps <- array(dim = c(360, 180, round(length(tas[1,1,])/10+1)))
    }
    else {
      standardized_temps <- array(dim = c(360, 180, round(length(tas[1,1,])/10)))
    }
    
    ## loop through each day
    day <- 1
    while (day < length(standardized_temps)/(180*360)+1) {
      ## extract temps for the day across all locations:
      daily_temps <-  expand.grid(lon, lat)
      colnames(daily_temps) <- c("lon", "lat") 
      daily_temps$temp <- as.vector(tas[,,day]) 
      
      ## resample temps onto 1x1 degree grid:
      interp <- interp(x = daily_temps$lon, y = daily_temps$lat, linear = FALSE,
                       z = daily_temps$temp, extrap = TRUE,
                       xo = seq(from = -179.5, to = 179.5), ## points represent centre of grid cells
                       yo = seq(from = -89.5, to = 89.5))
      
      ##interp_temps <- raster(interp)
      ##plot(interp_temps)
      
      ## add to new matrix:
      standardized_temps[,,day] <- interp[[3]]
      
      ## move to next day:
      day = day + 1
    }
    
    ## save standardized_temps to file:
    saveRDS(standardized_temps, paste(path, "time-chunk-", chunk, ".rds", sep = ""))
    
    ## move to next time chunk:
    chunk = chunk + 1
  }
  
  ## update lat and lon to reflect new standard values:
  lat <- seq(from = -89.5, to = 89.5)
  lon <- seq(from = -179.5, to = 179.5)
  
  ## clear unnecessary objects:
  rm(list = "tas")
  
  #####       REORGANIZE INTO SPATIAL CHUNKS      #####
  ## to detrend and perform sliding spectral analysis, need the whole time series for each location at once
  ## to satisfy this requirement while avoiding memory exhaustion, reorganize data
  ## go from time chunks spanning all of space to spatial chunks spanning all of time
  ## break into 8 x 60 degree lat x 60 degree lon chunks with data from 1850-2100
  lon_index = 1
  lat_index = 1
  count = 1
  while (lat_index < 180 & count <= 18) {
    
    ## get lat and lon bounds to extract in between:
    lon_bound1 <- lon_index
    lon_bound2 <- lon_index + 59
    lat_bound1 <- lat_index
    lat_bound2 <- lat_index + 59
    
    chunk = 1
    while (chunk < 11) {
      standardized_temps <- readRDS(paste(path, "time-chunk-",
                                          chunk, ".rds", sep = ""))
      
      ## extract temps within bounds and store:
      if(chunk == 1) {
        spatial_chunk <- standardized_temps[lon_bound1:lon_bound2,
                                            lat_bound1:lat_bound2, ]
      }
      else {
        spatial_chunk <- abind(spatial_chunk,
                               standardized_temps[lon_bound1:lon_bound2,
                                                  lat_bound1:lat_bound2, ],
                               along = 3)
      }
      
      chunk = chunk + 1
    }
    
    ## add name to list:
    if(count == 1) {
      filenames <- paste(path, "spatial-chunk_lon-", 
                         lon_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-", lat_bound2,
                         ".rds", sep = "")
    }
    else {
      filenames <- append(filenames, paste(path, "spatial-chunk_lon-", 
                                           lon_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-",
                                           lat_bound2,
                                           ".rds", sep = ""))
    }
    
    ## save spatial chunk:
    saveRDS(spatial_chunk, filenames[count])
    
    ## advance lat and lon indecies to move to next chunk
    if (count == 6) {
      lat_index <- lat_index + 60
      lon_index <- 1
    }
    else if (count == 12) {
      lat_index <- lat_index + 60
      lon_index <- 1
    }
    else {
      lon_index <- lon_index + 60
    }
    
    count = count + 1
  }
  
  ## returns list of spatial chunk filenames 
  return(filenames)
}

extract_and_organize_memory_conscious <- function(historical_filenames, rcp85_filenames, path) {
  
  #####      EXTRACTING DATA FROM FILES      #####
  filename <- paste(path, historical_filenames, sep = "")
  
  ## for first file from GCM, extract information about latitude and longitude
  raster = stack(filename)
  
  
  ncfile <- nc_open(filename)
  lat <- ncvar_get(ncfile, "lat_bnds")
  lon <- ncvar_get(ncfile, "lon_bnds") 
  
  ## define spatial chunks by dividing latitude and longitude 
  seq <- seq(from = 4, to = 10, by = 1)
  lat_div <- seq[first(which(ncfile$var$lat_bnds$size[[2]] %% seq == 0))]
  
  lat_dim <- seq(from = ncfile$var$lat_bnds$size[[2]]/lat_div, to = ncfile$var$lat_bnds$size[[2]], 
                 by = ncfile$var$lat_bnds$size[[2]]/lat_div)
  
  lon_div <- seq[first(which(ncfile$var$lon_bnds$size[[2]] %% seq == 0))]
  lon_dim <- seq(from = ncfile$var$lon_bnds$size[[2]]/lon_div, to = ncfile$var$lon_bnds$size[[2]], by = 
                   ncfile$var$lon_bnds$size[[2]]/lon_div)
  
  # close 
  nc_close(ncfile)
  
  lat_index <- 1
  lon_index <- 1
  count = 1
  while (lat_index < max(lat_dim)) {
    
    ## get lat and lon bounds to extract in between for this chunk:
    lon_bound1 <- lon_index
    lon_bound2 <- lon_index + lon_dim[1]-1
    lat_bound1 <- lat_index
    lat_bound2 <- lat_index + lat_dim[1]-1
    
    ## create vector of file names:
    filenames <- paste(path, append(historical_filenames, rcp85_filenames), sep = "")
    
    ## loop through GCM files and extract temps within bounds of the spatial chunk 
    file = 1
    while (file < length(historical_filenames)+length(rcp85_filenames)+1) {
      print(paste("On file:", file, sep = ""))
      ncfile <- nc_open(filenames[file])
      
      if(file == 1) {
        tas <- ncvar_get(ncfile, "tas")[lon_bound1:lon_bound2, lat_bound1:lat_bound2, ] ## air surface temperature in an array of [lon, lat, time]
        time <- ncvar_get(ncfile, "time_bnds")
      }
      else {
        tas <- abind(tas, ncvar_get(ncfile, "tas")[lon_bound1:lon_bound2, lat_bound1:lat_bound2, ],
                     along = 3) 
        time_modified <- ncvar_get(ncfile, "time_bnds")
        ## modify time variable to link time series together 
        time_modified[1,] <-  time_modified[1,] + time[1,ncol(time)]+1
        time_modified[2,] <-  time_modified[2,] + time[1,ncol(time)]+1
        time <- cbind(time, time_modified)
      }
      
      ## close the file
      nc_close(ncfile)
      
      ## move to next file
      file = file + 1
    }
    
    ## check that time series have 84,371 days
    if (length(tas[1,1,]) != 84371) {
      print("uh oh! wrong time series length")
      
      ## chop to correct length depending on GCM:
    }
    
    ## make lat and lon coords represent upper left corner  of grid cells:
    lat <- (lat[1,])
    lon <- (lon[1,])
    
    
    #   ## reorganize data so North America is left of Europe when plotted:
    #   ## update: no don't do this now
    #   ##lon[which(lon >= 180)] = lon[which(lon >= 180)] - 360
    
    
    #####     REMOVING OUTLIERS/INTERPOLATING MISSING TEMPS IN TIME    #####
    ## loop through cells
    x = 1
    while (x < nrow(tas) + 1) {
      y = 1
      while (y < ncol(tas) + 1) {
        ## get local time series 
        local_ts <- tas[x,y,] 
        
        ## temperature remove values > 60C (333.15K)
        local_ts[which(local_ts > 333.15)] <- NA
        
        ## interpolate if temps are missing
        if(length(is.na(local_ts)) != 0 & length(is.na(local_ts)) != length(local_ts)) {
          local_ts <- na_kalman(local_ts, smooth = TRUE, model = "StructTS")
          tas[x,y,] <- local_ts
        }
        
        y = y + 1
      }
      x = x + 1
    }
    
    ## add new filename to list:
    if(count == 1) {
      new_filenames <- paste(path, "tas_lon-", 
                         lon_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-", lat_bound2,
                         ".nc", sep = "")
    }
    else {
      new_filenames <- append(new_filenames, paste(path, "tas_lon-", 
                                           lon_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-",
                                           lat_bound2,
                                           ".rds", sep = ""))
    }
    
    ## save chunk as nc file:
    ncname <- new_filenames[count]
    dname <- "tas"
    
    # define dimensions
    xdim <- ncdim_def("lon", units = "degrees",
                      longname = "degrees longitude", as.double(lon[lon_bound1:lon_bound2]))
    ydim <- ncdim_def("lat", units = "degrees",
                      longname="degrees latitude", as.double(lat[lat_bound1:lat_bound2]))
    timedim <- ncdim_def("time", longname = "days since 1869-12-31", 
                         units = "days", as.double(time[1,]))
    
    # define surface air temperature variable
    fillvalue <- 1e30
    dlname <- "air temperature"
    tas_def <- ncvar_def("tas", units = 'degrees K', list(xdim,ydim,timedim),
                         fillvalue, dlname, prec ="double")
    
    ncout <- nc_create(ncname, vars = list(tas_def))
    
    # put variable
    ncvar_put(ncout,tas_def,tas)
  
    test <- nc_open(ncname)
    nc_close(test)
    
    saveRDS(tas, new_filenames[count])
    
    ## advance lat and long indecies to move to next chunk
    if (lon_index + lon_dim[1] > max(lon_dim)) {
      lat_index <- lat_index + lat_dim[1]
      lon_index <- 1
    }
    else {
      lon_index <- lon_index + lon_dim[1]
    }
    
    count = count + 1
  }
  
  ## get rid of big object 
  rm(list = "tas")
  
  
  #####      STANDARDIZING SPATIAL EXTENT      #####
  ## must resample temperatures so data is on a 1 degree x 1 degree grid of the same extent 
  ## doing it all in one object exhausts memory limit, so break into time chunks, each spanning entire spatial extent:
    
  ## set some variables for defining size of time chunk files:
  if (str_detect(rcp85_filenames[1], "CMCC-CESM")) {
    time_length =  round(84371/10)
    addition = 1
    chunks = 10
    times <- seq(from = 0, to = 84371, by = time_length) 
    times[length(times)] <- times[length(times)]+addition
  } 
  else {
    time_length =  round(84371/20)
    addition = -9
    chunks = 20
    times <- seq(from = 0, to = 84371+time_length, by = time_length) 
    times[length(times)] <- times[length(times)]+addition
  }
  
  chunk = 1
  while (chunk < chunks+1) {
    
    if(chunk == chunks) {
      ## create empty array to store new spatially-standardized temps
      standardized_temps <- array(dim = c(360, 180, time_length + addition)) ## deal with uneven division of time
      ## create new array to store temps while extracting from files:
      spatial_array <- array(dim = c(length(lon), length(lat), time_length + addition))
    }
    else {
      standardized_temps <- array(dim = c(360, 180, time_length))
      ## create new array to store temps:
      spatial_array <- array(dim = c(length(lon), length(lat), time_length))
    }
    
    ## loop through files and extract temperatures within time chunk: 
    lat_index <- 1
    lon_index <- 1
    file = 1 
    while (file < length(filenames)+1) {
      
      ## get lat and lon bounds to place temps between:
      lon_bound1 <- lon_index
      lon_bound2 <- lon_index + lon_dim[1]-1
      lat_bound1 <- lat_index
      lat_bound2 <- lat_index + lat_dim[1]-1
      
      ## read file:
      tas <- readRDS(filenames[file])[,,(times[chunk]+1):times[chunk+1]]
      
      ## extract temps:
      spatial_array[lon_bound1:lon_bound2, lat_bound1:lat_bound2,] <- tas
      
      ## advance lat and lon indecies to move to next chunk
      if (lon_index + lon_dim[1] > max(lon_dim)) {
        lat_index <- lat_index + lat_dim[1]
        lon_index <- 1
      }
      else {
        lon_index <- lon_index + lon_dim[1]
      }
      
      file = file + 1
    }
    
    rm(list = "tas")
    
    ## loop through each day of temps
    day <- 1
    while (day < length(standardized_temps)/(180*360)+1) {
      
      ## ensure reampling will not fail to interpolate edges of map by creating circularity in GCM
      c_spat_arr <- array(dim = c(length(lon) + 2, length(lat) + 2))
      c_spat_arr[2:(length(lon)+1), 2:(length(lat)+1)] <- spatial_array[,,day] ## reg map
      
      ## creating longitude circularity:
      ## copy right-most column to new empty left-most column, copy left-most column to new empty right-most column
      c_spat_arr[length(lon)+2, 2:(length(lat)+1)] <- spatial_array[1,,day] 
      c_spat_arr[1, 2:(length(lat)+1)] <- spatial_array[length(lon),,day] 
      ## creating latitude circulariy:
      ## copy bottom row to new empty top row, copy top row to new empty bottom row
      c_spat_arr[2:(length(lon)+1),length(lat) + 2] <- spatial_array[,1,day]
      c_spat_arr[2:(length(lon)+1), 1] <- spatial_array[,length(lat),day]  ## same idea here 
      
      ## fill in empty corners:
      c_spat_arr[1,1] <- spatial_array[length(lon),length(lat),day] 
      c_spat_arr[length(lon)+2,1] <- spatial_array[1,length(lat),day] 
      c_spat_arr[length(lon)+2,length(lat)+2] <- spatial_array[1,1,day] 
      c_spat_arr[1,length(lat)+2] <- spatial_array[length(lon),1,day]
      
      c_lon <- append(lon[1] - (lon[length(lon)] - lon[length(lon) - 1]),
                       append(lon, lon[length(lon)] + 
                                (lon[2] - lon[1])))
      c_lat <- append(lat[1] - (lat[length(lat)] - lat[length(lat) - 1]),
                      append(lat, lat[length(lat)] + 
                               (lat[2] - lat[1])))
      
      ## resample temps onto 1x1 degree grid by performing bilinear interpolation:
      interp <- bilinear.grid(x = c_lon,
                              y = c_lat,
                              z = c_spat_arr,
                              xlim = c(0.5, 359.5),
                              ylim = c(-89.5, 89.5),
                              dx = 1,
                              dy = 1)
      
      # b_df <- expand.grid(interp$x, interp$y)
      # colnames(b_df) <- c("lon", "lat")
      # b_df$temp <- as.vector(interp$z)
      # 
      # b_df$temp[which(b_df$temp == 0)] <- NA ## should be none since circularity was incorporated
      # 
      # b_df %>%
      #   ggplot() +
      #   geom_raster(aes(x = lon, y = lat, fill = temp))
      
      ## if zero values exist, print message that circularity did not work!
      if(length(which(interp[[3]] == 0)) != 0) {
        print("Oh no!!!! Introducing circularity was unsuccessful!!! :(")
      }
      
      ## add to new matrix:
      standardized_temps[,,day] <- interp[[3]]
      
      ## move to next day:
      print(paste("On day number: ", day, sep = ""))
      day = day + 1
    }
    
    ## save standardized_temps to file:
    saveRDS(standardized_temps, paste(path, "time-chunk-", chunk, ".rds", sep = ""))
    
    ## move to next time chunk:
    chunk = chunk + 1
  }
  
  
  ## update lat and lon to reflect new standard grid values:
  lat <- seq(from = -89.5, to = 89.5)
  lon <- seq(from = 0.5, to = 359.5)

  
  #####       REORGANIZE INTO SPATIAL CHUNKS      #####
  ## to detrend and perform sliding spectral analysis, need the whole time series for each location at once
  ## to satisfy this requirement while avoiding memory exhaustion, reorganize data
  ## go from time chunks spanning all of space to spatial chunks spanning all of time
  ## break into 8 x 60 degree lat x 60 degree lon chunks with data from 1850-2100
  lon_index = 1
  lat_index = 1
  count = 1
  while (lat_index < 180 & count <= 18) {
    
    ## get lat and lon bounds to extract in between:
    lon_bound1 <- lon_index
    lon_bound2 <- lon_index + 59
    lat_bound1 <- lat_index
    lat_bound2 <- lat_index + 59
    
    chunk = 1
    while (chunk < chunks+1) {
      standardized_temps <- readRDS(paste(path, "time-chunk-",
                                          chunk, ".rds", sep = ""))
      
      ## extract temps within bounds and store:
      if(chunk == 1) {
        spatial_chunk <- standardized_temps[lon_bound1:lon_bound2,
                                            lat_bound1:lat_bound2, ]
      }
      else {
        spatial_chunk <- abind(spatial_chunk,
                               standardized_temps[lon_bound1:lon_bound2,
                                                  lat_bound1:lat_bound2, ],
                               along = 3)
      }
      
      chunk = chunk + 1
    }
    
    ## add name to list:
    if(count == 1) {
      filenames <- paste(path, "spatial-chunk_long-", 
                         long_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-", lat_bound2,
                         ".rds", sep = "")
    }
    else {
      filenames <- append(filenames, paste(path, "spatial-chunk_long-", 
                                           lon_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-",
                                           lat_bound2,
                                           ".rds", sep = ""))
    }
    
    ## save spatial chunk:
    saveRDS(spatial_chunk, filenames[count])
    
    ## advance lat and lon indecies to move to next chunk
    if (count == 6) {
      lat_index <- lat_index + 60
      lon_index <- 1
    }
    else if (count == 12) {
      lat_index <- lat_index + 60
      lon_index <- 1
    }
    else {
      lon_index <- lon_index + 60
    }
    
    count = count + 1
  }
    
  ## returns list of spatial chunk filenames 
  return(filenames)
}


detrend_tas <- function(filenames) {
  ## get path
  path <- str_split_fixed(filenames[1], pattern = "spatial-chunk", n = 2)[1,1]
  
  lat <- seq(from = -89.5, to = 89.5)
  lon <- seq(from = -179.5, to = 179.5)
  
  ## create function for calculating date
  as_date <- function(x, origin = getOption("date_origin")){
    origin <- ifelse(is.null(origin), "1970-01-01", origin)
    as.Date(x, origin)
  }
  options(date_origin = "1869-12-31")
  
  lon_index = 1
  lat_index = 1
  count = 1
  while(count < length(filenames)+1) {
    filename <- filenames[count]
    
    ## retrieve spatial chunk:
    tas <- readRDS(filename)
    
    #####     LINEARLY AND SEASONALLY DETREND TIME SERIES IN EACH RASTER CELL OF CHUNK    #####
    l_detrended_tas <- tas
    s_detrended_tas <- tas
    
    ## get lat and lon bounds of chunk:
    lon_bound1 <- lon_index
    lon_bound2 <- lon_index + 59
    lat_bound1 <- lat_index
    lat_bound2 <- lat_index + 59
    
    x = 1 ## represents lonitude index
    while (x < 61) {
      y = 1 ## represents latitude index
      while (y < 61) {
        local_ts <- l_detrended_tas[x, y, ] ## get the local time series 
        
        ts_df <- data.frame(time = 1:length(local_ts), ## add simple time column representing days from 1850-01-01
                            temp = local_ts, 
                            date = as_date(1:length(local_ts))) %>% ## add a date column 
          mutate(year = str_split_fixed(.$date, pattern = "-", n = 2)[,1]) %>% ## add a year column
          group_by(year) %>% ## group by year
          do(mutate(., julian_date = seq(1:length(.$year)))) %>% ## add julian date column
          ungroup() %>%
          group_by(julian_date) %>%
          do(mutate(., temp_profile = mean(.$temp))) %>% ## compute mean temp for each day of year across all years
          ungroup() %>%
          mutate(s_detrended_temp = temp - temp_profile) ## create column represneting seasonally detrended data 
        
        ## run linear regression for grid cell
        l_output <- lm(ts_df, formula = temp ~ time)
        s_output <- lm(ts_df, formula = s_detrended_temp ~ time)
        
        ## extract residuals and add to detrended tas objects:
        l_detrended_tas[x, y, ] <- l_output$residuals
        s_detrended_tas[x, y, ] <- s_output$residuals
        
        print(paste("Detrending x = ",x, ", y = ", y, sep = ""))
        y = y + 1
      }
      x = x + 1
    }
    
    ## add name to list:
    if(count == 1) {
      l_filenames <- paste(path, "l-detrended-spatial-chunk_lon-", 
                           lon_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-", lat_bound2,
                           ".rds", sep = "")
      s_filenames <- paste(path, "s-detrended-spatial-chunk_lon-", 
                           lon_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-", lat_bound2,
                           ".rds", sep = "")
    }
    else {
      l_filenames <- append(l_filenames,
                            paste(path, "l-detrended-spatial-chunk_lon-", 
                                  lon_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-", lat_bound2,
                                  ".rds", sep = ""))
      s_filenames <- append(s_filenames,
                            paste(path, "s-detrended-spatial-chunk_lon-", 
                                  lon_bound1,"-", lon_bound2, "_lat-", lat_bound1, "-", lat_bound2,
                                  ".rds", sep = ""))
    }
    
    ## save detrended spatial chunks:
    saveRDS(l_detrended_tas, l_filenames[count])
    saveRDS(s_detrended_tas, s_filenames[count])
    
    ## advance lat and lon indecies to move to next chunk
    if (count == 6) {
      lat_index <- lat_index + 60
      lon_index <- 1
    }
    else if (count == 12) {
      lat_index <- lat_index + 60
      lon_index <- 1
    }
    else {
      lon_index <- lon_index + 60
    }
    
    count = count + 1
  }
  
  file_list <- list(l_filenames, s_filenames)
  
  return(file_list)
}

##############################################
##############################################
##############################################


## example data to test:
# 
# ## path to GCM files from my computer:
# path = "/Volumes/ADATA HV620/CMIP5-GCMs/01_CMCC-CESM/"
# 
# ## example file names:
# ## make list of filenames to open:
# cmcc_cesm_historical <- c("tas_day_CMCC-CESM_historical_r1i1p1_18700101-18741231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_18750101-18791231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_18800101-18841231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_18850101-18891231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_18900101-18941231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_18950101-18991231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19000101-19041231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19050101-19091231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19100101-19141231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19150101-19191231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19200101-19241231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19250101-19291231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19300101-19341231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19350101-19391231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19400101-19441231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19450101-19491231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19500101-19541231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19550101-19591231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19600101-19641231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19650101-19691231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19700101-19741231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19750101-19791231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19800101-19841231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19850101-19891231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19900101-19941231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_19950101-19991231.nc")
# 
# cmcc_cesm_rcp85 <- c("tas_day_CMCC-CESM_rcp85_r1i1p1_20000101-20041231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20050101-20051231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20060101-20151231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20160101-20251231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20260101-20351231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20360101-20451231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20460101-20551231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20560101-20651231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20660101-20751231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20760101-20851231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20860101-20951231.nc",
#                      "tas_day_CMCC-CESM_rcp85_r1i1p1_20960101-21001231.nc")
# 
# e_and_o <- extract_and_organize(path = path, historical_filenames = cmcc_cesm_historical, 
#                                 rcp85_filenames = cmcc_cesm_rcp85)
# 
# e_and_o <- readRDS("/Volumes/ADATA HV620/GCM-test/e_and_o.rds")
# 
# lat <- e_and_o[[1]]
# lon <- e_and_o[[2]]
# tas <- e_and_o[[3]]
# 
# d_tas <- detrend_tas(tas, lat, lon)
# 
# l_d_tas <- d_tas[[1]]
# s_d_tas <- d_tas[[2]]
# 
# filenames <- c("/Volumes/ADATA HV620/GCM-test/time-chunk-1.rds",
#                "/Volumes/ADATA HV620/GCM-test/time-chunk-2.rds",
#                "/Volumes/ADATA HV620/GCM-test/time-chunk-3.rds",
#                "/Volumes/ADATA HV620/GCM-test/time-chunk-4.rds",
#                "/Volumes/ADATA HV620/GCM-test/time-chunk-5.rds",
#                "/Volumes/ADATA HV620/GCM-test/time-chunk-6.rds",
#                "/Volumes/ADATA HV620/GCM-test/time-chunk-7.rds",
#                "/Volumes/ADATA HV620/GCM-test/time-chunk-8.rds",
#                "/Volumes/ADATA HV620/GCM-test/time-chunk-9.rds",
#                "/Volumes/ADATA HV620/GCM-test/time-chunk-10.rds"
#                )
