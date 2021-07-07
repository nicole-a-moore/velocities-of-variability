### extract and organize GCMs
require(tidyverse)
require(raster)
require(easyNCDF)

resample <- function(historical_filenames, rcp85_filenames, path) {
  
  #####      EXTRACTING DATA FROM FILES      #####
  ## create vector of file names:
  filenames <- append(historical_filenames, rcp85_filenames)
  paths <- paste(path, filenames, sep = "")
  ## raster object of desired resolution/extent
  r <- raster(ymx = 90, ymn = -90, xmn = 0, xmx = 360, res = 1) 
  
  ## for each file in the set of GCM files:
  file = 1
  dates <- c()
  while (file < (length(filenames)+1)) {
    
    ## create raster stack to store data 
    og_temps <- stack(paths[file])
    dates <- append(dates, names(og_temps))
    
    #####      STANDARDIZE TO 1X1 DEGREE GRID   #####
    ## resample temperatures so all GCMs are on a 1 degree x 1 degree grid of the same extent
    start <- Sys.time()
    new_temps <- resample(og_temps, r, method = 'bilinear',
                          filename = paste(path, "resampled_", filenames[file], sep = ""))
    ## estimating time left
    end <- Sys.time()
    elapsed <- end - start
    expected <- elapsed[[1]] * (length(filenames)-file) 
    print(paste0("remaining: ", round(expected,3), " minutes = ", round(expected/60,3), " hours"))
    
    file = file + 1
  }
  dates <- str_replace_all(dates, "X","")
  return(dates)
}

# ## create fake dates:
# m <- str_split_fixed(dates, '\\.', n = 3)[,2]
# m <- m[1:(365*3+366)]
# d <- str_split_fixed(dates, '\\.', n = 3)[,3]
# d <- d[1:(365*3+366)]
# y <- c(1870:2100)
# md <- paste(m, d, sep = ".")
# date <- c()
# i=1
# while (i <= 250) {
#   cur_y <- c(rep(y[i], each = 365), rep(y[i+1], each = 365), rep(y[i+2], each = 366),
#              rep(y[i+3], each = 365))
#   date <- c(date, paste(cur_y, md, sep = "."))
#   i = i+4
# }
# date <- str_replace_all(date, "X","")
# date <- date[1:84371]


reorganize <- function(historical_filenames, rcp85_filenames, path, dates) {
  
  ## get paths and filenames of resampled GCMs
  filenames <- paste("resampled_", append(historical_filenames, rcp85_filenames), sep = "")
  paths <- paste(path, filenames, sep = "")
  
  ## figure out which dates to remove:
  b4_1871 <- which(as.numeric(str_split_fixed(dates, "\\.",n=2)[,1]) <= 1870)
  ly_days <- which(str_detect(dates, "02.29"))
  remove <- c(b4_1871, ly_days)
  date_new <- dates[-remove]
  
  
  #####       REORGANIZE INTO SPATIAL CHUNKS      #####
  ## to detrend and perform sliding spectral analysis, need the whole time series for each location at once
  ## to satisfy this requirement while avoiding memory exhaustion, reorganize data
  ## go from time chunks spanning all of space to spatial chunks spanning all of time
  ## break into 8 x 60 degree lat x 60 degree lon chunks with data from 1850-2100
  lon_index = 0
  lat_index = 90
  count = 1
  sp_files <- c()
  while (count <= 18) {
    ## get lat and lon bounds to extract in between:
    lon_bound1 <- lon_index 
    lon_bound2 <- lon_index + 60
    lat_bound1 <- lat_index
    lat_bound2 <- lat_index - 60
    
    ## loop through resampled GCM files, extracting data within spatial chunk and appending 
    file = 1
    while (file < length(filenames)+1) {
      start <- Sys.time()
      temps <- stack(paths[file])
      
      ## crop to new extent within bounds of spatial chunk
      extent <- extent(lon_bound1, lon_bound2, lat_bound2, lat_bound1)
      cropped <- crop(temps, extent)
      
      ## add cropped raster to rest from the GCM:
      if(file == 1) {
        spatial_temps <- cropped
      }
      else {
        spatial_temps <- stack(spatial_temps, cropped)
      }
      ## estimating time left
      end <- Sys.time()
      elapsed <- end - start
      expected <- elapsed[[1]] * (length(filenames)-file) 
      print(paste0("remaining: ", round(expected,3), " minutes = ", round(expected/60,3), " hours"))
      file = file + 1
    }
    
    #writeRaster(c, paste(path, "spatial-temps-test.nc", sep = ""))
    spatial_temps = stack(paste(path, "spatial-temps-test.nc", sep = ""))
    ## turn into an array:
    temps_df <- as.array(spatial_temps) 
    
    ## remove dates:
    temps_df <- temps_df[,,-remove] 
    
    ## loop through cells in spatial chunk, removing outliers and interpolating missing vals in time series
    l_detrended_temps <- s_detrended_temps <- temps_df
    x = 1 ## longitude
    while (x < ncol(temps_df) + 1) {
      y = 1 ## latitude
      while (y < nrow(temps_df)+1) {
        ## get local time series 
        local_ts <- data.frame(time = 1:(length(temps_df[1,1,])), 
                               temp = temps_df[x,y,])
        
        if (!length(which(is.na(local_ts$temp))) == nrow(local_ts)) {
          #####     REMOVING OUTLIERS, INTERPOLATING MISSING TEMPS    #####
          ## temperature remove values > 60C (333.15K)
          local_ts$temp[which(local_ts$temp > 333.15)] <- NA
          
          ## interpolate if temps are missing
          if(length(is.na(local_ts$temp)) != 0 & length(is.na(local_ts$temp)) 
             != length(local_ts$temp)) {
            local_ts$temp <- na_kalman(local_ts$temp, smooth = TRUE, model = "StructTS")
            temps_df[x,y,] <- local_ts$temp ## save changed time series if it needed changing
          }
          
          print(paste("Done removing outliers & interpolating x ",x,  
                      " y ", y, " of chunk #", count,sep = ""))
          
          #####     LINEARLY AND SEASONALLY DETREND TIME SERIES IN EACH RASTER CELL    #####
          ## create empty objects
          local_ts$md <- str_split_fixed(date_new, "\\.", n=2)[,2]
          
          ts_df <- local_ts %>%
            group_by(md) %>%
            do(mutate(., temp_profile = mean(.$temp))) %>% ## compute temp climatology for each day of year
            ungroup() %>%
            mutate(s_detrended_temp = temp - temp_profile) ## create column representing seasonally
          
          ## run linear regression for grid cell
          l_output <- lm(ts_df, formula = temp ~ time)
          s_output <- lm(ts_df, formula = s_detrended_temp ~ time)
          
          ## extract residuals and add to detrended temps objects:
          l_detrended_temps[x,y,] <- l_output$residuals
          s_detrended_temps[x,y,] <- s_output$residuals
          print(paste("Done detrending x ",x,  " y ", y, " of chunk #", count,sep = ""))
        }
        y = y+1
      }
      x = x + 1
    }
    
    ## save:
    sp_files[count] <- paste(path, "spatial-temps_lon-", lon_bound1,"-", lon_bound2,
                             "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = "")
    ArrayToNc(temps_df, file_path = paste(path, "detrended_lon-", lon_bound1,"-", lon_bound2,
                                           "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(l_detrended_temps, file_path = paste(path, "l-detrended_lon-", lon_bound1,"-", lon_bound2,
                                          "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(s_detrended_temps, file_path = paste(path, "s-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    
    ## advance lat and lon indecies to move to next spatial chunk
    if (count == 6) {
      lat_index <- lat_index - 60
      lon_index <- 1
    }
    else if (count == 12) {
      lat_index <- lat_index - 60
      lon_index <- 1
    }
    else {
      lon_index <- lon_index + 60
    }
    
    ## move to next spatial chunk 
    count = count + 1
  }
  
  ## returns list of spatial chunk filenames 
  return(sp_files)
}




## for each element in the list of historical_filenames, rcp85_filenames, path:
i = 1
while (i < length(gcms)+1) {
  element = gcms[[i]]
  
  hfn <- element[[1]]
  rfn <- element[[2]]
  p <- element[[3]]
  
  d = resample(historical_filenames = hfn, rcp85_filenames = rcp, path = p)
  
  reorganize = resample(historical_filenames = hfn, rcp85_filenames = rcp, path = p, dates = d)
  
  i = i + 1
}
