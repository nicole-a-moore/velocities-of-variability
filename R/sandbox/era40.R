### analyzing era 40 reanalysis of historical observations of air and sea surface temperature 
library(tidyverse)
library(ncdf4)
library(easyNCDF)
library(raster)
library(broom)
library(foreach)
library(doParallel)

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
numCores

registerDoParallel(numCores)  # use multicore, set to the number of our cores


#################################################
###            resample observations           ## 
#################################################
r <- raster(ymx = 90, ymn = -90, xmn = 0, xmx = 360, res = 1) 
names_tas <- c()
names_sst <- c()
i=1
while (i <= length(filenames)) {
  ## open file, get time variable
  file = paste(path, filenames[i], sep = "")
  nc <- nc_open(file)
  time <- ncvar_get(nc, "initial_time0") 
  nc_close(nc)
  
  ## rasterize sst and tas:
  sst <- stack(file, varname="SSTK_GDS4_SFC")
  tas <- stack(file, varname="2T_GDS4_SFC")
  
  ## get rid of ice
  ocean <- data.frame(rasterToPoints(raster("data-processed/masks/cmip5-ocean.grd")))
  ocean$x <- ifelse(ocean$x >= 0, ocean$x - 180, ocean$x + 180) + 180
  ocean <- rasterFromXYZ(ocean)
  
  sst <- crop(sst, ocean)
  
  if (i == 0) {
    ## get rid of first day since not real for sst 
    sst <- sst[[5:nlayers(sst)]]
    tas <- tas[[5:nlayers(tas)]]
    time <- time[5:length(time)]
  }
  
  if (any(str_detect(time, "02/29"))) {
    sst <- sst[[-which(str_detect(time, "02/29"))]] ## get rid of leap year day
    tas <- tas[[-which(str_detect(time, "02/29"))]]
    time <- time[-which(str_detect(time, "02/29"))]
  }

  ## resample to mean daily temp by averaging across samples witin a day (4)
  indices <- rep(1:(length(time)/4),each = 4)
  sst_daily <- stackApply(sst, indices, fun = mean)
  tas_daily <- stackApply(tas, indices, fun = mean)
  
  #####      STANDARDIZE TO 1X1 DEGREE GRID   #####
  ## raster object of desired resolution/extent
  if ((i+1) %% 10 != 0) {
    if((i+1) %% 10 == 1){
      rtas <- tas_daily
      rsst <- sst_daily
    }
    else {
      rtas <- stack(rtas, tas_daily)
      rsst <- stack(rsst, sst_daily)
    }
   
  }
  else {
    rtas <- resample(stack(rtas, tas_daily), r, method = 'bilinear', 
                     filename = paste("data-processed/ERA-40/resampled_tas_", 
                                      filenames[i], sep = ""), overwrite = T)
    rsst <- resample(stack(rsst, sst_daily), r, method = 'bilinear',
                     filename = paste("data-processed/ERA-40/resampled_sst_", filenames[i], sep = ""),
                     overwrite = T)
    names_tas <- append(names_tas, paste("data-processed/ERA-40/resampled_tas_", filenames[i], sep = ""))
    names_sst <- append(names_sst, paste("data-processed/ERA-40/resampled_sst_", filenames[i], sep = ""))
  }
  
  print(paste0("Done file number: ", i), stdout())
  i = i+1
}

saveRDS(names_sst, 
        "data-processed/ERA-40/sst_filenames.rds")
saveRDS(names_tas, 
        "data-processed/ERA-40/tas_filenames.rds")


## get paths and filenames 
files_sst <- readRDS("data-processed/ERA-40/sst_filenames.rds")
files_tas <- readRDS("data-processed/ERA-40/tas_filenames.rds")

## make date vector
## first date: 09/01/1957
## last date: 08/31/2002
dates <- paste(rep(seq(1957, 2002), each = 365), ".", rep(seq(1:365), 45), sep = "")
dates <- dates[245:(45*365+242)]

reorganize_GCM(filenames = files_sst, type = "sst")
reorganize_GCM(filenames = files_tas, type = "tas")

## detrend and make spatial chunks 
reorganize_GCM <- function(filenames, type) {

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
      temps <- stack(filenames[file])
      
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
      print(paste0("file: ", file))
      file = file + 1
    }
    
    ## turn into an array:
    temps_df <- as.array(spatial_temps)  
    
    ## loop through cells in spatial chunk, removing outliers and interpolating missing vals in time series
    l_detrended_temps <- s_detrended_temps <- temps_df
    x = 1 ## longitude
    while (x < ncol(temps_df) + 1) {
      y = 1 ## latitude
      while (y < nrow(temps_df)+1) {
        ## get local time series 
        local_ts <- data.frame(time = 1:(length(temps_df[1,1,])), 
                               temp = temps_df[y,x,])
        
        # ggplot(local_ts, aes(x = time, y = temp)) + geom_line() +
        #   geom_smooth(method = "lm")
        
        if (!length(which(is.na(local_ts$temp))) == nrow(local_ts)) {
          #####     REMOVING OUTLIERS, INTERPOLATING MISSING TEMPS    #####
          ## temperature remove values > 60C (333.15K)
          local_ts$temp[which(local_ts$temp > 333.15)] <- NA
          
          ## interpolate if temps are missing
          if(length(which(is.na(local_ts$temp))) != 0 & length(which(is.na(local_ts$temp))) 
             != length(local_ts$temp)) {
            local_ts$temp <- na_kalman(local_ts$temp, smooth = TRUE, model = "StructTS")
            temps_df[y,x,] <- local_ts$temp ## save changed time series if it needed changing
          }
          
          print(paste("Done removing outliers & interpolating x ",x,  
                      " y ", y, " of chunk #", count,sep = ""))
          
          #####     LINEARLY AND SEASONALLY DETREND TIME SERIES IN EACH RASTER CELL    #####
          ## create empty objects
          local_ts$md <- str_split_fixed(dates, "\\.", n=2)[,2]
          
          ts_df <- local_ts %>%
            group_by(md) %>%
            do(mutate(., temp_profile = mean(.$temp))) %>% ## compute temp climatology for each day of year
            ungroup() %>%
            mutate(s_detrended_temp = temp - temp_profile) %>% ## create column representing seasonally detrended
            arrange(., time)
          
          ## run linear regression for grid cell
          l_output <- lm(ts_df, formula = temp ~ time)
          s_output <- lm(ts_df, formula = s_detrended_temp ~ time)
          
          ## extract residuals and add to detrended temps objects:
          l_detrended_temps[y,x,] <- l_output$residuals
          s_detrended_temps[y,x,] <- s_output$residuals
          
          print(paste("Done detrending x ", x,  " y ", y, " of chunk #", count,sep = ""))
        }
        y = y + 1
      }
      x = x + 1
    }
    
    ## save:
    sp_files[count] <- paste("data-processed/ERA-40/era-40_", type, "_spatial_temps_lon-", lon_bound1,"-", lon_bound2,
                             "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = "")
    # ArrayToNc(temps_df, file_path = paste(path, "spatial_temps_lon-", lon_bound1,"-", lon_bound2,
    #                                       "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(l_detrended_temps, file_path = paste("data-processed/ERA-40/era-40_", type, "_l-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(s_detrended_temps, file_path = paste("data-processed/ERA-40/era-40_", type, "_s-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    
    ## advance lat and lon indecies to move to next spatial chunk
    if (count %in% c(6, 12)) {
      lat_index <- lat_index - 60
      lon_index <- 0
    }
    else {
      lon_index <- lon_index + 60
    }
    
    ## move to next spatial chunk 
    count = count + 1
  }
  
  saveRDS(sp_files, paste("data-processed/ERA-40/era-40_", type, "_sp_files.rds", sep = ""))
  ## returns list of spatial chunk filenames 
  return(sp_files)
}


## mark cells where lots of the time series is missing
lon_index = 0
lat_index = 90
count = 1
sp_files <- c()
list_of_chunks <- list()
while (count <= 18) {
  ## get lat and lon bounds to extract in between:
  lon_bound1 <- lon_index 
  lon_bound2 <- lon_index + 60
  lat_bound1 <- lat_index
  lat_bound2 <- lat_index - 60
  
  ## loop through resampled GCM files, extracting data within spatial chunk and appending 
  spatial_temps <- foreach(file=1:length(filenames), .combine=stack) %dopar% {
    temps <- stack(filenames[file])
    
    ## crop to new extent within bounds of spatial chunk
    extent <- extent(lon_bound1, lon_bound2, lat_bound2, lat_bound1)
    cropped <- crop(temps, extent)
    
    cropped
  }
  
  ## turn into an array:
  temps_df <- as.array(spatial_temps)  
  
  ## loop through cells in spatial chunk, counting missing values 
  missing_vals <- temps_df
  x = 1 ## longitude
  while (x < ncol(temps_df) + 1) {
    y = 1 ## latitude
    while (y < nrow(temps_df)+1) {
      ## get local time series 
      local_ts <- data.frame(time = 1:(length(temps_df[1,1,])), 
                             temp = temps_df[y,x,])
      
      # ggplot(local_ts, aes(x = time, y = temp)) + geom_line() +
      #   geom_smooth(method = "lm")
      
      ## count missing temps and save
      missing_vals[y,x,] <- length(which(is.na(local_ts$temp)))
      print(paste("Done counting x ", x,  " y ", y, " of chunk #", count,sep = ""))
      
      y = y + 1
    }
    x = x + 1
  }
  
  ## save:
  temp <- raster(missing_vals[,,1])
  extent(temp) = extent(lon_bound1, lon_bound2, lat_bound2, lat_bound1)
  
  saveRDS(temp, paste("data-processed/ERA-40/temp_", count, ".rds", sep = ""))
  
  list_of_chunks[[count]] = temp
  
  ## advance lat and lon indecies to move to next spatial chunk
  if (count %in% c(6, 12)) {
    lat_index <- lat_index - 60
    lon_index <- 0
  }
  else {
    lon_index <- lon_index + 60
  }
  
  ## move to next spatial chunk 
  count = count + 1
}

i=1
mosaic <- list_of_chunks[[i]]
while (i < length(list_of_chunks)) {
  mosaic <- mosaic(mosaic, list_of_chunks[[i+1]], fun = mean)
  i = i+1
}
plot(mosaic)
saveRDS(mosaic, "data-processed/ERA-40/missing-data-count.rds")


## calculate spectral exponent after getting rid of sea in tas
sp_files_tas <- readRDS("data-processed/ERA-40/era-40_tas_sp_files.rds")
sp_files_sst <- readRDS("data-processed/ERA-40/era-40_sst_sp_files.rds")

##############################################
#####              FUNCTIONS:           ######
##############################################
## function that creates a parabolic window for a given series
## uses Eke 2000, Eq. 6: W(j) = 1 - (2j/(N+1) - 1)^2 for j = 1,...,N
parabolic_window <- function(series, N) {
  j = c(1:N)
  W = c()
  for (i in j) {
    W[i] = 1 - ((2*j[i])/(N+1) - 1)^2
  }
  #plot(W)
  return(W)
}

## function that bridge detrends a windowed series
## calculates line connecting the first and last points of the series
## then subtracts that line from data 
bridge_detrender <- function(windowed_series, N) {
  ## regress to get equation of line:
  data <- data.frame(x = c(1, N), y = windowed_series[c(1,N)])
  eq = lm(y ~ x, data = data) 
  coeffs = eq$coeff
  #plot(windowed_series)
  #abline(a = windowed_series[1], b = coeffs[2])
  
  ## subtract the line from the data
  df <- data.frame(x = 1:N)
  predictions <- predict(eq, df)
  windowed_series = windowed_series - predictions
  #plot(windowed_series) 
  
  return(windowed_series)
}

## function to calculate spectral exponent over a time series window uisng periodogram
spectral_exponent_calculator_PSD <- function(ts_window, l) {
  # Fourier transform the time series window: 
  dft <- fft(ts_window)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral <- data.frame(freq = freq, power = amp^2)
  
  ## plot spectrum:
  # spectral %>%
  #   ggplot(aes(x = freq, y = power)) + geom_line() +
  #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")
  
  ## get estimate of spectral exponent over time series window:
  ## fit slope to low freqeuncies
  model_output_low <- spectral %>%
    filter(freq < 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  ## fit slope to high frequencies
  model_output_high <- spectral %>%
    filter(freq >= 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  ## fit slope for species with generation times < 1year
  model_output_all <- spectral %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  return(list(-model_output_low$estimate, -model_output_high$estimate, -model_output_all$estimate))
}

## function to calculate spectral exponent over a time series window uisng average wavelet coefficient method
spectral_exponent_calculator_AWC <- function(ts_window, N) {
  ## a. compute wavelet transform
  wavelets <- biwavelet::wt(data.frame(time = 1:N, val = ts_window), do.sig = F)
  
  ## b. calculate arithmetic mean with respect to the translation coefficient (b)
  data <- data.frame(avg_wavelet_coeff = rowMeans(sqrt(wavelets$power), na.rm = T), period = wavelets$period)
  
  ## plot average coefficients versus period on a log–log plot
  # ggplot(data = data, aes(x = period, y = avg_wavelet_coeff)) + 
  #   geom_point() +
  #   scale_x_log10() + scale_y_log10() + 
  #   theme_light() +
  #   labs(x = "Period", y = "Average wavelet coefficient")
  
  ## c. calculate beta by calculating slope 
  lm <- lm(log(avg_wavelet_coeff) ~ log(period), 
           data = data)
  ## slope equal to H + 1/2
  b = 2*(as.numeric(lm$coefficients[2]) - 1/2) + 1
  
  return(b)
}
## function to calculate spectral exponent over a time series within sliding windows of varying widths (from 5-10 years) with a time step of one year 
## takes as input:
##      - a detrended time series
##      - its latitude and longitude vectors 
## returns a list containing: 
##      - matrix of change in spectral exponents using each window width
##      - data frame of spectral exponents within each sliding window across all locations and for all window widths 
sliding_window_spec_exp <- function(names, type) {

  l_filenames <- str_replace_all(names, "spatial_temps", 'l-detrended')
  s_filenames <-  str_replace_all(names, "spatial_temps", 's-detrended')
  
  ## read in the ocean coordinates
  ocean <- read.csv("data-processed/masks/cmip5-ocean-coords.csv")
  ocean_coords <- paste(ocean$lat, ocean$lon)
  
  lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
  lon <-seq(from = 0.5, to = 359.5, length.out = 360) 
  
  lon_index = 0
  lat_index = 0
  count = 1
  while(count < length(l_filenames)+1) {
    
    ## retrieve spatial chunk from nc file
    l_open = nc_open(l_filenames[count])
    l_detrended_tas = ncvar_get(l_open, "var1_1")
    nc_close(l_open)
    
    s_open = nc_open(s_filenames[count])
    s_detrended_tas = ncvar_get(s_open, "var1_1")
    nc_close(s_open)
    
    spec_exp_list <- list()
    element <- 1
    x = 1 ## longitude
    while (x < ncol(l_detrended_tas)+1) {
      y = 1 ## latitude
      
      ## parallelize:
      spec_exp_list[[x]] <-  foreach(y=1:nrow(l_detrended_tas), combine=rbind)  %dopar% {
        ## if not in the ocean
        if (!paste(lat[y + lat_index], lon[x + lon_index]) %in% ocean_coords | type == "sst") {
          l_local_ts <- l_detrended_tas[y,x,] ## get the local detrended time series
          s_local_ts <- s_detrended_tas[y,x,]
          
          ## if no time series, skip to next latitude
          if (length(which(is.na(l_local_ts))) == 16423) {
            rep(NA, 13)
          }
          else {
            local_ts <- data.frame(time = 1:length(l_local_ts), ## add integer time (days from 1871.01.01)
                                   l_temp = l_local_ts,
                                   s_temp = s_local_ts,
                                   date = dates) %>% ## add a date column
              mutate(year = str_split_fixed(.$date, 
                                            pattern = "\\.", n = 2)[,1]) %>% ## add a year column
              group_by(year) ## group by year
            
            #########################################
            ##        SENSITIVITY ANALYSIS:        ##
            #########################################
            ## calculate spectral exponent using FFT over n year windows
            ## store spectral exponents and calculate slope
            all_windows <- list()
            n = 5
            while (n < 11) {
              year_start <- 1958
              year_stop <- 1958 + n - 1
              
              while (year_start <= (2002 - n)) {
                ## extract temps within time window
                ts_chunk <- filter(local_ts, year %in% year_start:year_stop)
                
                ## get length
                L = nrow(ts_chunk)
                
                ## preprocess the time series:
                ## a. subtracting mean
                ts_s <- ts_chunk$s_temp - mean(ts_chunk$s_temp)
                
                ## b. windowing - multiply by a parabolic window 
                window_s <- parabolic_window(series = ts_s, N = L)
                ts_s <- ts_s*window_s
                
                ## c. bridge detrending (endmatching)
                ## ie. subtracting from the data the line connecting the first and last points of the series
                ts_s <- bridge_detrender(windowed_series = ts_s, N = L)
                
                ## calculate spectral exponent in window using PSD and AWC methods
                s_exp_PSD <- spectral_exponent_calculator_PSD(ts_s, l = L)
                
                s_exp_PSD_low <- s_exp_PSD[[1]]
                s_exp_PSD_high <- s_exp_PSD[[2]]
                s_exp_PSD_all <- s_exp_PSD[[3]]
                
                if (n == 10) {
                  ## calculate on linearly-detrended time series
                  ts_l <- ts_chunk$l_temp - mean(ts_chunk$l_temp)
                  window_l <- parabolic_window(series = ts_l, N = L)
                  ts_l <- ts_l*window_l
                  ts_l <- bridge_detrender(windowed_series = ts_l, N = L)
                  
                  l_exp_PSD <- spectral_exponent_calculator_PSD(ts_l, l = L)
                  
                  l_exp_PSD_low <- l_exp_PSD[[1]]
                  l_exp_PSD_high <- l_exp_PSD[[2]]
                  l_exp_PSD_all <-  l_exp_PSD[[3]]
                  
                  ## and perform wavelet analysis
                  l_exp_AWC <- spectral_exponent_calculator_AWC(ts_l, N = L)
                  s_exp_AWC <- spectral_exponent_calculator_AWC(ts_s, N = L)
                  
                  ## store:
                  all_windows[[element]] <- c(l_exp_PSD_low, s_exp_PSD_low, l_exp_AWC, s_exp_AWC, 
                                                l_exp_PSD_high, s_exp_PSD_high, l_exp_PSD_all, s_exp_PSD_all, 
                                                year_start, year_stop, 
                                                lat[y + lat_index],
                                                lon[x + lon_index], paste(n, "years"))
                } 
                else {
                  ## store:
                  all_windows[[element]] <- c(NA, s_exp_PSD_low, NA, NA,
                                                NA, s_exp_PSD_high, NA, s_exp_PSD_all,
                                                year_start, year_stop, 
                                                lat[y + lat_index],
                                                lon[x + lon_index], paste(n, "years"))
                }
                
                ## move to next window
                year_start = year_stop + 1
                year_stop = year_stop + n 
                
                element = element + 1
              }
              
              ## move to next window width
              n = n + 1
            }
            all_windows <- do.call("rbind", all_windows)
            all_windows
          }
        }
        else {
          print(paste("Skipping ocean (lat = ", lat[y + lat_index],
                      ", lon = ", lon[x + lon_index], ")", sep = ""))
          
          all_windows <- rep(NA, 13)
          all_windows
        }
      }
      ## get rid of NA rows:
      temp <- do.call("rbind", spec_exp_list[[x]])
      spec_exp_list[[x]] = temp[rowSums(is.na(temp)) != ncol(temp), ]
      
      print(paste("Calculating spectral exponent for lon = ", lon[x + lon_index], sep = ""))
      
      x = x + 1
    }
    
    ## bind rows in list into data frame
    spec_exp_df <- data.frame(do.call(rbind, spec_exp_list), stringsAsFactors = FALSE)
    colnames(spec_exp_df) <- c("l_spec_exp_PSD_low", "s_spec_exp_PSD_low", "l_spec_exp_AWC", "s_spec_exp_AWC",
                               "l_spec_exp_PSD_high", "s_spec_exp_PSD_high", "l_spec_exp_PSD_all", "s_spec_exp_PSD_all",
                               "window_start_year",
                               "window_stop_year", "lat", "lon", "time_window_width")
    
    ## convert numbers to numeric
    spec_exp_df[,1:12] <- sapply(spec_exp_df[,1:12], as.numeric)
    
    ## regress spectral exponent and extract slope representing change in spectral exponent over time for each location and window width
    l_model_output_PSD_low <- spec_exp_df %>%
      filter(time_window_width == "10 years") %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = l_spec_exp_PSD_low ~ window_start_year))) %>%
      filter(term == "window_start_year")
    
    colnames(l_model_output_PSD_low)[5:8] <- paste("l", colnames(l_model_output_PSD_low)[5:8], "PSD_low", sep = "_")
    
    s_model_output_PSD_low <- spec_exp_df %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = s_spec_exp_PSD_low ~ window_start_year))) %>%
      filter(term == "window_start_year")
    
    colnames(s_model_output_PSD_low)[5:8] <- paste("s", colnames(s_model_output_PSD_low)[5:8], "PSD_low", sep = "_")
    
    l_model_output_PSD_high <- spec_exp_df %>%
      filter(time_window_width == "10 years") %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = l_spec_exp_PSD_high ~ window_start_year))) %>%
      filter(term == "window_start_year")
    
    colnames(l_model_output_PSD_high)[5:8] <- paste("l", colnames(l_model_output_PSD_high)[5:8], "PSD_high", sep = "_")
    
    s_model_output_PSD_high <- spec_exp_df %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = s_spec_exp_PSD_high ~ window_start_year))) %>%
      filter(term == "window_start_year")
    
    colnames(s_model_output_PSD_high)[5:8] <- paste("s", colnames(s_model_output_PSD_high)[5:8], "PSD_high", sep = "_")
    
    l_model_output_PSD_all <- spec_exp_df %>%
      filter(time_window_width == "10 years") %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = l_spec_exp_PSD_all ~ window_start_year))) %>%
      filter(term == "window_start_year")
    
    colnames(l_model_output_PSD_all)[5:8] <- paste("l", colnames(l_model_output_PSD_all)[5:8], "PSD_all", sep = "_")
    
    s_model_output_PSD_all <- spec_exp_df %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = s_spec_exp_PSD_all ~ window_start_year))) %>%
      filter(term == "window_start_year")
    
    colnames(s_model_output_PSD_all)[5:8] <- paste("s", colnames(s_model_output_PSD_all)[5:8], "PSD_all", sep = "_")
    
    l_model_output_AWC <- spec_exp_df %>%
      filter(time_window_width == "10 years") %>%
      group_by(lat, lon) %>%
      do(tidy(lm(., formula = l_spec_exp_AWC ~ window_start_year))) %>%
      filter(term == "window_start_year")
    
    colnames(l_model_output_AWC)[4:7] <- paste("l", colnames(l_model_output_AWC)[4:7], "AWC", sep = "_")
    l_model_output_AWC$time_window_width = "10 years"
    
    s_model_output_AWC <- spec_exp_df %>%
      group_by(lat, lon) %>%
      do(tidy(lm(., formula = s_spec_exp_AWC ~ window_start_year))) %>%
      filter(term == "window_start_year")
    
    colnames(s_model_output_AWC)[4:7] <- paste("s", colnames(s_model_output_AWC)[4:7], "AWC", sep = "_")
    s_model_output_AWC$time_window_width = "10 years"
    
    
    ## bind model output columns to spectral exponent data:
    spec_exp_df <- left_join(spec_exp_df, l_model_output_PSD_low) %>%
      left_join(., s_model_output_PSD_low) %>%
      left_join(., l_model_output_AWC) %>%
      left_join(., s_model_output_AWC) %>%
      left_join(., l_model_output_PSD_high) %>%
      left_join(., s_model_output_PSD_high) %>%
      left_join(., l_model_output_PSD_all) %>%
      left_join(., s_model_output_PSD_all)
    
    ## add filename to list:
    if(count == 1) {
      se_filenames <- paste("data-processed/ERA-40/era-40_",type,"_spec-exp_long-", 
                            lon_index,"-", lon_index + 60, "_lat-",
                            90-lat_index,"-", 90-lat_index-60,
                            ".csv", sep = "")
    }
    else {
      se_filenames <- append(se_filenames,
                             paste("data-processed/ERA-40/era-40_",type,"_spec-exp_long-", 
                                   lon_index,"-", lon_index + 60,"_lat-",
                                   90-lat_index,"-", 90-lat_index-60,  
                                   ".csv", sep = ""))
    }
    
    ## save spectral exponent for this chunk:
    write.csv(spec_exp_df, se_filenames[count], row.names = FALSE)
    
    ## advance lat and long indecies to move to next chunk
    if (count == 6) {
      lat_index <- lat_index + 60
      lon_index <- 0
    }
    else if (count == 12) {
      lat_index <- lat_index + 60
      lon_index <- 0
    }
    else {
      lon_index <- lon_index + 60
    }
    
    count = count + 1
  }
  
  saveRDS(se_filenames, paste("data-processed/ERA-40/era-40_", type, "_se_filenames.rds", sep = ""))
  
  return(se_filenames)
}

se_filenames_tas <- readRDS("data-processed/ERA-40/era-40_tas_se_filenames.rds")
se_filenames_sst <- readRDS("data-processed/ERA-40/era-40_sst_se_filenames.rds")

####################################################################
#####    transform spectral exponent data into a rasterStack  ######
####################################################################
## function to convert output from script 04 into 6 rasterStacks (one for each sliding window width, from 5-10 years) that can be used with the gVoCC functions 
create_rasterStack <- function(se_filenames, type) {
  
  ## combine all spectral exponent csvs into one big dataframe
  file = 1
  while (file < length(se_filenames) + 1) {
    if (file == 1) {
      spec_exp <- read.csv(se_filenames[file])
    }
    else {
      spec_exp <- rbind(spec_exp, read.csv(se_filenames[file]))
    }
    print(paste("Reading file #", file, "/", length(se_filenames), sep = ""))
    file = file + 1
  }
  
  ## change lon so it matches sea surface temperature 
  spec_exp <- spec_exp %>%
    mutate(lon = ifelse(lon >= 180, lon - 180, lon + 178))
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
              res = 1)
  
  ## reorder time_window_width so list elements are in order of increasing time window widths:
  spec_exp$time_window_width <- factor(spec_exp$time_window_width, levels = 
                                         c("5 years", "6 years", "7 years", "8 years",
                                           "9 years", "10 years"))
  
  ## split spectral exponent data by window widths:
  ww_split <- split(spec_exp, spec_exp$time_window_width)
  
  ## for each window width:
  i = 1 
  while (i < length(ww_split) + 1) {
    
    ## split by window_start_year (each window step), loop through them and create raster layer for each
    ww <- ww_split[[i]]
    step_split <- split(ww, ww$window_start_year)
    step = length(step_split) ## loop backwards since rasterStacks add layers to beginning
    while (step >= 1) {
      
      ## pull out each different slope estimate
      l_sp_df_PSD_low <- step_split[[step]] %>%
        select(lon, lat, l_spec_exp_PSD_low)
      
      s_sp_df_PSD_low <- step_split[[step]] %>%
        select(lon, lat, s_spec_exp_PSD_low) 
      
      l_sp_df_AWC <- step_split[[step]] %>%
        select(lon, lat, l_spec_exp_AWC)
      
      s_sp_df_AWC <- step_split[[step]] %>%
        select(lon, lat, s_spec_exp_AWC) 
      
      l_sp_df_PSD_high <- step_split[[step]] %>%
        select(lon, lat, l_spec_exp_PSD_high)
      
      s_sp_df_PSD_high <- step_split[[step]] %>%
        select(lon, lat, s_spec_exp_PSD_high) 
      
      l_sp_df_PSD_all <- step_split[[step]] %>%
        select(lon, lat, l_spec_exp_PSD_all)
      
      s_sp_df_PSD_all <- step_split[[step]] %>%
        select(lon, lat, s_spec_exp_PSD_all) 
      
      ## create raster layer:
      l_layer_PSD_low <- rasterFromXYZ(l_sp_df_PSD_low)
      s_layer_PSD_low <- rasterFromXYZ(s_sp_df_PSD_low)
      l_layer_AWC <- rasterFromXYZ(l_sp_df_AWC)
      s_layer_AWC <- rasterFromXYZ(s_sp_df_AWC)
      l_layer_PSD_high <- rasterFromXYZ(l_sp_df_PSD_high)
      s_layer_PSD_high <- rasterFromXYZ(s_sp_df_PSD_high)
      l_layer_PSD_all <- rasterFromXYZ(l_sp_df_PSD_all)
      s_layer_PSD_all <- rasterFromXYZ(s_sp_df_PSD_all)
      #plot(l_layer_PSD_high)
      
      ## add to temporary rasterstack:
      if (step == length(step_split)) {
        l_temp_stack_PSD_low <- l_layer_PSD_low
        s_temp_stack_PSD_low <- s_layer_PSD_low
        l_temp_stack_AWC <- l_layer_AWC
        s_temp_stack_AWC <- s_layer_AWC
        l_temp_stack_PSD_high <- l_layer_PSD_high
        s_temp_stack_PSD_high <- s_layer_PSD_high
        l_temp_stack_PSD_all <- l_layer_PSD_all
        s_temp_stack_PSD_all <- s_layer_PSD_all
      }
      else {
        l_temp_stack_PSD_low <- addLayer(l_layer_PSD_low, l_temp_stack_PSD_low)
        s_temp_stack_PSD_low <- addLayer(s_layer_PSD_low, s_temp_stack_PSD_low)
        l_temp_stack_AWC <- addLayer(l_layer_AWC, l_temp_stack_AWC)
        s_temp_stack_AWC <- addLayer(s_layer_AWC, s_temp_stack_AWC)
        l_temp_stack_PSD_high <- addLayer(l_layer_PSD_high, l_temp_stack_PSD_high)
        s_temp_stack_PSD_high <- addLayer(s_layer_PSD_high, s_temp_stack_PSD_high)
        l_temp_stack_PSD_all <- addLayer(l_layer_PSD_all, l_temp_stack_PSD_all)
        s_temp_stack_PSD_all <- addLayer(s_layer_PSD_all, s_temp_stack_PSD_all)
      }
      
      ## move to nested for loop
      step = step - 1
    }
    
    names(l_temp_stack_PSD_low) <- names(l_temp_stack_PSD_high) <- names(s_temp_stack_PSD_low) <- 
      names(s_temp_stack_PSD_high) <- names(l_temp_stack_AWC) <- names(s_temp_stack_AWC) <- 
      names(l_temp_stack_PSD_all) <- names(s_temp_stack_PSD_all)  <- 
      paste("window", 1:nlayers(l_temp_stack_PSD_low), sep = "_")
    
    ## save temporary raster stack: 
    if (i == 1) {
      l_stack_list_PSD_low <- list(l_temp_stack_PSD_low)
      s_stack_list_PSD_low <- list(s_temp_stack_PSD_low)
      l_stack_list_AWC <- list(l_temp_stack_AWC)
      s_stack_list_AWC <- list(s_temp_stack_AWC)
      l_stack_list_PSD_high <- list(l_temp_stack_PSD_high)
      s_stack_list_PSD_high <- list(s_temp_stack_PSD_high)
      l_stack_list_PSD_all <- list(l_temp_stack_PSD_all)
      s_stack_list_PSD_all <- list(s_temp_stack_PSD_all)
    }
    else {
      l_stack_list_PSD_low <- append(l_stack_list_PSD_low, l_temp_stack_PSD_low)
      s_stack_list_PSD_low <- append(s_stack_list_PSD_low, s_temp_stack_PSD_low)
      l_stack_list_AWC <- append(l_stack_list_AWC, l_temp_stack_AWC)
      s_stack_list_AWC <- append(s_stack_list_AWC, s_temp_stack_AWC)
      l_stack_list_PSD_high <- append(l_stack_list_PSD_high, l_temp_stack_PSD_high)
      s_stack_list_PSD_high <- append(s_stack_list_PSD_high, s_temp_stack_PSD_high)
      l_stack_list_PSD_all <- append(l_stack_list_PSD_all, l_temp_stack_PSD_all)
      s_stack_list_PSD_all <- append(s_stack_list_PSD_all, s_temp_stack_PSD_all)
    }
    
    ## move to next time window width
    i = i + 1
  }
  
  ## name the list items 
  names(l_stack_list_PSD_low) <- names(l_stack_list_PSD_high) <- names(s_stack_list_PSD_low) <-
    names(s_stack_list_PSD_high) <- names(l_stack_list_AWC) <- names(s_stack_list_PSD_high) <-
    names(s_stack_list_PSD_all) <- names(s_stack_list_PSD_all) <-
    names(ww_split)
  
  ## save the rasterstack 
  saveRDS(l_stack_list_PSD_low, paste("data-processed/ERA-40/era-40_", type, "_l_stack_list_PSD_low.rds", sep = ""))
  saveRDS(s_stack_list_PSD_low, paste("data-processed/ERA-40/era-40_", type, "_s_stack_list_PSD_low.rds", sep = ""))
  saveRDS(l_stack_list_AWC, paste("data-processed/ERA-40/era-40_", type, "_l_stack_list_AWC.rds", sep = ""))
  saveRDS(s_stack_list_AWC, paste("data-processed/ERA-40/era-40_", type, "_s_stack_list_AWC.rds", sep = ""))
  saveRDS(l_stack_list_PSD_high, paste("data-processed/ERA-40/era-40_", type, "_l_stack_list_PSD_high.rds", sep = ""))
  saveRDS(s_stack_list_PSD_high, paste("data-processed/ERA-40/era-40_", type, "_s_stack_list_PSD_high.rds", sep = ""))
  saveRDS(l_stack_list_PSD_all, paste("data-processed/ERA-40/era-40_", type, "_l_stack_list_PSD_all.rds", sep = ""))
  saveRDS(s_stack_list_PSD_all, paste("data-processed/ERA-40/era-40_", type, "_s_stack_list_PSD_all.rds", sep = ""))

  stacks <- list(l_stack_list_PSD_low, s_stack_list_PSD_low, 
                 l_stack_list_AWC, s_stack_list_AWC,
                 l_stack_list_PSD_high, s_stack_list_PSD_high,
                 l_stack_list_PSD_all, s_stack_list_PSD_all)
  
  ## return the 6 lists of rasterStacks
  return(stacks)
}



spec_exp %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", c(s_spec_exp_PSD_high, s_spec_exp_PSD_all, s_spec_exp_PSD_low,
                                                 s_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  group_by(window_start_year, exp_type) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_start_year, mean_spec_exp) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") + 
  geom_smooth(method = "lm")  

spec_exp$dataset <- "ERA40"

spec_exp_berk 
spec_exp_berk$dataset = "Berkley Earth"
spec_exp_berk = select(spec_exp_berk, -lat_lon)

spec= rbind(spec_exp_berk, spec_exp)

spec %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", c(s_spec_exp_PSD_high, s_spec_exp_PSD_all, s_spec_exp_PSD_low,
                                                 s_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  group_by(window_start_year, exp_type, dataset) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_start_year, mean_spec_exp, dataset) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = exp_type, shape = dataset)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent", shape = "", colour = "") + 
  geom_smooth(method = "lm")  + 
  scale_colour_discrete(labels = c("AWC method", "All frequencies", 
                                   "High frequencies", "Low frequencies")) 

spec %>%
  gather(key = "exp_type", value = "spec_exp", c(s_spec_exp_PSD_high, s_spec_exp_PSD_all, s_spec_exp_PSD_low,
                                                 s_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  group_by(window_start_year, exp_type, time_window_width, dataset) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_start_year, mean_spec_exp, dataset, time_window_width) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = exp_type, shape = dataset)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent", shape = "", colour = "") + 
  geom_smooth(method = "lm") + facet_wrap(~time_window_width) +
  scale_colour_discrete(labels = c("AWC method", "All frequencies", 
                                   "High frequencies", "Low frequencies")) 


#################################################
###                 filenames                  ## 
#################################################
path = "data-raw/ERA-40/"

filenames <- rev(c("e4oper.an.sfc.200208.moore543260.nc",
               "e4oper.an.sfc.200207.moore543260.nc",
               "e4oper.an.sfc.200206.moore543260.nc",
               "e4oper.an.sfc.200205.moore543260.nc",
               "e4oper.an.sfc.200204.moore543260.nc",
               "e4oper.an.sfc.200203.moore543260.nc",
               "e4oper.an.sfc.200202.moore543260.nc",
               "e4oper.an.sfc.200201.moore543260.nc",
               "e4oper.an.sfc.200112.moore543260.nc",
               "e4oper.an.sfc.200111.moore543260.nc",
               "e4oper.an.sfc.200110.moore543260.nc",
               "e4oper.an.sfc.200109.moore543260.nc",
               "e4oper.an.sfc.200108.moore543260.nc",
               "e4oper.an.sfc.200107.moore543260.nc",
               "e4oper.an.sfc.200106.moore543260.nc",
               "e4oper.an.sfc.200105.moore543260.nc",
               "e4oper.an.sfc.200104.moore543260.nc",
               "e4oper.an.sfc.200103.moore543260.nc",
               "e4oper.an.sfc.200102.moore543260.nc",
               "e4oper.an.sfc.200101.moore543260.nc",
               "e4oper.an.sfc.200012.moore543260.nc",
               "e4oper.an.sfc.200011.moore543260.nc",
               "e4oper.an.sfc.200010.moore543260.nc",
               "e4oper.an.sfc.200009.moore543260.nc",
               "e4oper.an.sfc.200008.moore543260.nc",
               "e4oper.an.sfc.200007.moore543260.nc",
               "e4oper.an.sfc.200006.moore543260.nc",
               "e4oper.an.sfc.200005.moore543260.nc",
               "e4oper.an.sfc.200004.moore543260.nc",
               "e4oper.an.sfc.200003.moore543260.nc",
               "e4oper.an.sfc.200002.moore543260.nc",
               "e4oper.an.sfc.200001.moore543260.nc",
               "e4oper.an.sfc.199912.moore543260.nc",
               "e4oper.an.sfc.199911.moore543260.nc",
               "e4oper.an.sfc.199910.moore543260.nc",
               "e4oper.an.sfc.199909.moore543260.nc",
               "e4oper.an.sfc.199908.moore543260.nc",
               "e4oper.an.sfc.199907.moore543260.nc",
               "e4oper.an.sfc.199906.moore543260.nc",
               "e4oper.an.sfc.199905.moore543260.nc",
               "e4oper.an.sfc.199904.moore543260.nc",
               "e4oper.an.sfc.199903.moore543260.nc",
               "e4oper.an.sfc.199902.moore543260.nc",
               "e4oper.an.sfc.199901.moore543260.nc",
               "e4oper.an.sfc.199812.moore543260.nc",
               "e4oper.an.sfc.199811.moore543260.nc",
               "e4oper.an.sfc.199810.moore543260.nc",
               "e4oper.an.sfc.199809.moore543260.nc",
               "e4oper.an.sfc.199808.moore543260.nc",
               "e4oper.an.sfc.199807.moore543260.nc",
               "e4oper.an.sfc.199806.moore543260.nc",
               "e4oper.an.sfc.199805.moore543260.nc",
               "e4oper.an.sfc.199804.moore543260.nc",
               "e4oper.an.sfc.199803.moore543260.nc",
               "e4oper.an.sfc.199802.moore543260.nc",
               "e4oper.an.sfc.199801.moore543260.nc",
               "e4oper.an.sfc.199712.moore543260.nc",
               "e4oper.an.sfc.199711.moore543260.nc",
               "e4oper.an.sfc.199710.moore543260.nc",
               "e4oper.an.sfc.199709.moore543260.nc",
               "e4oper.an.sfc.199708.moore543260.nc",
               "e4oper.an.sfc.199707.moore543260.nc",
               "e4oper.an.sfc.199706.moore543260.nc",
               "e4oper.an.sfc.199705.moore543260.nc",
               "e4oper.an.sfc.199704.moore543260.nc",
               "e4oper.an.sfc.199703.moore543260.nc",
               "e4oper.an.sfc.199702.moore543260.nc",
               "e4oper.an.sfc.199701.moore543260.nc",
               "e4oper.an.sfc.199612.moore543260.nc",
               "e4oper.an.sfc.199611.moore543260.nc",
               "e4oper.an.sfc.199610.moore543260.nc",
               "e4oper.an.sfc.199609.moore543260.nc",
               "e4oper.an.sfc.199608.moore543260.nc",
               "e4oper.an.sfc.199607.moore543260.nc",
               "e4oper.an.sfc.199606.moore543260.nc",
               "e4oper.an.sfc.199605.moore543260.nc",
               "e4oper.an.sfc.199604.moore543260.nc",
               "e4oper.an.sfc.199603.moore543260.nc",
               "e4oper.an.sfc.199602.moore543260.nc",
               "e4oper.an.sfc.199601.moore543260.nc",
               "e4oper.an.sfc.199512.moore543260.nc",
               "e4oper.an.sfc.199511.moore543260.nc",
               "e4oper.an.sfc.199510.moore543260.nc",
               "e4oper.an.sfc.199509.moore543260.nc",
               "e4oper.an.sfc.199508.moore543260.nc",
               "e4oper.an.sfc.199507.moore543260.nc",
               "e4oper.an.sfc.199506.moore543260.nc",
               "e4oper.an.sfc.199505.moore543260.nc",
               "e4oper.an.sfc.199504.moore543260.nc",
               "e4oper.an.sfc.199503.moore543260.nc",
               "e4oper.an.sfc.199502.moore543260.nc",
               "e4oper.an.sfc.199501.moore543260.nc",
               "e4oper.an.sfc.199412.moore543260.nc",
               "e4oper.an.sfc.199411.moore543260.nc",
               "e4oper.an.sfc.199410.moore543260.nc",
               "e4oper.an.sfc.199409.moore543260.nc",
               "e4oper.an.sfc.199408.moore543260.nc",
               "e4oper.an.sfc.199407.moore543260.nc",
               "e4oper.an.sfc.199406.moore543260.nc",
               "e4oper.an.sfc.199405.moore543260.nc",
               "e4oper.an.sfc.199404.moore543260.nc",
               "e4oper.an.sfc.199403.moore543260.nc",
               "e4oper.an.sfc.199402.moore543260.nc",
               "e4oper.an.sfc.199401.moore543260.nc",
               "e4oper.an.sfc.199312.moore543260.nc",
               "e4oper.an.sfc.199311.moore543260.nc",
               "e4oper.an.sfc.199310.moore543260.nc",
               "e4oper.an.sfc.199309.moore543260.nc",
               "e4oper.an.sfc.199308.moore543260.nc",
               "e4oper.an.sfc.199307.moore543260.nc",
               "e4oper.an.sfc.199306.moore543260.nc",
               "e4oper.an.sfc.199305.moore543260.nc",
               "e4oper.an.sfc.199304.moore543260.nc",
               "e4oper.an.sfc.199303.moore543260.nc",
               "e4oper.an.sfc.199302.moore543260.nc",
               "e4oper.an.sfc.199301.moore543260.nc",
               "e4oper.an.sfc.199212.moore543260.nc",
               "e4oper.an.sfc.199211.moore543260.nc",
               "e4oper.an.sfc.199210.moore543260.nc",
               "e4oper.an.sfc.199209.moore543260.nc",
               "e4oper.an.sfc.199208.moore543260.nc",
               "e4oper.an.sfc.199207.moore543260.nc",
               "e4oper.an.sfc.199206.moore543260.nc",
               "e4oper.an.sfc.199205.moore543260.nc",
               "e4oper.an.sfc.199204.moore543260.nc",
               "e4oper.an.sfc.199203.moore543260.nc",
               "e4oper.an.sfc.199202.moore543260.nc",
               "e4oper.an.sfc.199201.moore543260.nc",
               "e4oper.an.sfc.199112.moore543260.nc",
               "e4oper.an.sfc.199111.moore543260.nc",
               "e4oper.an.sfc.199110.moore543260.nc",
               "e4oper.an.sfc.199109.moore543260.nc",
               "e4oper.an.sfc.199108.moore543260.nc",
               "e4oper.an.sfc.199107.moore543260.nc",
               "e4oper.an.sfc.199106.moore543260.nc",
               "e4oper.an.sfc.199105.moore543260.nc",
               "e4oper.an.sfc.199104.moore543260.nc",
               "e4oper.an.sfc.199103.moore543260.nc",
               "e4oper.an.sfc.199102.moore543260.nc",
               "e4oper.an.sfc.199101.moore543260.nc",
               "e4oper.an.sfc.199012.moore543260.nc",
               "e4oper.an.sfc.199011.moore543260.nc",
               "e4oper.an.sfc.199010.moore543260.nc",
               "e4oper.an.sfc.199009.moore543260.nc",
               "e4oper.an.sfc.199008.moore543260.nc",
               "e4oper.an.sfc.199007.moore543260.nc",
               "e4oper.an.sfc.199006.moore543260.nc",
               "e4oper.an.sfc.199005.moore543260.nc",
               "e4oper.an.sfc.199004.moore543260.nc",
               "e4oper.an.sfc.199003.moore543260.nc",
               "e4oper.an.sfc.199002.moore543260.nc",
               "e4oper.an.sfc.199001.moore543260.nc",
               "e4oper.an.sfc.198912.moore543260.nc",
               "e4oper.an.sfc.198911.moore543260.nc",
               "e4oper.an.sfc.198910.moore543260.nc",
               "e4oper.an.sfc.198909.moore543260.nc",
               "e4oper.an.sfc.198908.moore543260.nc",
               "e4oper.an.sfc.198907.moore543260.nc",
               "e4oper.an.sfc.198906.moore543260.nc",
               "e4oper.an.sfc.198905.moore543260.nc",
               "e4oper.an.sfc.198904.moore543260.nc",
               "e4oper.an.sfc.198903.moore543260.nc",
               "e4oper.an.sfc.198902.moore543260.nc",
               "e4oper.an.sfc.198901.moore543260.nc",
               "e4oper.an.sfc.198812.moore543260.nc",
               "e4oper.an.sfc.198811.moore543260.nc",
               "e4oper.an.sfc.198810.moore543260.nc",
               "e4oper.an.sfc.198809.moore543260.nc",
               "e4oper.an.sfc.198808.moore543260.nc",
               "e4oper.an.sfc.198807.moore543260.nc",
               "e4oper.an.sfc.198806.moore543260.nc",
               "e4oper.an.sfc.198805.moore543260.nc",
               "e4oper.an.sfc.198804.moore543260.nc",
               "e4oper.an.sfc.198803.moore543260.nc",
               "e4oper.an.sfc.198802.moore543260.nc",
               "e4oper.an.sfc.198801.moore543260.nc",
               "e4oper.an.sfc.198712.moore543260.nc",
               "e4oper.an.sfc.198711.moore543260.nc",
               "e4oper.an.sfc.198710.moore543260.nc",
               "e4oper.an.sfc.198709.moore543260.nc",
               "e4oper.an.sfc.198708.moore543260.nc",
               "e4oper.an.sfc.198707.moore543260.nc",
               "e4oper.an.sfc.198706.moore543260.nc",
               "e4oper.an.sfc.198705.moore543260.nc",
               "e4oper.an.sfc.198704.moore543260.nc",
               "e4oper.an.sfc.198703.moore543260.nc",
               "e4oper.an.sfc.198702.moore543260.nc",
               "e4oper.an.sfc.198701.moore543260.nc",
               "e4oper.an.sfc.198612.moore543260.nc",
               "e4oper.an.sfc.198611.moore543260.nc",
               "e4oper.an.sfc.198610.moore543260.nc",
               "e4oper.an.sfc.198609.moore543260.nc",
               "e4oper.an.sfc.198608.moore543260.nc",
               "e4oper.an.sfc.198607.moore543260.nc",
               "e4oper.an.sfc.198606.moore543260.nc",
               "e4oper.an.sfc.198605.moore543260.nc",
               "e4oper.an.sfc.198604.moore543260.nc",
               "e4oper.an.sfc.198603.moore543260.nc",
               "e4oper.an.sfc.198602.moore543260.nc",
               "e4oper.an.sfc.198601.moore543260.nc",
               "e4oper.an.sfc.198512.moore543260.nc",
               "e4oper.an.sfc.198511.moore543260.nc",
               "e4oper.an.sfc.198510.moore543260.nc",
               "e4oper.an.sfc.198509.moore543260.nc",
               "e4oper.an.sfc.198508.moore543260.nc",
               "e4oper.an.sfc.198507.moore543260.nc",
               "e4oper.an.sfc.198506.moore543260.nc",
               "e4oper.an.sfc.198505.moore543260.nc",
               "e4oper.an.sfc.198504.moore543260.nc",
               "e4oper.an.sfc.198503.moore543260.nc",
               "e4oper.an.sfc.198502.moore543260.nc",
               "e4oper.an.sfc.198501.moore543260.nc",
               "e4oper.an.sfc.198412.moore543260.nc",
               "e4oper.an.sfc.198411.moore543260.nc",
               "e4oper.an.sfc.198410.moore543260.nc",
               "e4oper.an.sfc.198409.moore543260.nc",
               "e4oper.an.sfc.198408.moore543260.nc",
               "e4oper.an.sfc.198407.moore543260.nc",
               "e4oper.an.sfc.198406.moore543260.nc",
               "e4oper.an.sfc.198405.moore543260.nc",
               "e4oper.an.sfc.198404.moore543260.nc",
               "e4oper.an.sfc.198403.moore543260.nc",
               "e4oper.an.sfc.198402.moore543260.nc",
               "e4oper.an.sfc.198401.moore543260.nc",
               "e4oper.an.sfc.198312.moore543260.nc",
               "e4oper.an.sfc.198311.moore543260.nc",
               "e4oper.an.sfc.198310.moore543260.nc",
               "e4oper.an.sfc.198309.moore543260.nc",
               "e4oper.an.sfc.198308.moore543260.nc",
               "e4oper.an.sfc.198307.moore543260.nc",
               "e4oper.an.sfc.198306.moore543260.nc",
               "e4oper.an.sfc.198305.moore543260.nc",
               "e4oper.an.sfc.198304.moore543260.nc",
               "e4oper.an.sfc.198303.moore543260.nc",
               "e4oper.an.sfc.198302.moore543260.nc",
               "e4oper.an.sfc.198301.moore543260.nc",
               "e4oper.an.sfc.198212.moore543260.nc",
               "e4oper.an.sfc.198211.moore543260.nc",
               "e4oper.an.sfc.198210.moore543260.nc",
               "e4oper.an.sfc.198209.moore543260.nc",
               "e4oper.an.sfc.198208.moore543260.nc",
               "e4oper.an.sfc.198207.moore543260.nc",
               "e4oper.an.sfc.198206.moore543260.nc",
               "e4oper.an.sfc.198205.moore543260.nc",
               "e4oper.an.sfc.198204.moore543260.nc",
               "e4oper.an.sfc.198203.moore543260.nc",
               "e4oper.an.sfc.198202.moore543260.nc",
               "e4oper.an.sfc.198201.moore543260.nc",
               "e4oper.an.sfc.198112.moore543260.nc",
               "e4oper.an.sfc.198111.moore543260.nc",
               "e4oper.an.sfc.198110.moore543260.nc",
               "e4oper.an.sfc.198109.moore543260.nc",
               "e4oper.an.sfc.198108.moore543260.nc",
               "e4oper.an.sfc.198107.moore543260.nc",
               "e4oper.an.sfc.198106.moore543260.nc",
               "e4oper.an.sfc.198105.moore543260.nc",
               "e4oper.an.sfc.198104.moore543260.nc",
               "e4oper.an.sfc.198103.moore543260.nc",
               "e4oper.an.sfc.198102.moore543260.nc",
               "e4oper.an.sfc.198101.moore543260.nc",
               "e4oper.an.sfc.198012.moore543260.nc",
               "e4oper.an.sfc.198011.moore543260.nc",
               "e4oper.an.sfc.198010.moore543260.nc",
               "e4oper.an.sfc.198009.moore543260.nc",
               "e4oper.an.sfc.198008.moore543260.nc",
               "e4oper.an.sfc.198007.moore543260.nc",
               "e4oper.an.sfc.198006.moore543260.nc",
               "e4oper.an.sfc.198005.moore543260.nc",
               "e4oper.an.sfc.198004.moore543260.nc",
               "e4oper.an.sfc.198003.moore543260.nc",
               "e4oper.an.sfc.198002.moore543260.nc",
               "e4oper.an.sfc.198001.moore543260.nc",
               "e4oper.an.sfc.197912.moore543260.nc",
               "e4oper.an.sfc.197911.moore543260.nc",
               "e4oper.an.sfc.197910.moore543260.nc",
               "e4oper.an.sfc.197909.moore543260.nc",
               "e4oper.an.sfc.197908.moore543260.nc",
               "e4oper.an.sfc.197907.moore543260.nc",
               "e4oper.an.sfc.197906.moore543260.nc",
               "e4oper.an.sfc.197905.moore543260.nc",
               "e4oper.an.sfc.197904.moore543260.nc",
               "e4oper.an.sfc.197903.moore543260.nc",
               "e4oper.an.sfc.197902.moore543260.nc",
               "e4oper.an.sfc.197901.moore543260.nc",
               "e4oper.an.sfc.197812.moore543260.nc",
               "e4oper.an.sfc.197811.moore543260.nc",
               "e4oper.an.sfc.197810.moore543260.nc",
               "e4oper.an.sfc.197809.moore543260.nc",
               "e4oper.an.sfc.197808.moore543260.nc",
               "e4oper.an.sfc.197807.moore543260.nc",
               "e4oper.an.sfc.197806.moore543260.nc",
               "e4oper.an.sfc.197805.moore543260.nc",
               "e4oper.an.sfc.197804.moore543260.nc",
               "e4oper.an.sfc.197803.moore543260.nc",
               "e4oper.an.sfc.197802.moore543260.nc",
               "e4oper.an.sfc.197801.moore543260.nc",
               "e4oper.an.sfc.197712.moore543260.nc",
               "e4oper.an.sfc.197711.moore543260.nc",
               "e4oper.an.sfc.197710.moore543260.nc",
               "e4oper.an.sfc.197709.moore543260.nc",
               "e4oper.an.sfc.197708.moore543260.nc",
               "e4oper.an.sfc.197707.moore543260.nc",
               "e4oper.an.sfc.197706.moore543260.nc",
               "e4oper.an.sfc.197705.moore543260.nc",
               "e4oper.an.sfc.197704.moore543260.nc",
               "e4oper.an.sfc.197703.moore543260.nc",
               "e4oper.an.sfc.197702.moore543260.nc",
               "e4oper.an.sfc.197701.moore543260.nc",
               "e4oper.an.sfc.197612.moore543260.nc",
               "e4oper.an.sfc.197611.moore543260.nc",
               "e4oper.an.sfc.197610.moore543260.nc",
               "e4oper.an.sfc.197609.moore543260.nc",
               "e4oper.an.sfc.197608.moore543260.nc",
               "e4oper.an.sfc.197607.moore543260.nc",
               "e4oper.an.sfc.197606.moore543260.nc",
               "e4oper.an.sfc.197605.moore543260.nc",
               "e4oper.an.sfc.197604.moore543260.nc",
               "e4oper.an.sfc.197603.moore543260.nc",
               "e4oper.an.sfc.197602.moore543260.nc",
               "e4oper.an.sfc.197601.moore543260.nc",
               "e4oper.an.sfc.197512.moore543260.nc",
               "e4oper.an.sfc.197511.moore543260.nc",
               "e4oper.an.sfc.197510.moore543260.nc",
               "e4oper.an.sfc.197509.moore543260.nc",
               "e4oper.an.sfc.197508.moore543260.nc",
               "e4oper.an.sfc.197507.moore543260.nc",
               "e4oper.an.sfc.197506.moore543260.nc",
               "e4oper.an.sfc.197505.moore543260.nc",
               "e4oper.an.sfc.197504.moore543260.nc",
               "e4oper.an.sfc.197503.moore543260.nc",
               "e4oper.an.sfc.197502.moore543260.nc",
               "e4oper.an.sfc.197501.moore543260.nc",
               "e4oper.an.sfc.197412.moore543260.nc",
               "e4oper.an.sfc.197411.moore543260.nc",
               "e4oper.an.sfc.197410.moore543260.nc",
               "e4oper.an.sfc.197409.moore543260.nc",
               "e4oper.an.sfc.197408.moore543260.nc",
               "e4oper.an.sfc.197407.moore543260.nc",
               "e4oper.an.sfc.197406.moore543260.nc",
               "e4oper.an.sfc.197405.moore543260.nc",
               "e4oper.an.sfc.197404.moore543260.nc",
               "e4oper.an.sfc.197403.moore543260.nc",
               "e4oper.an.sfc.197402.moore543260.nc",
               "e4oper.an.sfc.197401.moore543260.nc",
               "e4oper.an.sfc.197312.moore543260.nc",
               "e4oper.an.sfc.197311.moore543260.nc",
               "e4oper.an.sfc.197310.moore543260.nc",
               "e4oper.an.sfc.197309.moore543260.nc",
               "e4oper.an.sfc.197308.moore543260.nc",
               "e4oper.an.sfc.197307.moore543260.nc",
               "e4oper.an.sfc.197306.moore543260.nc",
               "e4oper.an.sfc.197305.moore543260.nc",
               "e4oper.an.sfc.197304.moore543260.nc",
               "e4oper.an.sfc.197303.moore543260.nc",
               "e4oper.an.sfc.197302.moore543260.nc",
               "e4oper.an.sfc.197301.moore543260.nc",
               "e4oper.an.sfc.197212.moore543260.nc",
               "e4oper.an.sfc.197211.moore543260.nc",
               "e4oper.an.sfc.197210.moore543260.nc",
               "e4oper.an.sfc.197209.moore543260.nc",
               "e4oper.an.sfc.197208.moore543260.nc",
               "e4oper.an.sfc.197207.moore543260.nc",
               "e4oper.an.sfc.197206.moore543260.nc",
               "e4oper.an.sfc.197205.moore543260.nc",
               "e4oper.an.sfc.197204.moore543260.nc",
               "e4oper.an.sfc.197203.moore543260.nc",
               "e4oper.an.sfc.197202.moore543260.nc",
               "e4oper.an.sfc.197201.moore543260.nc",
               "e4oper.an.sfc.197112.moore543260.nc",
               "e4oper.an.sfc.197111.moore543260.nc",
               "e4oper.an.sfc.197110.moore543260.nc",
               "e4oper.an.sfc.197109.moore543260.nc",
               "e4oper.an.sfc.197108.moore543260.nc",
               "e4oper.an.sfc.197107.moore543260.nc",
               "e4oper.an.sfc.197106.moore543260.nc",
               "e4oper.an.sfc.197105.moore543260.nc",
               "e4oper.an.sfc.197104.moore543260.nc",
               "e4oper.an.sfc.197103.moore543260.nc",
               "e4oper.an.sfc.197102.moore543260.nc",
               "e4oper.an.sfc.197101.moore543260.nc",
               "e4oper.an.sfc.197012.moore543260.nc",
               "e4oper.an.sfc.197011.moore543260.nc",
               "e4oper.an.sfc.197010.moore543260.nc",
               "e4oper.an.sfc.197009.moore543260.nc",
               "e4oper.an.sfc.197008.moore543260.nc",
               "e4oper.an.sfc.197007.moore543260.nc",
               "e4oper.an.sfc.197006.moore543260.nc",
               "e4oper.an.sfc.197005.moore543260.nc",
               "e4oper.an.sfc.197004.moore543260.nc",
               "e4oper.an.sfc.197003.moore543260.nc",
               "e4oper.an.sfc.197002.moore543260.nc",
               "e4oper.an.sfc.197001.moore543260.nc",
               "e4oper.an.sfc.196912.moore543260.nc",
               "e4oper.an.sfc.196911.moore543260.nc",
               "e4oper.an.sfc.196910.moore543260.nc",
               "e4oper.an.sfc.196909.moore543260.nc",
               "e4oper.an.sfc.196908.moore543260.nc",
               "e4oper.an.sfc.196907.moore543260.nc",
               "e4oper.an.sfc.196906.moore543260.nc",
               "e4oper.an.sfc.196905.moore543260.nc",
               "e4oper.an.sfc.196904.moore543260.nc",
               "e4oper.an.sfc.196903.moore543260.nc",
               "e4oper.an.sfc.196902.moore543260.nc",
               "e4oper.an.sfc.196901.moore543260.nc",
               "e4oper.an.sfc.196812.moore543260.nc",
               "e4oper.an.sfc.196811.moore543260.nc",
               "e4oper.an.sfc.196810.moore543260.nc",
               "e4oper.an.sfc.196809.moore543260.nc",
               "e4oper.an.sfc.196808.moore543260.nc",
               "e4oper.an.sfc.196807.moore543260.nc",
               "e4oper.an.sfc.196806.moore543260.nc",
               "e4oper.an.sfc.196805.moore543260.nc",
               "e4oper.an.sfc.196804.moore543260.nc",
               "e4oper.an.sfc.196803.moore543260.nc",
               "e4oper.an.sfc.196802.moore543260.nc",
               "e4oper.an.sfc.196801.moore543260.nc",
               "e4oper.an.sfc.196712.moore543260.nc",
               "e4oper.an.sfc.196711.moore543260.nc",
               "e4oper.an.sfc.196710.moore543260.nc",
               "e4oper.an.sfc.196709.moore543260.nc",
               "e4oper.an.sfc.196708.moore543260.nc",
               "e4oper.an.sfc.196707.moore543260.nc",
               "e4oper.an.sfc.196706.moore543260.nc",
               "e4oper.an.sfc.196705.moore543260.nc",
               "e4oper.an.sfc.196704.moore543260.nc",
               "e4oper.an.sfc.196703.moore543260.nc",
               "e4oper.an.sfc.196702.moore543260.nc",
               "e4oper.an.sfc.196701.moore543260.nc",
               "e4oper.an.sfc.196612.moore543260.nc",
               "e4oper.an.sfc.196611.moore543260.nc",
               "e4oper.an.sfc.196610.moore543260.nc",
               "e4oper.an.sfc.196609.moore543260.nc",
               "e4oper.an.sfc.196608.moore543260.nc",
               "e4oper.an.sfc.196607.moore543260.nc",
               "e4oper.an.sfc.196606.moore543260.nc",
               "e4oper.an.sfc.196605.moore543260.nc",
               "e4oper.an.sfc.196604.moore543260.nc",
               "e4oper.an.sfc.196603.moore543260.nc",
               "e4oper.an.sfc.196602.moore543260.nc",
               "e4oper.an.sfc.196601.moore543260.nc",
               "e4oper.an.sfc.196512.moore543260.nc",
               "e4oper.an.sfc.196511.moore543260.nc",
               "e4oper.an.sfc.196510.moore543260.nc",
               "e4oper.an.sfc.196509.moore543260.nc",
               "e4oper.an.sfc.196508.moore543260.nc",
               "e4oper.an.sfc.196507.moore543260.nc",
               "e4oper.an.sfc.196506.moore543260.nc",
               "e4oper.an.sfc.196505.moore543260.nc",
               "e4oper.an.sfc.196504.moore543260.nc",
               "e4oper.an.sfc.196503.moore543260.nc",
               "e4oper.an.sfc.196502.moore543260.nc",
               "e4oper.an.sfc.196501.moore543260.nc",
               "e4oper.an.sfc.196412.moore543260.nc",
               "e4oper.an.sfc.196411.moore543260.nc",
               "e4oper.an.sfc.196410.moore543260.nc",
               "e4oper.an.sfc.196409.moore543260.nc",
               "e4oper.an.sfc.196408.moore543260.nc",
               "e4oper.an.sfc.196407.moore543260.nc",
               "e4oper.an.sfc.196406.moore543260.nc",
               "e4oper.an.sfc.196405.moore543260.nc",
               "e4oper.an.sfc.196404.moore543260.nc",
               "e4oper.an.sfc.196403.moore543260.nc",
               "e4oper.an.sfc.196402.moore543260.nc",
               "e4oper.an.sfc.196401.moore543260.nc",
               "e4oper.an.sfc.196312.moore543260.nc",
               "e4oper.an.sfc.196311.moore543260.nc",
               "e4oper.an.sfc.196310.moore543260.nc",
               "e4oper.an.sfc.196309.moore543260.nc",
               "e4oper.an.sfc.196308.moore543260.nc",
               "e4oper.an.sfc.196307.moore543260.nc",
               "e4oper.an.sfc.196306.moore543260.nc",
               "e4oper.an.sfc.196305.moore543260.nc",
               "e4oper.an.sfc.196304.moore543260.nc",
               "e4oper.an.sfc.196303.moore543260.nc",
               "e4oper.an.sfc.196302.moore543260.nc",
               "e4oper.an.sfc.196301.moore543260.nc",
               "e4oper.an.sfc.196212.moore543260.nc",
               "e4oper.an.sfc.196211.moore543260.nc",
               "e4oper.an.sfc.196210.moore543260.nc",
               "e4oper.an.sfc.196209.moore543260.nc",
               "e4oper.an.sfc.196208.moore543260.nc",
               "e4oper.an.sfc.196207.moore543260.nc",
               "e4oper.an.sfc.196206.moore543260.nc",
               "e4oper.an.sfc.196205.moore543260.nc",
               "e4oper.an.sfc.196204.moore543260.nc",
               "e4oper.an.sfc.196203.moore543260.nc",
               "e4oper.an.sfc.196202.moore543260.nc",
               "e4oper.an.sfc.196201.moore543260.nc",
               "e4oper.an.sfc.196112.moore543260.nc",
               "e4oper.an.sfc.196111.moore543260.nc",
               "e4oper.an.sfc.196110.moore543260.nc",
               "e4oper.an.sfc.196109.moore543260.nc",
               "e4oper.an.sfc.196108.moore543260.nc",
               "e4oper.an.sfc.196107.moore543260.nc",
               "e4oper.an.sfc.196106.moore543260.nc",
               "e4oper.an.sfc.196105.moore543260.nc",
               "e4oper.an.sfc.196104.moore543260.nc",
               "e4oper.an.sfc.196103.moore543260.nc",
               "e4oper.an.sfc.196102.moore543260.nc",
               "e4oper.an.sfc.196101.moore543260.nc",
               "e4oper.an.sfc.196012.moore543260.nc",
               "e4oper.an.sfc.196011.moore543260.nc",
               "e4oper.an.sfc.196010.moore543260.nc",
               "e4oper.an.sfc.196009.moore543260.nc",
               "e4oper.an.sfc.196008.moore543260.nc",
               "e4oper.an.sfc.196007.moore543260.nc",
               "e4oper.an.sfc.196006.moore543260.nc",
               "e4oper.an.sfc.196005.moore543260.nc",
               "e4oper.an.sfc.196004.moore543260.nc",
               "e4oper.an.sfc.196003.moore543260.nc",
               "e4oper.an.sfc.196002.moore543260.nc",
               "e4oper.an.sfc.196001.moore543260.nc",
               "e4oper.an.sfc.195912.moore543260.nc",
               "e4oper.an.sfc.195911.moore543260.nc",
               "e4oper.an.sfc.195910.moore543260.nc",
               "e4oper.an.sfc.195909.moore543260.nc",
               "e4oper.an.sfc.195908.moore543260.nc",
               "e4oper.an.sfc.195907.moore543260.nc",
               "e4oper.an.sfc.195906.moore543260.nc",
               "e4oper.an.sfc.195905.moore543260.nc",
               "e4oper.an.sfc.195904.moore543260.nc",
               "e4oper.an.sfc.195903.moore543260.nc",
               "e4oper.an.sfc.195902.moore543260.nc",
               "e4oper.an.sfc.195901.moore543260.nc",
               "e4oper.an.sfc.195812.moore543260.nc",
               "e4oper.an.sfc.195811.moore543260.nc",
               "e4oper.an.sfc.195810.moore543260.nc",
               "e4oper.an.sfc.195809.moore543260.nc",
               "e4oper.an.sfc.195808.moore543260.nc",
               "e4oper.an.sfc.195807.moore543260.nc",
               "e4oper.an.sfc.195806.moore543260.nc",
               "e4oper.an.sfc.195805.moore543260.nc",
               "e4oper.an.sfc.195804.moore543260.nc",
               "e4oper.an.sfc.195803.moore543260.nc",
               "e4oper.an.sfc.195802.moore543260.nc",
               "e4oper.an.sfc.195801.moore543260.nc",
               "e4oper.an.sfc.195712.moore543260.nc",
               "e4oper.an.sfc.195711.moore543260.nc",
               "e4oper.an.sfc.195710.moore543260.nc",
               "e4oper.an.sfc.195709.moore543260.nc"))
