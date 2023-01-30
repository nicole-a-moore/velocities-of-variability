### analyzing berkeley earth historical observations of air temperature 
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
###                 filenames                  ## 
#################################################
path = "data-raw/BerkeleyEarth/"

filenames <- c("Complete_TAVG_Daily_LatLong1_1880.nc",
               "Complete_TAVG_Daily_LatLong1_1890.nc",
               "Complete_TAVG_Daily_LatLong1_1900.nc",
               "Complete_TAVG_Daily_LatLong1_1910.nc",
               "Complete_TAVG_Daily_LatLong1_1920.nc",
               "Complete_TAVG_Daily_LatLong1_1930.nc",
               "Complete_TAVG_Daily_LatLong1_1940.nc",
               "Complete_TAVG_Daily_LatLong1_1950.nc",
               "Complete_TAVG_Daily_LatLong1_1960.nc",
               "Complete_TAVG_Daily_LatLong1_1970.nc",
               "Complete_TAVG_Daily_LatLong1_1980.nc",
               "Complete_TAVG_Daily_LatLong1_1990.nc",
               "Complete_TAVG_Daily_LatLong1_2000.nc",
               "Complete_TAVG_Daily_LatLong1_2010.nc",
               "Complete_TAVG_Daily_LatLong1_2020.nc")

#################################################
###            resample observations           ## 
#################################################
r <- raster(ymx = 90, ymn = -90, xmn = 0, xmx = 360, res = 1) 
names_tas <- c()
i=1
while (i < length(filenames)+1) {
  # open file, get time variable
  file = paste(path, filenames[i], sep = "")
  nc <- nc_open(file)
  year <- ncvar_get(nc, "year")
  month <- ncvar_get(nc, "month")
  day <- ncvar_get(nc, "day")
  nc_close(nc)

  time <- paste(year, month, day, sep = "/")
  ## rasterize tas:
  tas <- stack(file, varname="temperature")
  clim <- stack(file, varname="climatology")

  extent(tas) <- c(0, 360, -90, 90)
  extent(clim) <- c(0, 360, -90, 90)

  if (any(str_detect(time, "2/29"))) {
    ly_days <- which(str_detect(time, "2/29"))
    ly_days <- ly_days[!ly_days %in% which(str_detect(time, "12/29"))]
    tas <- tas[[-ly_days]]
    time <- time[-ly_days]
  }
  
  ## add climatology to anomaly
  tas <- tas + clim

  #####      STANDARDIZE TO 1X1 DEGREE GRID   #####
  ## raster object of desired resolution/extent
  rtas <- resample(tas, r, method = 'bilinear',
                     filename = paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/resampled_tas_",
                                      filenames[i], sep = ""), overwrite = T)
  names_tas <- append(names_tas, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/resampled_tas_", filenames[i], 
                                       sep = ""))
  
  
  print(paste0("Done file number: ", i), stdout())
  i = i+1
}

saveRDS(names_tas, "/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/tas_filenames.rds")


## get paths and filenames 
files_tas <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/tas_filenames.rds")

## make date vector
## first date: 1880/1/1
## last date: 2021/12/31
dates <- paste(rep(seq(1880, 2021), each = 365), ".", rep(seq(1:365), 140), sep = "")

reorganize_GCM(filenames = files_tas)

## detrend and make spatial chunks 
reorganize_GCM <- function(filenames) {

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
          #####     REMOVING OUTLIERS    #####  
          ## temperature remove values > 60C (333.15K)
          local_ts$temp[which(local_ts$temp > 60)] <- NA
          
          #####     LINEARLY AND SEASONALLY DETREND TIME SERIES IN EACH RASTER CELL    #####
          ## create empty objects
          local_ts$md <- str_split_fixed(dates, "\\.", n=2)[,2]
          
          ts_df <- local_ts %>%
            group_by(md) %>%
            do(mutate(., temp_profile = mean(.$temp, na.rm=TRUE))) %>% ## compute temp climatology for each day of year
            ungroup() %>%
            mutate(s_detrended_temp = temp - temp_profile) %>% ## create column representing seasonally detrended
            arrange(., time)
          
          ## run linear regression for grid cell
          l_output <- lm(ts_df, formula = temp ~ time)
          s_output <- lm(ts_df, formula = s_detrended_temp ~ time)
          
          ## extract residuals and add to detrended temps objects:
          na_indecies <- which(is.na(local_ts$temp))
          na_indecies <- data.frame(num = local_ts$time) %>%
            mutate(is_na = ifelse(num %in% na_indecies, TRUE, FALSE))
          
          not_na <- filter(na_indecies, !is_na)
          not_na$l_temps <- l_output$residuals
          not_na$s_temps <- s_output$residuals
          
          detrended <- left_join(na_indecies, not_na, by = c("num", "is_na"))
          s_detrended = detrended$s_temps
          l_detrended = detrended$l_temps
          og_temps = local_ts$temp
          
          #####  INTERPOLATING MISSING TEMPS   #####
          ## interpolate if temps are missing
          if(length(which(is.na(s_detrended))) != 0 & length(which(is.na(s_detrended)))!= length(s_detrended)) {
            l_detrended <- imputeTS::na_kalman(l_detrended, smooth = TRUE, model = "StructTS")
            s_detrended <- imputeTS::na_kalman(s_detrended, smooth = TRUE, model = "StructTS")
            og_temps <- imputeTS::na_kalman(og_temps, smooth = TRUE, model = "StructTS")
            
            ## if first values are missing, get rid of poorly interpolated temps
            if(first(which(is.na(detrended$s_temps))) == 1) {
              l_detrended[first(which(is.na(detrended$s_temps))):(first(which(!is.na(detrended$s_temps)))-1)] = NA
              s_detrended[first(which(is.na(detrended$s_temps))):(first(which(!is.na(detrended$s_temps)))-1)] = NA
              og_temps[first(which(is.na(local_ts$temp))):(first(which(!is.na(local_ts$temp)))-1)] = NA
              
            }
          }
          
          ## save time series 
          temps_df[y,x,] <- og_temps
          l_detrended_temps[y,x,] <- l_detrended
          s_detrended_temps[y,x,] <- s_detrended 
          
          print(paste("Done detrending and interpolating x ", x,  " y ", y, " of chunk #", count,sep = ""))
        }
        y = y + 1
      }
      x = x + 1
    }
    
    ## save:
    sp_files[count] <- paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_spatial_temps_lon-", lon_bound1,"-", lon_bound2,
                             "_lat-", lat_bound1, "-", lat_bound2, ".nc", sep = "")
    ArrayToNc(temps_df, file_path = paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_spatial_temps_lon-", lon_bound1,"-", lon_bound2,
                                          "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(l_detrended_temps, file_path = paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_l-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(s_detrended_temps, file_path = paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_s-detrended_lon-", lon_bound1,"-", lon_bound2,
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
  
  saveRDS(sp_files, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_sp_files.rds", sep = ""))
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
  
  saveRDS(temp, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/temp_", count, ".rds", sep = ""))
  
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
saveRDS(mosaic, "/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/missing-data-count.rds")

## calculate spectral exponent after getting rid of sea in tas
sp_files_tas <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_sp_files.rds")

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

## function to calculate spectral exponent over a time series window using periodogram
spectral_exponent_calculator_PSD <- function(ts_window, l) {
  # Fourier transform the time series window: 
  dft <- fft(ts_window)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral <- data.frame(freq = freq, power = amp^2)
  
  # plot spectrum:
  # gg = spectral %>%
  #   filter(freq < 1/8*max(spectral$freq)) %>%
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

## function to calculate spectral exponent over a time series window using average wavelet coefficient method
spectral_exponent_calculator_AWC <- function(ts_window, N) {
  ## a. compute wavelet transform
  wavelets <- biwavelet::wt(data.frame(time = 1:N, val = ts_window), do.sig = F)
  
  ## b. calculate arithmetic mean with respect to the translation coefficient (b)
  data <- data.frame(avg_wavelet_coeff = rowMeans(sqrt(wavelets$power)), period = wavelets$period)
  
  # # plot average coefficients versus period on a log–log plot
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
sliding_window_spec_exp <- function(names) {

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
    x = 1 ## longitude
    while (x < ncol(l_detrended_tas)+1) {
      y = 1 ## latitude
      
      ## parallelize:
      spec_exp_list[[x]] <- foreach(y=1:60, combine=rbind)  %dopar% {
        l_local_ts <- l_detrended_tas[y,x,] ## get the local detrended time series
        s_local_ts <- s_detrended_tas[y,x,]
        
        ## if no time series, skip to next latitude
        if (length(which(is.na(l_local_ts))) == 51830) {
          rep(NA, 19)
        }
        else {
          local_ts <- data.frame(time = 1:length(l_local_ts), ## add integer time 
                                 l_temp = l_local_ts,
                                 s_temp = s_local_ts,
                                 date = dates) %>% ## add a date column
            mutate(year = str_split_fixed(.$date, 
                                          pattern = "\\.", n = 2)[,1]) %>% ## add a year column
            group_by(year) ## group by year
          
          ## calculate spectral exponent across all years
          ## get length of non- NA values
          L = length(which(!is.na(local_ts$l_temp)))
          first = first(which(!is.na(local_ts$l_temp)))

          ## preprocess the time series:
          ## a. subtracting mean
          ts_s_all <- local_ts$s_temp[first:nrow(local_ts)] - mean(local_ts$s_temp[first:nrow(local_ts)])
     
          ## b. windowing - multiply by a parabolic window
          window_s <- parabolic_window(series = ts_s_all, N = L)
          ts_s_all <- ts_s_all*window_s
      
          ## c. bridge detrending (endmatching)
          ## ie. subtracting from the data the line connecting the first and last points of the series
          ts_s_all <- bridge_detrender(windowed_series = ts_s_all, N = L)
       
          ## calculate spectral exponent in window using PSD and AWC methods
          s_exp_PSD_all <- spectral_exponent_calculator_PSD(ts_s_all, l = L)
      
          s_exp_PSD_low_all <- s_exp_PSD_all[[1]]
          s_exp_PSD_high_all <- s_exp_PSD_all[[2]]
          s_exp_PSD_all_all <- s_exp_PSD_all[[3]]
        
          #########################################
          ##        SENSITIVITY ANALYSIS:        ##
          #########################################
          ## calculate spectral exponent using FFT over n year windows
          ## store spectral exponents and calculate slope
          all_windows <- list()
          n = 5
          element <- 1
          while (n < 11) {
            year_start <- 1880
            year_stop <- 1880 + n - 1
            
            while (year_start <= (2020 - n)) {
              ## extract temps within time window
              ts_chunk <- filter(local_ts, year %in% year_start:year_stop)
              
              ## if there are no missing temps, calculate spectral exponent in window
              ## otherwise skip
              if(!any(is.na(ts_chunk$l_temp))) {
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
                  all_windows[[element]] <- c(l_exp_PSD_low, s_exp_PSD_low,
                                              l_exp_AWC, s_exp_AWC, 
                                              l_exp_PSD_high, s_exp_PSD_high, 
                                              l_exp_PSD_all, s_exp_PSD_all, 
                                              year_start, year_stop, 
                                              lat[y + lat_index],
                                              lon[x + lon_index], paste(n, "years"))
                } 
                else {
                  ## store:
                  all_windows[[element]] <- c(NA, s_exp_PSD_low, 
                                              NA, NA, 
                                              NA, s_exp_PSD_high, 
                                              NA, s_exp_PSD_all,
                                              year_start, year_stop, 
                                              lat[y + lat_index],
                                              lon[x + lon_index], paste(n, "years"))
                }
              }
              
              ## move to next window
              year_start = year_stop + 1
              year_stop = year_stop + n 
              
              element = element + 1
            }
            
            ## move to next window width
            n = n + 1
          }
          all_windows <- as.data.frame(do.call("rbind", all_windows))
          all_windows$s_exp_PSD_low_all <- rep(s_exp_PSD_low_all, nrow(all_windows))
          all_windows$s_exp_PSD_high_all <- rep(s_exp_PSD_high_all,  nrow(all_windows))
          all_windows$s_exp_PSD_all_all <- rep(s_exp_PSD_all_all,  nrow(all_windows))
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
    colnames(spec_exp_df) <- c("l_spec_exp_PSD_low", "s_spec_exp_PSD_low", 
                               "l_spec_exp_AWC", "s_spec_exp_AWC",
                               "l_spec_exp_PSD_high", "s_spec_exp_PSD_high", 
                               "l_spec_exp_PSD_all", "s_spec_exp_PSD_all", 
                               "window_start_year",
                               "window_stop_year", "lat", "lon", "time_window_width",
                               "s_spec_exp_PSD_low_all", "s_spec_exp_PSD_high_all",
                               "s_spec_exp_PSD_all_all")
    
    ## convert numbers to numeric
    spec_exp_df[,c(1:12, 14:16)] <- sapply(spec_exp_df[,c(1:12, 14:16)], as.numeric)
    
    ## regress spectral exponent and extract slope representing change in spectral exponent over time for each location and window width
    l_model_output_PSD_low <- spec_exp_df %>%
      filter(time_window_width == "10 years") %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = l_spec_exp_PSD_low ~ window_start_year, na.rm = TRUE))) %>%
      filter(term == "window_start_year") %>%
      select(-term)
    
    colnames(l_model_output_PSD_low)[4:7] <- paste("l", colnames(l_model_output_PSD_low)[4:7], "PSD_low", sep = "_")
    
    s_model_output_PSD_low <- spec_exp_df %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = s_spec_exp_PSD_low ~ window_start_year, na.rm = TRUE))) %>%
      filter(term == "window_start_year") %>%
      select(-term)
    
    colnames(s_model_output_PSD_low)[4:7] <- paste("s", colnames(s_model_output_PSD_low)[4:7], "PSD_low", sep = "_")
    
    l_model_output_PSD_high <- spec_exp_df %>%
      filter(time_window_width == "10 years") %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = l_spec_exp_PSD_high ~ window_start_year, na.rm = TRUE))) %>%
      filter(term == "window_start_year")  %>%
      select(-term)
    
    colnames(l_model_output_PSD_high)[4:7] <- paste("l", colnames(l_model_output_PSD_high)[4:7], "PSD_high", sep = "_")
    
    s_model_output_PSD_high <- spec_exp_df %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = s_spec_exp_PSD_high ~ window_start_year, na.rm = TRUE))) %>%
      filter(term == "window_start_year")  %>%
      select(-term)
    
    colnames(s_model_output_PSD_high)[4:7] <- paste("s", colnames(s_model_output_PSD_high)[4:7], "PSD_high", sep = "_")
    
    l_model_output_PSD_all <- spec_exp_df %>%
      filter(time_window_width == "10 years") %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = l_spec_exp_PSD_all ~ window_start_year, na.rm = TRUE))) %>%
      filter(term == "window_start_year")  %>%
      select(-term)
    
    colnames(l_model_output_PSD_all)[4:7] <- paste("l", colnames(l_model_output_PSD_all)[4:7], "PSD_all", sep = "_")
    
    s_model_output_PSD_all <- spec_exp_df %>%
      group_by(lat, lon, time_window_width) %>%
      do(tidy(lm(., formula = s_spec_exp_PSD_all ~ window_start_year, na.rm = TRUE))) %>%
      filter(term == "window_start_year")  %>%
      select(-term)
    
    colnames(s_model_output_PSD_all)[4:7] <- paste("s", colnames(s_model_output_PSD_all)[4:7], "PSD_all", sep = "_")
    
    l_model_output_AWC <- spec_exp_df %>%
      filter(time_window_width == "10 years") %>%
      group_by(lat, lon) %>%
      do(tidy(lm(., formula = l_spec_exp_AWC ~ window_start_year, na.rm = TRUE))) %>%
      filter(term == "window_start_year") %>%
      select(-term)
    
    colnames(l_model_output_AWC)[3:6] <- paste("l", colnames(l_model_output_AWC)[3:6], "AWC", sep = "_")
    l_model_output_AWC$time_window_width = "10 years"
    
    s_model_output_AWC <- spec_exp_df %>%
      group_by(lat, lon) %>%
      do(tidy(lm(., formula = s_spec_exp_AWC ~ window_start_year, na.rm = TRUE))) %>%
      filter(term == "window_start_year") %>%
      select(-term)
    
    colnames(s_model_output_AWC)[3:6] <- paste("s", colnames(s_model_output_AWC)[3:6], "AWC", sep = "_")
    s_model_output_AWC$time_window_width = "10 years"
    
    
    ## bind model output columns to spectral exponent data:
    spec_exp_df <- left_join(spec_exp_df, s_model_output_PSD_low) %>%
      left_join(., l_model_output_PSD_low) %>%
      left_join(., l_model_output_AWC) %>%
      left_join(., s_model_output_AWC) %>%
      left_join(., l_model_output_PSD_high) %>%
      left_join(., s_model_output_PSD_high) %>%
      left_join(., l_model_output_PSD_all) %>%
      left_join(., s_model_output_PSD_all)
    
    ## add filename to list:
    if(count == 1) {
      se_filenames <- paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_spec-exp_long-", 
                            lon_index,"-", lon_index + 60, "_lat-",
                            90-lat_index,"-", 90-lat_index-60,
                            ".csv", sep = "")
    }
    else {
      se_filenames <- append(se_filenames,
                             paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_spec-exp_long-", 
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
    
    spec_exp_df <- c()
    
    count = count + 1
  }
  
  saveRDS(se_filenames, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/se_filenames_temp.rds", sep = ""))
  
  return(se_filenames)
}

se_filenames_tas <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/se_filenames.rds")

## plot some:
data <- read.csv(se_filenames_tas[[10]]) %>%
  mutate(lat_lon = paste(lat, lon, sep = "_")) %>%
  filter(time_window_width == "5 years") %>%
  filter(lat_lon %in% unique(lat_lon)[c(6,32,82,56,99,101,85,23,66,
                                        8,302,45,74,80,232,145,2,4,55,66,7)])

######################################################################################################################
#####    make data frame that has lat, lon, first ts noise colour, whole ts noise colour, and change in colour  ######
######################################################################################################################
se_filenames = se_filenames_tas
file = 1
while (file <= length(se_filenames)) {
  
  if(file.exists(se_filenames[[file]])) {
    if(file == 1){
      data <- read.csv(se_filenames[[file]]) %>%
        filter(window_start_year %in% c(1880, 2015)) %>%
        select(lat, lon, time_window_width, window_start_year,
               s_spec_exp_PSD_high, s_spec_exp_PSD_low, s_spec_exp_PSD_all,
               s_spec_exp_PSD_high_all, s_spec_exp_PSD_low_all, s_spec_exp_PSD_all_all,
               s_estimate_PSD_high, s_estimate_PSD_low, s_estimate_PSD_all) %>%
        filter(time_window_width == "5 years") %>%
        unique() 
    }
    else {
      data <- read.csv(se_filenames[[file]]) %>%
        filter(window_start_year %in% c(1880, 2015)) %>%
        select(lat, lon, time_window_width, window_start_year,
               s_spec_exp_PSD_high, s_spec_exp_PSD_low, s_spec_exp_PSD_all,
               s_spec_exp_PSD_high_all, s_spec_exp_PSD_low_all, s_spec_exp_PSD_all_all,
               s_estimate_PSD_high, s_estimate_PSD_low, s_estimate_PSD_all) %>%
        filter(time_window_width == "5 years") %>%
        unique() %>%
        rbind(., data)
    }
  }
  
  file = file + 1
}

saveRDS(data, "/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_noise-colour.rds")

data <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarth_noise-colour.rds")

## find mean spectral colour across terrestrial environments in 2015 
data %>%
  filter(window_start_year == "2015") %>%
  select(lat, lon, s_spec_exp_PSD_low) %>%
  unique()%>%
  summarise(mean(s_spec_exp_PSD_low))

## 0.8284292

## find 95% quantiles
data %>%
  filter(window_start_year == "2015") %>%
  select(lat, lon, s_spec_exp_PSD_low) %>%
  unique()%>%
  ggplot(aes(x = s_spec_exp_PSD_low)) + geom_histogram()

data %>%
  filter(window_start_year == "2015") %>%
  select(lat, lon, s_spec_exp_PSD_low) %>%
  unique()%>%
  summarise(quantile(s_spec_exp_PSD_low, 0.05),
            quantile(s_spec_exp_PSD_low, 0.95))

## min = 0.07505786
## max = 1.506056 

## find min and max spectral change across terrestrial environments 
data %>%
 select(lat, lon, s_estimate_PSD_low) %>%
  unique()%>%
  summarise(min(s_estimate_PSD_low),
            max(s_estimate_PSD_low))

## min = -0.03194815 per 140 years = -0.1250417/200000 days 
## max = 0.0420622 per 140 years = 0.164627/200000 days 



##################################
#####    probably garbage   ######
##################################

####################################################################
#####    transform spectral exponent data into a rasterStack  ######
####################################################################
## function to convert output from script 04 into 6 rasterStacks (one for each sliding window width, from 5-10 years) that can be used with the gVoCC functions 
create_rasterStack <- function(se_filenames) {
  
  ## combine all spectral exponent csvs into one big dataframe
  file = 1
  while (file < length(se_filenames) + 1) {
    if(file.exists(se_filenames[file])) { 
      if (file == 1) {
        spec_exp <- read.csv(se_filenames[file])
      }
      else {
        spec_exp <- rbind(spec_exp, read.csv(se_filenames[file]))
      }
    }
    
    print(paste("Reading file #", file, "/", length(se_filenames), sep = ""))
    file = file + 1
  }
 
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
              res = 1)
  
  ## reorder time_window_width so list elements are in order of increasing time window widths:
  spec_exp$time_window_width <- factor(spec_exp$time_window_width, levels = 
                                         c("5 years", "6 years", "7 years", "8 years",
                                           "9 years", "10 years"))
  
  ## get rid of ones with large chunks of ts missing
  mosaic <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/missing-data-count.rds")
  ## make a mask to get rid of ones with more than 10000 missing values
  gt_10000 <- mosaic
  gt_10000[mosaic >= 10000] <- 1
  gt_10000[mosaic < 10000] <- NA
  plot(gt_10000)
  
  ## get lats and lons 
  gt_10000 <- data.frame(rasterToPoints(gt_10000))
  
  ## filter spec exp to exclude cells in gt_10000
  spec_exp <- mutate(spec_exp, lat_lon = paste(lat, lon))
  gt_10000 <- mutate(gt_10000, lat_lon = paste(y, x))
  
  spec_exp <- filter(spec_exp, !lat_lon %in% gt_10000$lat_lon)
  
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
        ## make sure extents are the same
        ## (won't be for layers with time windows that had all data missing)
        l_layer_PSD_low = extend(l_layer_PSD_low, l_temp_stack_PSD_low)
        s_layer_PSD_low = extend(s_layer_PSD_low, s_temp_stack_PSD_low)
        l_layer_AWC = extend(l_layer_AWC, l_temp_stack_AWC)
        s_layer_AWC = extend(s_layer_AWC, s_temp_stack_AWC)
        l_layer_PSD_high = extend(l_layer_PSD_high, l_temp_stack_PSD_high)
        s_layer_PSD_high = extend(s_layer_PSD_high, s_temp_stack_PSD_high)
        l_layer_PSD_all = extend(l_layer_PSD_all, l_temp_stack_PSD_all)
        s_layer_PSD_all = extend(s_layer_PSD_all, s_temp_stack_PSD_all)

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
    names(s_stack_list_PSD_high) <- names(l_stack_list_AWC) <- names(s_stack_list_AWC) <-
    names(l_stack_list_PSD_all) <- names(s_stack_list_PSD_all) <-
    names(ww_split)
  
  ## save the rasterstack 
  saveRDS(l_stack_list_PSD_low, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarthl_stack_list_PSD_low.rds", sep = ""))
  saveRDS(s_stack_list_PSD_low, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarths_stack_list_PSD_low.rds", sep = ""))
  saveRDS(l_stack_list_AWC, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarthl_stack_list_AWC.rds", sep = ""))
  saveRDS(s_stack_list_AWC, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarths_stack_list_AWC.rds", sep = ""))
  saveRDS(l_stack_list_PSD_high, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarthl_stack_list_PSD_high.rds", sep = ""))
  saveRDS(s_stack_list_PSD_high, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarths_stack_list_PSD_high.rds", sep = ""))
  saveRDS(l_stack_list_PSD_all, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarthl_stack_list_PSD_all.rds", sep = ""))
  saveRDS(s_stack_list_PSD_all, paste("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/BerkeleyEarths_stack_list_PSD_all.rds", sep = ""))

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

spec_exp %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", c(l_spec_exp_PSD_high, l_spec_exp_PSD_all, l_spec_exp_PSD_low,
                                                 l_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  group_by(window_start_year, exp_type) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window
  select(window_start_year, mean_spec_exp) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") +
  geom_smooth(method = "lm")

spec_exp %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", c(l_spec_exp_PSD_high, l_spec_exp_PSD_all, l_spec_exp_PSD_low,
                                                 l_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  mutate(lat_lon = paste(lat, lon)) %>%
  filter(lat_lon %in% unique(.$lat_lon)[c(1,20990,2000,8000,10000)]) %>%
  ggplot(., aes(x = window_start_year, y = spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spectral exponent") +
  geom_smooth(method = "lm", se = F) + facet_wrap(~lat_lon)

spec_exp %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", c(s_spec_exp_PSD_high, s_spec_exp_PSD_all, s_spec_exp_PSD_low,
                                                 s_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  mutate(lat_lon = paste(lat, lon)) %>%
  filter(lat_lon %in% unique(.$lat_lon)[c(1,20990,2000,8000,10000)]) %>%
  ggplot(., aes(x = window_start_year, y = spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spectral exponent") +
  geom_smooth(method = "lm", se = F) + facet_wrap(~lat_lon)

## lat -83.5, lon 161.5
## count = 5,  x = 42, y = 54
## example of ts with first half missing

## look at missingness
mosaic <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/missing-data-count.rds")
plot(mosaic)
hist(values(mosaic))
lots <- mosaic
lots[mosaic == 51830] <- NA
lots[mosaic < 10000 & mosaic != 51830] <- 0
lots[mosaic >= 10000 & mosaic != 51830] <- 1
plot(lots)

little <- mosaic
little[mosaic > 10000] <- 51830
little[little == 51830] <- NA

plot(little)

## make a mask to get rid of ones with more than 10000 missing values
gt_10000 <- mosaic
gt_10000[mosaic >= 10000] <- 1
gt_10000[mosaic < 10000] <- NA
plot(gt_10000)

## get lats and lons 
gt_10000 <- data.frame(rasterToPoints(gt_10000))

## filter spec exp to exclude cells in gt_10000
spec_exp <- mutate(spec_exp, lat_lon = paste(lat, lon))
gt_10000 <- mutate(gt_10000, lat_lon = paste(y, x))

sub_spec_exp <- filter(spec_exp, !lat_lon %in% gt_10000$lat_lon)

sub_spec_exp %>%
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

sub_spec_exp %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", c(l_spec_exp_PSD_high, l_spec_exp_PSD_all, l_spec_exp_PSD_low,
                                                 l_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  group_by(window_start_year, exp_type) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window
  select(window_start_year, mean_spec_exp) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") +
  geom_smooth(method = "lm")

sub_spec_exp %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", c(s_spec_exp_PSD_high, s_spec_exp_PSD_all, s_spec_exp_PSD_low,
                                                 s_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  mutate(lat_lon = paste(lat, lon)) %>%
  filter(lat_lon %in% unique(.$lat_lon)[c(1,40,2000,900,10000)]) %>%
  ggplot(., aes(x = window_start_year, y = spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spectral exponent") +
  geom_smooth(method = "lm", se = F) + facet_wrap(~lat_lon)

write.csv(sub_spec_exp, "/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/spec_exp_qc.csv", row.names = F)


## example plots of frankenstein series
first_spec_exp <- data %>%
  select(lat, lon, time_window_width, s_spec_exp_PSD_low, window_start_year) %>%
  filter(time_window_width == "5 years") %>%
  group_by(lat, lon) %>%
  filter(window_start_year == min(window_start_year)) %>%
  ungroup() %>%
  select(lat, lon, s_spec_exp_PSD_low) %>%
  unique() %>%
  rename("first_spec_exp" = s_spec_exp_PSD_low)

new <- data %>%
  filter(time_window_width == "5 years") %>%
  select(lat, lon, s_spec_exp_PSD_low, window_start_year) %>%
  left_join(., first_spec_exp) 

new %>%
  mutate(lat_lon = paste(lat, lon)) %>%
  filter(lat_lon %in% unique(.$lat_lon)[sample(1:907, 20)]) %>%
  ggplot(aes(x = window_start_year, y = s_spec_exp_PSD_low)) + geom_point(colour = "red", size = 0.5) +
  facet_wrap(~lat_lon) +
  geom_point(aes(y = first_spec_exp), colour = "lightblue", size = 0.5) +
  labs(x = "Year", y = "Spectral exponent (low freq)")

new %>%
  mutate(lat_lon = paste(lat, lon)) %>%
  filter(lat_lon %in% unique(.$lat_lon)[1]) %>%
  ggplot(aes(x = window_start_year, y = s_spec_exp_PSD_low)) + geom_point(colour = "red") +
  facet_wrap(~lat_lon)  +
  labs(x = "Year", y = "Spectral exponent (low freq)")



