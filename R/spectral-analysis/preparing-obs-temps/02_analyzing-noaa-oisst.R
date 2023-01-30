### analyzing noaa oisst reanalysis of sea surface temperature 
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
names_sst <- c()
i=1
while (i < length(filenames)+1) {
  # open file, get time variable
  file = paste(path, filenames[i], sep = "")
  nc <- nc_open(file)
  time <- ncvar_get(nc, "time") ## days since 1800-01-01
  nc_close(nc)

  ## calculate date in year, month, day format
  time = as.Date(time, origin=as.Date("1800-01-01"))
  time = str_replace_all(time, "\\-", "/")

  ## rasterize sst:
  sst <- stack(file, varname="sst")
  extent(sst) <- c(0, 360, -90, 90)

  if (any(str_detect(time, "2/29"))) {
    ly_days <- which(str_detect(time, "2/29"))
    ly_days <- ly_days[!ly_days %in% which(str_detect(time, "12/29"))]
    if(length(ly_days) !=0) {
      sst <- sst[[-ly_days]]
      time <- time[-ly_days]
    }
  }

  #####      STANDARDIZE TO 1X1 DEGREE GRID   #####
  ## raster object of desired resolution/extent
  rsst <- resample(sst, r, method = 'bilinear',
                     filename = paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/resampled_sst_",
                                      filenames[i], sep = ""), overwrite = T)
  names_sst <- append(names_sst, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/resampled_sst_", filenames[i], 
                                       sep = ""))
  
  
  print(paste0("Done file number: ", i), stdout())
  i = i+1
}

saveRDS(names_sst, "/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/resamp_filenames.rds")


## get paths and filenames 
files_sst <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/resamp_filenames.rds")

## make date vector
## first date: 1981/09/01
## last date: 2022/3/13
dates <- paste(rep(seq(1981, 2022), each = 365), ".", rep(seq(1:365), 41), sep = "")
dates <- dates[244:(41*365+72)]


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
          
          #####  INTERPOLATING MISSING TEMPS   #####
          ## interpolate if temps are missing
          if(length(which(is.na(s_detrended))) != 0 & length(which(is.na(s_detrended)))!= length(s_detrended)) {
            l_detrended <- imputeTS::na_kalman(l_detrended, smooth = TRUE, model = "StructTS")
            s_detrended <- imputeTS::na_kalman(s_detrended, smooth = TRUE, model = "StructTS")
            
            ## if first values are missing, get rid of poorly interpolated temps
            if(first(which(is.na(detrended$s_temps))) == 1) {
              l_detrended[first(which(is.na(detrended$l_temps))):first(which(!is.na(detrended$l_temps)))-1] = NA
              s_detrended[first(which(is.na(detrended$s_temps))):first(which(!is.na(detrended$s_temps)))-1] = NA
            }
          }
          l_detrended_temps[y,x,] <- l_detrended
          s_detrended_temps[y,x,] <- s_detrended ## save time series 
          
          print(paste("Done detrending and interpolating x ", x,  " y ", y, " of chunk #", count,sep = ""))
        }
        y = y + 1
      }
      x = x + 1
    }
    
    ## save:
    sp_files[count] <- paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_spatial_temps_lon-", lon_bound1,"-", lon_bound2,
                             "_lat-", lat_bound1, "-", lat_bound2, ".nc", sep = "")
    ArrayToNc(temps_df, file_path = paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_spatial_temps_lon-", lon_bound1,"-", lon_bound2,
                                          "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(l_detrended_temps, file_path = paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_l-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(s_detrended_temps, file_path = paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_s-detrended_lon-", lon_bound1,"-", lon_bound2,
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
  
  saveRDS(sp_files, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_sp_files.rds", sep = ""))
  ## returns list of spatial chunk filenames 
  return(sp_files)
}

reorganize_GCM(filenames = files_sst)


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
  
  saveRDS(temp, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/temp_", count, ".rds", sep = ""))
  
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
saveRDS(mosaic, "/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/missing-data-count.rds")

## calculate spectral exponent after getting rid of sea in sst
sp_files_sst <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_sp_files.rds")

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

## function to calculate spectral exponent over a time series window using average wavelet coefficient method
spectral_exponent_calculator_AWC <- function(ts_window, N) {
  ## a. compute wavelet transform
  wavelets <- biwavelet::wt(data.frame(time = 1:N, val = ts_window), do.sig = F)
  
  ## b. calculate arithmetic mean with respect to the translation coefficient (b)
  data <- data.frame(avg_wavelet_coeff = rowMeans(sqrt(wavelets$power), na.rm = TRUE), period = wavelets$period)
  
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
  
  lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
  lon <-seq(from = 0.5, to = 359.5, length.out = 360) 
  
  lon_index = 0
  lat_index = 0
  count = 1
  while(count < length(l_filenames)+1) {
    
    ## retrieve spatial chunk from nc file
    l_open = nc_open(l_filenames[count])
    l_detrended_sst = ncvar_get(l_open, "var1_1")
    nc_close(l_open)
    
    s_open = nc_open(s_filenames[count])
    s_detrended_sst = ncvar_get(s_open, "var1_1")
    nc_close(s_open)
    
    spec_exp_list <- list()
    x = 1 ## longitude
    while (x < ncol(l_detrended_sst)+1) {
      y = 1 ## latitude
      
      ## parallelize:
      spec_exp_list[[x]] <-  foreach(y=1:nrow(l_detrended_sst), combine=rbind)  %dopar% {
        l_local_ts <- l_detrended_sst[y,x,] ## get the local detrended time series
        s_local_ts <- s_detrended_sst[y,x,]
        
        ## if no time series, skip to next latitude
        if (length(which(is.na(l_local_ts))) == 14794) {
          rep(NA, 13)
        }
        else {
          local_ts <- data.frame(time = 1:length(l_local_ts), ## add integer time 
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
          element <- 1
          while (n < 11) {
            year_start <- 1981
            year_stop <- 1981 + n - 1
            
            while (year_start <= (2022 - n)) {
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
      se_filenames <- paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_spec-exp_long-", 
                            lon_index,"-", lon_index + 60, "_lat-",
                            90-lat_index,"-", 90-lat_index-60,
                            ".csv", sep = "")
    }
    else {
      se_filenames <- append(se_filenames,
                             paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_spec-exp_long-", 
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
  
  saveRDS(se_filenames, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/se_filenames.rds", sep = ""))
  
  return(se_filenames)
}

se_filenames_sst <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/se_filenames.rds")

######################################################################################################################
#####    make data frame that has lat, lon, first ts noise colour, whole ts noise colour, and change in colour  ######
######################################################################################################################
se_filenames = se_filenames_sst
file = 1
while (file <= length(se_filenames)) {
  
  if(file == 1){
    data <- read.csv(se_filenames[[file]]) %>%
      filter(window_start_year %in% c(1981, 2016)) %>%
      select(lat, lon, time_window_width, window_start_year,
             s_spec_exp_PSD_high, s_spec_exp_PSD_low, s_spec_exp_PSD_all,
             s_estimate_PSD_high, s_estimate_PSD_low, s_estimate_PSD_all) %>%
      filter(time_window_width == "5 years") %>%
      unique() 
  }
  else {
    data <- read.csv(se_filenames[[file]]) %>%
      filter(window_start_year %in% c(1981, 2016)) %>%
      select(lat, lon, time_window_width, window_start_year,
             s_spec_exp_PSD_high, s_spec_exp_PSD_low, s_spec_exp_PSD_all,
             s_estimate_PSD_high, s_estimate_PSD_low, s_estimate_PSD_all) %>%
      filter(time_window_width == "5 years") %>%
      unique() %>%
      rbind(., data)
  }
  
  file = file + 1
}

saveRDS(data, "/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_noise-colour.rds")

data <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_noise-colour.rds")

## find mean spectral colour across marine environments in 2016 
data %>%
  filter(window_start_year == "2016") %>%
  select(lat, lon, s_spec_exp_PSD_low) %>%
  unique()%>%
  summarise(mean(s_spec_exp_PSD_low))

## 1.335845

## find 95% quantiles
data %>%
  filter(window_start_year == "2016") %>%
  select(lat, lon, s_spec_exp_PSD_low) %>%
  unique()%>%
  ggplot(aes(x = s_spec_exp_PSD_low)) + geom_histogram()

data %>%
  filter(window_start_year == "2016") %>%
  select(lat, lon, s_spec_exp_PSD_low) %>%
  unique()%>%
  summarise(quantile(s_spec_exp_PSD_low, 0.05),
            quantile(s_spec_exp_PSD_low, 0.95))

## min = 0.8837043
## max = 1.91154

## find min and max spectral change across marine environments 
data %>%
  select(lat, lon, s_estimate_PSD_low) %>%
  unique()%>%
  summarise(min(s_estimate_PSD_low),
            max(s_estimate_PSD_low))

## min = -0.03409852 per 40 years = -0.467103/200000 days 
## max = 0.04619845 per 40 years = 0.6328555/200000 days 



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
    if (file == 1) {
      spec_exp <- read.csv(se_filenames[file])
    }
    else {
      spec_exp <- rbind(spec_exp, read.csv(se_filenames[file]))
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
  mosaic <- readRDS("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/missing-data-count.rds")
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
  saveRDS(l_stack_list_PSD_low, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_l_stack_list_PSD_low.rds", sep = ""))
  saveRDS(s_stack_list_PSD_low, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_s_stack_list_PSD_low.rds", sep = ""))
  saveRDS(l_stack_list_AWC, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_l_stack_list_AWC.rds", sep = ""))
  saveRDS(s_stack_list_AWC, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_s_stack_list_AWC.rds", sep = ""))
  saveRDS(l_stack_list_PSD_high, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_l_stack_list_PSD_high.rds", sep = ""))
  saveRDS(s_stack_list_PSD_high, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_s_stack_list_PSD_high.rds", sep = ""))
  saveRDS(l_stack_list_PSD_all, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_l_stack_list_PSD_all.rds", sep = ""))
  saveRDS(s_stack_list_PSD_all, paste("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/NOAA-OISST_s_stack_list_PSD_all.rds", sep = ""))

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

write.csv(spec_exp, "/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/spec_exp_qc.csv", row.names = F)

### comparing resutls gcm to results from observations
spec_exp_noaa <- read.csv("/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/spec_exp_qc.csv")
spec_exp_noaa$dataset <- "NOAA-OISST"

spec_exp_berk <- read.csv("data-processed/BerkeleyEarth/spec_exp_qc.csv")
spec_exp_berk$dataset = "Berkley Earth"

spec_exp_gcm <- read.csv("data-processed/GCM_spec_exp.csv")
spec_exp_gcm$dataset = "GCM"

## get rid of air temperature for places in the ocean:
spec_exp_gcm$lat_lon <- paste(spec_exp_gcm$lat, spec_exp_gcm$lon)

spec_exp_gcm <- spec_exp_gcm %>%
  filter(lat_lon %in% spec_exp_berk$lat_lon)

## rasterize and plot to make sure it worked:
spec_exp_gcm %>%
  filter(time_window_width == "5 years", window_start_year == "1871") %>%
  unique() %>%
  select(lon, lat, l_spec_exp_PSD_low) %>%
  rasterFromXYZ() %>%
  plot()

spec_exp_berk = select(spec_exp_berk, -lat_lon)
spec_exp_gcm = select(spec_exp_gcm, -lat_lon)
spec_exp_noaa = select(spec_exp_noaa, -lat_lon)

spec_exp_berk <- spec_exp_berk[-which(!colnames(spec_exp_berk) %in% colnames(spec_exp_gcm))]
spec_exp_noaa <- spec_exp_noaa[-which(!colnames(spec_exp_noaa) %in% colnames(spec_exp_gcm))]

spec = rbind(spec_exp_berk, spec_exp_noaa) %>%
  rbind(., spec_exp_gcm)

#### ADD PSD_ALL LATER

## add in the era 40 data
spec_exp_era_sst <- read.csv("data-processed/ERA-40/spec_exp_qc_sst.csv")
spec_exp_era_sst$dataset = "ERA-40 sst"

spec_exp_era_tas <- read.csv("data-processed/ERA-40/spec_exp_qc_tas.csv")
spec_exp_era_tas$dataset = "ERA-40 tas"

spec_exp_era_tas <- spec_exp_era_tas[-which(!colnames(spec_exp_era_tas) %in% colnames(spec_exp_gcm))]
spec_exp_era_sst <- spec_exp_era_sst[-which(!colnames(spec_exp_era_sst) %in% colnames(spec_exp_gcm))]

spec = rbind(spec_exp_era_sst, spec) %>%
  rbind(., spec_exp_era_tas) 

spec <- spec %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", 
         c(s_spec_exp_PSD_high, s_spec_exp_PSD_low,s_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) 

## only berk, gcm, noaa
spec %>%
  filter(!dataset %in% c("ERA-40 sst", "ERA-40 tas")) %>%
  group_by(window_start_year, exp_type, dataset) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_start_year, mean_spec_exp, dataset) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = dataset, shape = dataset)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent", shape = "", colour = "") + 
  geom_smooth(method = "lm") +
  facet_wrap(~exp_type)

## all
spec %>%
  group_by(window_start_year, exp_type, dataset) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_start_year, mean_spec_exp, dataset) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = dataset, shape = dataset)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent", shape = "", colour = "") + 
  geom_smooth(method = "lm") +
  facet_wrap(~exp_type)




#################################################
###                 filenames                  ## 
#################################################
path = "data-raw/NOAA-OISST/"

filenames <- c("sst.day.mean.1981.nc",
               "sst.day.mean.1982.nc",
               "sst.day.mean.1983.nc",
               "sst.day.mean.1984.nc",
               "sst.day.mean.1985.nc",
               "sst.day.mean.1986.nc",
               "sst.day.mean.1987.nc",
               "sst.day.mean.1988.nc",
               "sst.day.mean.1989.nc",
               "sst.day.mean.1990.nc",
               "sst.day.mean.1991.nc",
               "sst.day.mean.1992.nc",
               "sst.day.mean.1993.nc",
               "sst.day.mean.1994.nc",
               "sst.day.mean.1995.nc",
               "sst.day.mean.1996.nc",
               "sst.day.mean.1997.nc",
               "sst.day.mean.1998.nc",
               "sst.day.mean.1999.nc",
               "sst.day.mean.2000.nc",
               "sst.day.mean.2001.nc",
               "sst.day.mean.2002.nc",
               "sst.day.mean.2003.nc",
               "sst.day.mean.2004.nc",
               "sst.day.mean.2005.nc",
               "sst.day.mean.2006.nc",
               "sst.day.mean.2007.nc",
               "sst.day.mean.2008.nc",
               "sst.day.mean.2009.nc",
               "sst.day.mean.2010.nc",
               "sst.day.mean.2011.nc",
               "sst.day.mean.2012.nc",
               "sst.day.mean.2013.nc",
               "sst.day.mean.2014.nc",
               "sst.day.mean.2015.nc",
               "sst.day.mean.2016.nc",
               "sst.day.mean.2017.nc",
               "sst.day.mean.2018.nc",
               "sst.day.mean.2019.nc",
               "sst.day.mean.2020.nc",
               "sst.day.mean.2021.nc",
               "sst.day.mean.2022.nc")
