## calculating spectral exponent over sliding windows of time
library(tidyverse)
library(broom)
library(evobiR)
library(sp)
library(raster)
select <- dplyr::select


##############################################
#####              FUNCTIONS:           ######
##############################################

## function to calculate spectral exponent over a time series window
spectral_exponent_calculator <- function(ts_window) {
  
  l <- length(ts_window)
  
  # Fourier transform the time series window: 
  dft <- fft(ts_window)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral <- data.frame(freq = freq, power = amp^2)
  
  # ## plot spectrum:
  # spectral %>% 
  #   ggplot(aes(x = freq, y = power)) + geom_line() +
  #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")
  
  ## get estimate of spectral exponent over time series window:
  model_output <- lm(spectral, formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  return(model_output$estimate)
}

## function to calculate spectral exponent over a time series within sliding windows of varying widths (from 5-10 years) with a time step of one year 
## takes as input:
##      - a detrended time series
##      - its latitude and longitude vectors 
## returns a list containing: 
##      - matrix of change in spectral exponents using each window width
##      - data frame of spectral exponents within each sliding window across all locations and for all window widths 
sliding_window_spec_exp <- function(path) {
 
  ## read in spatial chunk file names:
  filepath = paste(path, "sp_files.rds", sep = "")
  names = readRDS(filepath)
  
  l_filenames <- str_replace_all(names, "spatial_temps", 'l-detrended')
  s_filenames <-  str_replace_all(names, "spatial_temps", 's-detrended')
 
  ## create function for calculating date
  as_date <- function(x, origin = getOption("date_origin")){
    origin <- ifelse(is.null(origin), "1970-01-01", origin)
    as.Date(x, origin)
  }
  options(date_origin = "1869-12-31")
  
  long_index = 1
  lat_index = 1
  count = 1
  while(count < length(l_filenames)+1) {
    
    ## retrieve spatial chunk
    l_detrended_tas <- readRDS(l_filenames[count])
    s_detrended_tas <- readRDS(s_filenames[count])
    
    ## get lat and long bounds of chunk:
    long_bound1 <- long_index
    long_bound2 <- long_index + 59
    lat_bound1 <- lat_index
    lat_bound2 <- lat_index + 59
    
    spec_exp_list <- list()
    element <- 1
    x = 1 ## represents longitude index
    while (x < 61) {
      y = 1 ## represents latitude index
      while (y < 61) {
        l_local_ts <- l_detrended_tas[x, y, ] ## get the local detrended time series 
        s_local_ts <- s_detrended_tas[x, y, ]
        
        local_ts <- data.frame(time = 1:length(l_local_ts), ## add simple time column representing days from 1850-01-01
                                 l_temp = l_local_ts,
                                 s_temp = s_local_ts,
                                 date = as_date(1:length(l_local_ts))) %>% ## add a date column 
          mutate(year = str_split_fixed(.$date, pattern = "-", n = 2)[,1]) %>% ## add a year column
          group_by(year) ## group by year
        
        #########################################
        ##        SENSITIVITY ANALYSIS:        ##
        #########################################
        ## calculate spectral exponent using FFT over n year windows
        ## store spectral exponents and calculate slope
        n = 5
        while (n < 11) {
          year_start <- 1870
          year_stop <- 1870 + n
          
          while (year_start <= (2100 - n)) {
            
            ## extract temps within time window 
            ts_chunk <- filter(local_ts, year %in% year_start:year_stop)
    
            ## calcualte spectral exponent in window 
            l_exp <- spectral_exponent_calculator(ts_chunk$l_temp)
            s_exp <- spectral_exponent_calculator(ts_chunk$s_temp)
            
            ## store: 
            spec_exp_list[[element]] <- c(l_exp, s_exp, year_start, year_stop, lat[lat_index+y-1],
                                       long[long_index+x-1], paste(n, "years"))
            
            ## move to next window
            year_start = year_stop
            year_stop = year_stop + n
            
            element = element + 1
          }
          
          # # plot:
          # spec_exp_all[,1:6] <- sapply(spec_exp_all[,1:6], as.numeric) %>%
          #   ggplot(., aes (x = window_start_year, y = l_spec_exp)) + geom_point() +
          #   geom_smooth(method = "lm") + labs(x = "Time window start index",
          #                                     y = "Spectral exponent over window")
          # spec_exp_all[,1:6] <- sapply(spec_exp_all[,1:6], as.numeric) %>%
          #   ggplot(., aes (x = window_start_year, y = s_spec_exp)) + geom_point() +
          #   geom_smooth(method = "lm") + labs(x = "Time window start index",
          #                                     y = "Spectral exponent over window")
          
          ## move to next window width
          n = n + 1
        }
        
        print(paste("Calculating spectral exponent for x = ",x, ", y = ", y, sep = ""))
        y = y + 1
      }
      
      x = x + 1
    }
    
    ## bind rows in list into data frame
    spec_exp_df <- data.frame(do.call(rbind, spec_exp_list), stringsAsFactors = FALSE)
    colnames(spec_exp_df) <- c("l_spec_exp", "s_spec_exp", "window_start_year",
                                "window_stop_year", "lat", "long", "time_window_width")
    
    ## convert numbers to numeric 
    spec_exp_df[,1:6] <- sapply(spec_exp_df[,1:6], as.numeric) 
    
    ## regress spectral exponent and extract slope representing change in spectral exponent over time for each location and window width
    l_model_output <- spec_exp_df %>%
      group_by(lat, long, time_window_width) %>%
      do(tidy(lm(., formula = l_spec_exp ~ window_start_year))) %>%
      filter(term == "window_start_year")
    
    colnames(l_model_output)[5:8] <- paste("l", colnames(l_model_output)[5:8], sep = "_")
    
    s_model_output <- spec_exp_df %>%
        group_by(lat, long, time_window_width) %>%
        do(tidy(lm(., formula = s_spec_exp ~ window_start_year))) %>%
        filter(term == "window_start_year")
    
    colnames(s_model_output)[5:8] <- paste("s", colnames(s_model_output)[5:8], sep = "_")
    
    ## bind model output columns to spectral exponent data:
    spec_exp_df <- left_join(spec_exp_df, l_model_output) %>%
      left_join(., s_model_output)
    
    
    ## add filename to list:
    if(count == 1) {
      se_filenames <- paste("data-processed/", path, "/spec-exp_long-", 
                             long_bound1,"-", long_bound2, "_lat-", lat_bound1, "-", lat_bound2,
                             ".csv", sep = "")
    }
    else {
      se_filenames <- append(se_filenames,
                             paste("data-processed/", path, "/spec-exp_long-", 
                                   long_bound1,"-", long_bound2, "_lat-", lat_bound1, "-", lat_bound2,
                                   ".csv", sep = ""))
    }
    
    ## save spectral exponent for this chunk:
    write.csv(spec_exp_df, se_filenames[count], row.names = FALSE)
    
    ## advance lat and long indecies to move to next chunk
    if (count == 6) {
      lat_index <- lat_index + 60
      long_index <- 1
    }
    else if (count == 12) {
      lat_index <- lat_index + 60
      long_index <- 1
    }
    else {
      long_index <- long_index + 60
    }
    
    count = count + 1
  }
  
  
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
  
  spec_exp_filename <- paste("data-processed/large-files/", path, "_spec-exp.csv", sep = "")
  
  ## write out to one file: 
  write.csv(spec_exp, spec_exp_filename, row.names = FALSE)
  
  return(spec_exp_filename)
}


## function to convert output from sliding_window_spec_exp into 6 rasterStacks (one for each window width, from 5-10 years) that can be used with the gVoCC functions 
create_rasterStack <- function(spec_exp_filename) {
  
  spec_exp <- read.csv(spec_exp_filename)
  
  ## get path folder
  path <- str_split_fixed(spec_exp_filename, pattern = "spec-exp.csv", n = 2)[1,1]
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
              res = 1)
  
  ## reorder time_window_width factor so list elements are in order:
  spec_exp$time_window_width <- factor(spec_exp$time_window_width, levels = 
                                         c("5 years", "6 years", "7 years", "8 years",
                                           "9 years", "10 years"))
  
  ## split spectral exponent data by window width:
  ww_split <- split(spec_exp, spec_exp$time_window_width)
  
  ## for each window width:
  i = 1 
  while (i < length(ww_split) + 1) {
    
    ## split by window_start_year (each window step), loop through them and create raster layer for each
    ww <- ww_split[[i]]
    step_split <- split(ww, ww$window_start_year)
    step = length(step_split) ## loop backwards since rasterStacks add layers to beginning
    while (step >= 1) {
      
      l_sp_df <- step_split[[step]] %>%
        select(long, lat, l_spec_exp)
        
      s_sp_df <- step_split[[step]] %>%
        select(long, lat, s_spec_exp) 
      
      ## create raster layer:
      l_layer <- rasterFromXYZ(l_sp_df)
      s_layer <- rasterFromXYZ(s_sp_df)
      #plot(l_layer)
      #plot(s_layer)
      
      ## add to temporary rasterstack:
      if (step == length(step_split)) {
        l_temp_stack <- l_layer
        s_temp_stack <- s_layer
      }
      else {
        l_temp_stack <- addLayer(l_layer, l_temp_stack)
        s_temp_stack <- addLayer(s_layer, s_temp_stack)
      }
      
      ## move to nested for loop
      step = step - 1
    }
    
    names(l_temp_stack) <- paste("window", 1:nlayers(l_temp_stack), sep = "_")
    names(s_temp_stack) <- paste("window", 1:nlayers(s_temp_stack), sep = "_")
    
    ## save temporary raster stack: 
    if (i == 1) {
      l_stack_list <- list(l_temp_stack)
      s_stack_list <- list(s_temp_stack)
    }
    else {
      l_stack_list <- append(l_stack_list, l_temp_stack)
      s_stack_list <- append(s_stack_list, s_temp_stack)
    }
    
    ## move to next time window width
    i = i + 1
  }
  
  ## name the list items 
  names(l_stack_list) <- names(ww_split)
  names(s_stack_list) <- names(ww_split)
  
  ## save the rasterstack 
  saveRDS(l_stack_list, paste(path, "l_stack_list.rds", sep = ""))
  saveRDS(s_stack_list, paste(path, "s_stack_list.rds", sep = ""))
  
  stacks <- list(l_stack_list, s_stack_list)
  
  ## return the two lists of rasterStacks
  return(stacks)
}

##############################################
##############################################
##############################################


# ## example use:
# 
# ## filenames of GCM data:
# cmcc_cesm_historical <- c("tas_day_CMCC-CESM_historical_r1i1p1_18500101-18541231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_18550101-18591231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_18600101-18641231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_18650101-18691231.nc",
#                           "tas_day_CMCC-CESM_historical_r1i1p1_18700101-18741231.nc",
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
# ## filename vectors:
# l_filenames <- c("/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-1-60_lat-1-60.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-61-120_lat-1-60.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-121-180_lat-1-60.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-181-240_lat-1-60.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-241-300_lat-1-60.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-301-360_lat-1-60.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-1-60_lat-61-120.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-61-120_lat-61-120.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-121-180_lat-61-120.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-181-240_lat-61-120.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-241-300_lat-61-120.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-301-360_lat-61-120.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-1-60_lat-121-180.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-61-120_lat-121-180.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-121-180_lat-121-180.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-181-240_lat-121-180.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-241-300_lat-121-180.rds",
#                 "/Volumes/ADATA HV620/GCM-test/l-detrended-spatial-chunk_long-301-360_lat-121-180.rds")
# 
# s_filenames <- c("/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-1-60_lat-1-60.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-61-120_lat-1-60.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-121-180_lat-1-60.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-181-240_lat-1-60.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-241-300_lat-1-60.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-301-360_lat-1-60.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-1-60_lat-61-120.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-61-120_lat-61-120.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-121-180_lat-61-120.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-181-240_lat-61-120.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-241-300_lat-61-120.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-301-360_lat-61-120.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-1-60_lat-121-180.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-61-120_lat-121-180.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-121-180_lat-121-180.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-181-240_lat-121-180.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-241-300_lat-121-180.rds",
#                  "/Volumes/ADATA HV620/GCM-test/s-detrended-spatial-chunk_long-301-360_lat-121-180.rds")
# 
# file_list <- list(l_filenames, s_filenames)
# 
# se_filenames <- c("data-processed/spec-exp_long-1-60_lat-1-60.csv",
#                  "data-processed/spec-exp_long-61-120_lat-1-60.csv",
#                  "data-processed/spec-exp_long-121-180_lat-1-60.csv",
#                  "data-processed/spec-exp_long-181-240_lat-1-60.csv",
#                  "data-processed/spec-exp_long-241-300_lat-1-60.csv",
#                  "data-processed/spec-exp_long-301-360_lat-1-60.csv",
#                  "data-processed/spec-exp_long-1-60_lat-61-120.csv",
#                  "data-processed/spec-exp_long-61-120_lat-61-120.csv",
#                  "data-processed/spec-exp_long-121-180_lat-61-120.csv",
#                  "data-processed/spec-exp_long-181-240_lat-61-120.csv",
#                  "data-processed/spec-exp_long-241-300_lat-61-120.csv",
#                  "data-processed/spec-exp_long-301-360_lat-61-120.csv",
#                  "data-processed/spec-exp_long-1-60_lat-121-180.csv",
#                  "data-processed/spec-exp_long-61-120_lat-121-180.csv",
#                  "data-processed/spec-exp_long-121-180_lat-121-180.csv",
#                  "data-processed/spec-exp_long-181-240_lat-121-180.csv",
#                  "data-processed/spec-exp_long-241-300_lat-121-180.csv",
#                  "data-processed/spec-exp_long-301-360_lat-121-180.csv")
# 
# ## path to GCM files from my computer:
# path = "/Volumes/ADATA HV620/GCM-test/"












## garbage;
## sliding window version: 
# {
# while(count < length(filenames)+1) {
#   
#   ## retrieve spatial chunk
#   l_detrended_tas <- readRDS(l_filenames[count])
#   s_detrended_tas <- readRDS(s_filenames[count])
#   
#   ## get lat and long bounds of chunk:
#   long_bound1 <- long_index
#   long_bound2 <- long_index + 59
#   lat_bound1 <- lat_index
#   lat_bound2 <- lat_index + 59
#   
#   x = 1 ## represents longitude index
#   while (x < 61) {
#     y = 1 ## represents latitude index
#     while (y < 61) {
#       l_local_ts <- l_detrended_tas[x, y, ] ## get the local detrended time series 
#       s_local_ts <- s_detrended_tas[x, y, ]
#       
#       #########################################
#       ##        SENSITIVITY ANALYSIS:        ##
#       #########################################
#       ## calculate spectral exponent using FFT over n year windows
#       ## store spectral exponents and calculate slope
#       n = 5
#       while (n < 11) {
#         ## find width of window in days (365*n)
#         window_width <- round(365*n, digits = 0)
#         
#         ## calcualte spectral exponent in each window of n years, moving one window width at a time
#         ## note: how to deal with leap years here?
#         l_exp <- SlidingWindow(FUN = spectral_exponent_calculator, data = l_local_ts, step = 365, 
#                                window = window_width)
#         s_exp <- SlidingWindow(FUN = spectral_exponent_calculator, data = s_local_ts, step = 365, 
#                                window = window_width)
#         
#         lse <- l_exp %>%
#           as.data.frame() %>%
#           mutate(time_window_start = seq(from = 1, by = 365, 
#                                          to = (floor(length(l_local_ts)-window_width)))) %>%
#           rename(spec_exp = ".") %>%
#           mutate(lat = lat[lat_index+y-1]) %>%
#           mutate(long = long[long_index+x-1]) %>%
#           mutate(time_window_width = paste(n, "years")) %>%
#           mutate(time_step = "1y")
#         
#         sse <- s_exp %>%
#           as.data.frame() %>%
#           mutate(time_window_start = seq(from = 1, by = 365, 
#                                          to = (floor(length(s_local_ts)-window_width)))) %>%
#           rename(spec_exp = ".") %>%
#           mutate(lat = lat[y]) %>%
#           mutate(long = long[x]) %>%
#           mutate(time_window_width = paste(n, "years")) %>%
#           mutate(time_step = "1y")
#         
#         ## plot:
#         # lse %>%
#         #   ggplot(., aes (x = time_window_start, y = spec_exp)) + geom_point() +
#         #   geom_smooth(method = "lm") + labs(x = "Time window start index",
#         #                                     y = "Spectral exponent over window")
#         # sse %>%
#         #   ggplot(., aes (x = time_window_start, y = spec_exp)) + geom_point() +
#         #   geom_smooth(method = "lm") + labs(x = "Time window start index",
#         #                                     y = "Spectral exponent over window")
#         
#         ## regress spectral exponent and extract slope representing change in spectral exponent over time
#         l_model_output <- lm(lse, formula = spec_exp ~ time_window_start) %>%
#           tidy(.) %>%
#           filter(term == "time_window_start")
#         s_model_output <- lm(lse, formula = spec_exp ~ time_window_start) %>%
#           tidy(.) %>%
#           filter(term == "time_window_start")
#         
#         lse <- cbind(lse, l_model_output)
#         sse <- cbind(sse, s_model_output)
#         
#         ## store in database 
#         if(x == 1 & y == 1 & n == 5) {
#           all_lse <- lse 
#           all_sse <- sse 
#         }
#         else {
#           all_lse <- rbind(all_lse,lse) 
#           all_sse <- rbind(all_sse,sse)  
#         }
#         
#         n = n + 1
#       }
#       
#       print(paste("Calculating spectral exponent for x = ",x, ", y = ", y, sep = ""))
#       y = y + 1
#     }
#     
#     x = x + 1
#   }
#   
#   ## add name to list:
#   if(count == 1) {
#     lse_filenames <- paste("l-spec-exp_long-", 
#                            long_bound1,"-", long_bound2, "_lat-", lat_bound1, "-", lat_bound2,
#                            ".rds", sep = "")
#     sse_filenames <- paste("s-spec-exp_long-", 
#                            long_bound1,"-", long_bound2, "_lat-", lat_bound1, "-", lat_bound2,
#                            ".rds", sep = "")
#   }
#   else {
#     lse_filenames <- append(l_filenames,
#                             paste("l-spec-exp_long-", 
#                                   long_bound1,"-", long_bound2, "_lat-", lat_bound1, "-", lat_bound2,
#                                   ".rds", sep = ""))
#     sse_filenames <- append(s_filenames,
#                             paste("s-spec-exp_long-", 
#                                   long_bound1,"-", long_bound2, "_lat-", lat_bound1, "-", lat_bound2,
#                                   ".rds", sep = ""))
#   }
#   
#   ## save detrended spatial chunks:
#   saveRDS(lse, lse_filenames[count])
#   saveRDS(sse, sse_filenames[count])
#   
#   ## advance lat and long indecies to move to next chunk
#   if (count == 6) {
#     lat_index <- lat_index + 60
#     long_index <- 1
#   }
#   else if (count == 12) {
#     lat_index <- lat_index + 60
#     long_index <- 1
#   }
#   else {
#     long_index <- long_index + 60
#   }
#   
#   count = count + 1
# }
# 
# 
# 
# 
# x = 1
# while (x < length(long)+1) {
#   y = 1 
#   while (y < length(lat)+1) {
#     local_ts <- detrended_tas[x, y, ] ## get the local detrended time series 
#     
#     #########################################
#     ##        SENSITIVITY ANALYSIS:        ##
#     #########################################
#     ## calculate spectral exponent using FFT over n year windows
#     ## store spectral exponents and calculate slope
#     n = 5
#     while (n < 11) {
#       ## find width of window in days (365.25*n)
#       window_width <- round(365*n, digits = 0)
#       
#       ## calcualte spectral exponent in each window of n years, moving one window width at a time
#       exp <- SlidingWindow(FUN = spectral_exponent_calculator, data = local_ts, step = 365, 
#                            window = window_width)
#       
#       spec_exp <- exp %>%
#         as.data.frame() %>%
#         mutate(time_window_start = seq(from = 1, by = 365, 
#                                        to = (floor(length(local_ts)-window_width)))) %>%
#         rename(spec_exp = ".") %>%
#         mutate(lat = lat[y]) %>%
#         mutate(long = long[x]) %>%
#         mutate(time_window_width = paste(n, "years")) %>%
#         mutate(time_step = "1y")
#       
#       # ## plot:
#       # spec_exp %>%
#       #   ggplot(., aes (x = time_window_start, y = spec_exp)) + geom_point() +
#       #   geom_smooth(method = "lm") + labs(x = "Time window start index",
#       #                                     y = "Spectral exponent over window")
#       
#       ## regress spectral exponent and extract slope representing change in spectral exponent over time
#       model_output <- lm(spec_exp, formula = spec_exp ~ time_window_start) %>%
#         tidy(.) %>%
#         filter(term == "time_window_start")
#       
#       spec_exp <- cbind(spec_exp, model_output)
#       
#       ## store in database 
#       if(x == 1 & y == 1 & n == 5) {
#         all_spec_exp <- spec_exp 
#       }
#       else {
#         all_spec_exp <- rbind(all_spec_exp, spec_exp)
#       }
#       
#       ## store in array
#       spec3D[y, x, n-4] <- model_output$estimate
#       
#       ## advance to next time window width
#       print(paste("Calculating spectral exponent for x = ",x, ", y = ", y, ", window width = ", n, " years.", sep = ""))
#       n = n + 1
#     }
#     
#     
#     ## advance to next longitude
#     y = y + 1
#   }
#   ## advance to next latitude
#   x = x + 1
# }
# 
# dimnames(spec3D) <- NULL
# 
# return(list(spec3D, all_spec_exp))
# }
# 
# 