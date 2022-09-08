## trying inverse fourier transform 


## 1. read in un-detrended Berkeley Earth temps
## 2. Fourier transform them
## 3. inverse Fourier transform them and save results 

filenames <- readRDS("data-processed/BerkeleyEarth/be_sp_files.rds")

## read in the ocean coordinates
ocean <- read.csv("data-processed/masks/cmip5-ocean-coords.csv")
ocean_coords <- paste(ocean$lat, ocean$lon)

lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 

lon_index = 0
lat_index = 0
count = 1
while(count < length(filenames)+1) {
  
  ## retrieve spatial chunk from nc file
  open = nc_open(filenames[count])
  detrended_tas = ncvar_get(open, "var1_1")
  nc_close(open)
  
  x = 1 ## longitude
  while (x < ncol(detrended_tas)+1) {
    y = 1 ## latitude
    
    while(y < nrow(detrended_tas)+1) {
      local_ts <- detrended_tas[y,x,] ## get the local time series
      
      ## if no time series, skip to next
      if (length(which(is.na(local_ts))) == 51830) {
        y = y + 1
      }
      else {
        
        ## Fourier transform the raw temperatures 
        dft <- fft(local_ts)/l
        amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
        amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
        freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
        
        ## create periodogram data by squaring amplitude of FFT output
        spectral <- data.frame(freq = freq, power = amp^2)
        
        #plot spectrum:
        spectral %>%
          filter(freq < 1/8*max(spectral$freq)) %>%
          ggplot(aes(x = freq, y = power)) + geom_line() +
          scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")
        
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
      
        y = y + 1
      }
    }
    
    print(paste("Recreating time series for lat = ", lat[y + lat_index], 
                " and lon = ", lon[x + lon_index], sep = ""))
    
    x = x + 1
  }
  count = count + 1
}


## see if change in spectral exponent differs 


true_ts = local_ts

## analyze change 
local_ts <- data.frame(time = 1:length(local_ts), ## add integer time 
                       true_ts = true_ts,
                       date = dates) %>% ## add a date column
  mutate(year = str_split_fixed(.$date, 
                                pattern = "\\.", n = 2)[,1]) %>% ## add a year column
  group_by(year) ## group by year


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
    if(!any(is.na(ts_chunk$true_ts))) {
    
      ## get length
      L = nrow(ts_chunk)
      
      ## preprocess the time series:
      ## a. subtracting mean
      ts_true <- ts_chunk$true_ts - mean(ts_chunk$true_ts)
      # ts_noise <- noise - mean(noise)
      
      ## b. windowing - multiply by a parabolic window 
      window_true <- parabolic_window(series = ts_true, N = L)
      ts_true <- ts_true*window_true
      # window_noise <- parabolic_window(series = ts_noise, N = L)
      # ts_noise <- ts_noise*window_noise
      
      ## c. bridge detrending (endmatching)
      ## ie. subtracting from the data the line connecting the first and last points of the series
      ts_true <- bridge_detrender(windowed_series = ts_true, N = L)
      # ts_noise <- bridge_detrender(windowed_series = ts_noise, N = L)
      
      ## calculate spectral exponent in window using PSD and AWC methods
      true_exp_PSD <- spectral_exponent_calculator_PSD(ts_true, l = L)
      # noise_exp_PSD <- spectral_exponent_calculator_PSD(ts_noise, l = L)
      
      true_exp_PSD_low <- true_exp_PSD[[1]]
      true_exp_PSD_high <- true_exp_PSD[[2]]
      true_exp_PSD_all <- true_exp_PSD[[3]]
      # noise_exp_PSD_low <- noise_exp_PSD[[1]]
      # noise_exp_PSD_high <- noise_exp_PSD[[2]]
      # noise_exp_PSD_all <- noise_exp_PSD[[3]]
      
      if(is.na(noise)[[1]]) {
        ## save first time chunk 
        noise <- ts_chunk$true_ts
        
        ts_noise <- noise - mean(noise)
        window_noise <- parabolic_window(series = ts_noise, N = L)
        ts_noise <- ts_noise*window_noise
        ts_noise <- bridge_detrender(windowed_series = ts_noise, N = L)
        
        noise_exp_PSD <- spectral_exponent_calculator_PSD(ts_noise, l = L)
        
        noise_exp_PSD_low <- noise_exp_PSD[[1]]
        noise_exp_PSD_high <- noise_exp_PSD[[2]]
        noise_exp_PSD_all <- noise_exp_PSD[[3]]
      }
      
      all_windows[[element]] <- c(true_exp_PSD_low, true_exp_PSD_high, true_exp_PSD_all,
                                  noise_exp_PSD_low, noise_exp_PSD_high, noise_exp_PSD_all,
                                    year_start, year_stop, paste(n, "years"))
      
    }
    
    ## move to next window
    year_start = year_stop + 1
    year_stop = year_stop + n 
    
    element = element + 1
  }
  
  ## move to next window width
  n = n + 1
  noise = NA
}

all_windows <- as.data.frame(do.call("rbind", all_windows))

colnames(all_windows) <- c("true_exp_PSD_low", "true_exp_PSD_high", "true_exp_PSD_all",
                           "noise_exp_PSD_low", "noise_exp_PSD_high", "noise_exp_PSD_all",
                           "window_start_year","window_stop_year","time_window_width")

all_windows[,1:8] <- sapply(all_windows[,1:8], as.numeric)

true_model_output_PSD_low <- all_windows %>%
  group_by(time_window_width) %>%
  do(tidy(lm(., formula = true_exp_PSD_low ~ window_start_year))) %>%
  filter(term == "window_start_year") %>%
  select(-term)

noise_model_output_PSD_low <- all_windows %>%
  group_by(time_window_width) %>%
  do(tidy(lm(., formula = noise_exp_PSD_low ~ window_start_year))) %>%
  filter(term == "window_start_year") %>%
  select(-term)

## visualize 
gather(all_windows, key = "ts_type", value = "spec_exp", c(true_exp_PSD_low, noise_exp_PSD_low)) %>%
  ggplot(aes(x = window_start_year, y = spec_exp, colour = ts_type)) + geom_point() +
  facet_wrap(~time_window_width)


## next: in Berkeley Earth forcing analysis, force also with fake time series repeating first time window 
## but also with time series that is a bunch of small time series with that beta appended to each other
## or sampled around a mean? ie. take real ts_chunks, ample randomly around mean = first colour, conglomerate them
  
