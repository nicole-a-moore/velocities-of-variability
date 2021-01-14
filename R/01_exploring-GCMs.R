## checking out sample GCM - "CMCC_CESM"
library(tidyverse)
library(ncdf4)
library(abind)
library(sp)
library(sf)
library(raster)
library(broom)
library(evobiR)
library(gridExtra)
library(lubridate)
select <- dplyr::select

## make list of filenames to open:
cmcc_cesm_historical <- c("18500101-18541231.nc",
                          "18550101-18591231.nc",
                          "18600101-18641231.nc",
                          "18650101-18691231.nc",
                          "18700101-18741231.nc",
                          "18750101-18791231.nc",
                          "18800101-18841231.nc",
                          "18850101-18891231.nc",
                          "18900101-18941231.nc",
                          "18950101-18991231.nc",
                          "19000101-19041231.nc",
                          "19050101-19091231.nc",
                          "19100101-19141231.nc",
                          "19150101-19191231.nc",
                          "19200101-19241231.nc",
                          "19250101-19291231.nc",
                          "19300101-19341231.nc",
                          "19350101-19391231.nc",
                          "19400101-19441231.nc",
                          "19450101-19491231.nc",
                          "19500101-19541231.nc",
                          "19550101-19591231.nc",
                          "19600101-19641231.nc",
                          "19650101-19691231.nc",
                          "19700101-19741231.nc",
                          "19750101-19791231.nc",
                          "19800101-19841231.nc",
                          "19850101-19891231.nc",
                          "19900101-19941231.nc",
                          "19950101-19991231.nc")

cmcc_cesm_rcp85 <- c("20000101-20041231.nc",
                     "20050101-20051231.nc",
                     "20060101-20151231.nc",
                     "20160101-20251231.nc",
                     "20260101-20351231.nc",
                     "20360101-20451231.nc",
                     "20460101-20551231.nc",
                     "20560101-20651231.nc",
                     "20660101-20751231.nc",
                     "20760101-20851231.nc",
                     "20860101-20951231.nc",
                     "20960101-21001231.nc")

i = 1
while (i < length(cmcc_cesm_historical)+length(cmcc_cesm_rcp85)+1) {
  if (i < length(cmcc_cesm_historical)+1) {
    filename <- paste("/Volumes/ADATA HV620/GCM-test/tas_day_CMCC-CESM_historical_r1i1p1_",
                      cmcc_cesm_historical[i], sep = "")
  }
  else {
    filename <- paste("/Volumes/ADATA HV620/GCM-test/tas_day_CMCC-CESM_rcp85_r1i1p1_",
                      cmcc_cesm_rcp85[i-length(cmcc_cesm_historical)], sep = "")
  }
  
  ncfile <- nc_open(filename)
  
  ## for first file, extract latitude, longitude, air temp and time 
  if (i == 1) {
    lat <- ncvar_get(ncfile, "lat_bnds") 
    long <- ncvar_get(ncfile, "lon_bnds") 
    tas <- ncvar_get(ncfile, "tas") ## air surface temperature in an array of [long, lat, time]
    time <- ncvar_get(ncfile, "time_bnds")
  }
  ## for other files, append air temp and time onto existing variables 
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

dimnames(tas) <- NULL

## check that time series have 91,676 days
length(tas[1,1,]) ## it does

## make lat and long coords represent centre of grid cells:
lat <- (lat[1,] + lat[2,]) / 2
long <- (long[1,] + long[2,]) / 2 

## reorganize data so North America is left of Europe when plotted:
long[which(long >= 180)] = long[which(long >= 180)] - 360

## try plotting air temps at one time point:
tp <-  expand.grid(long, lat)
colnames(tp) <- c("longitude", "latitude") 
tp$temp <- as.vector(tas[,,1])

## plot in base R:
r = raster(nrows = 48, ncols = 96, xmn=-180, xmx=176.25, ymn=-87.65043, ymx= 87.65043, 
           crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

raster <- rasterize(tp[,1:2], r, tp[,3])
plot(raster, asp = 1)

## plot in base ggplot:
gg1 <- tp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = temp)) + 
  ggtitle("Mean air surface temperature on Jan 1, 1850") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Temperature (K)")




## for now, don't worry about standardizing grid size (will have to standardize across GCMs, though)
## just try FFT 

## create function for making dates since 1850-01-01:
as.date <- function(x, origin = getOption("date.origin")){
  origin <- ifelse(is.null(origin), "1970-01-01", origin)
  as.Date(x, origin)
}
options(date.origin = "1849-12-31")

## 1. Linearly detrend temp time series for each spatial cell separately:
l_detrended_tas <- tas
x = 1 ## represents longitude index
while (x < length(long)+1) {
  y = 1 ## represents latitude index
  while (y < length(lat)+1) {
    local_ts <- l_detrended_tas[x, y, ] ## get the local time series 
    ts_df <- data.frame(time = 1:length(local_ts), ## add simple time column representing days from 1850-01-01
                        temp = local_ts, 
                        date = as.date(1:length(local_ts))) ## add a date column  
    
    ## run linear regression for grid cell 
    output <- lm(ts_df, formula = temp ~ time)
    
    ## extract residuals and add to l_detrended_tas:
    l_detrended_tas[x, y, ] <- output$residuals
    
    y = y + 1
  }
  x = x + 1
}


## save detrended time series object:
##saveRDS(l_detrended_tas, "/Volumes/ADATA HV620/GCM-test/l-detrended-tas.rds")
l_detrended_tas <- readRDS("/Volumes/ADATA HV620/GCM-test/l-detrended-tas.rds")

## try plotting one time point of the new set of detrended temp time series:
detrended_tp <-  expand.grid(long, lat)
colnames(detrended_tp) <- c("longitude", "latitude") 
detrended_tp$temp <- as.vector(l_detrended_tas[,,1])

## plot in base R:
raster <- rasterize(detrended_tp[,1:2], r, detrended_tp[,3])
plot(raster, asp = 1)

## plot in ggplot:
gg2 <- detrended_tp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = temp)) + 
  ggtitle("Linearly detrended air surface temperature (Jan 1, 1850)") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Temperature residual (K)")


## now try with seasonally detrended data
## 2. Seasonally detrend temp time series for each spatial cell separately:
s_detrended_tas <- tas
x = 1 ## represents longitude index
while (x < length(long)+1) {
  y = 1 ## represents latitude index
  while (y < length(lat)+1) {
    local_ts <- s_detrended_tas[x, y, ] ## get the local time series 
    
    ts_df <- data.frame(time = 1:length(local_ts), ## add simple time column representing days from 1850-01-01
                        temp = local_ts, 
                        date = as.date(1:length(local_ts))) ## add a date column
    
    ## compute the within-year temperature profile at each geographical location and subtract from temp
    ts_df <- ts_df %>%
      mutate(year = str_split_fixed(.$date, pattern = "-", n = 2)[,1]) %>% ## add a year column
      group_by(year) %>% ## group by year
      do(mutate(., julian_date = seq(1:length(.$year)))) %>% ## add julian date column
      ungroup() %>%
      group_by(julian_date) %>%
      do(mutate(., temp_profile = mean(.$temp))) %>% ## compute mean temp for each day of year across all years
      ungroup() %>%
      mutate(s_detrended_temp = temp - temp_profile)
    
    ## run linear regression for grid cell 
    output <- lm(ts_df, formula = s_detrended_temp ~ time)
    
    ## extract residuals and add to detrended_tas:
    s_detrended_tas[x, y, ] <- output$residuals
    
    print(paste("On x = ",x, ", y = ", y, sep = ""))
    y = y + 1
  }
  x = x + 1
}


## save detrended time series object:
##saveRDS(s_detrended_tas, "/Volumes/ADATA HV620/GCM-test/s-detrended-tas.rds")
s_detrended_tas <- readRDS("/Volumes/ADATA HV620/GCM-test/s-detrended-tas.rds")


## 3. Compute periodogram using FFT across all years (to start) at each location 
spec2D <- array(dim = c(48, 96)) ## 2D array to store spec exp in 
x = 1
while (x < length(long)+1) {
  y = 1 
  while (y < length(lat)+1) {
    local_ts <- l_detrended_tas[x, y, ] ## get the local detrended time series 
    
    ## calculate periodogram using FFT
    L <- length(local_ts)
    
    # Fourier transform the series: 
    dft <- fft(local_ts)/L
    amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
    amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
    freq <- 1:(L/2)/L ## frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
    
    ## create periodogram data by squaring amplitude of FFT output
    spectral <- data.frame(freq = freq, power = amp^2)
    
    # ## plot spectrum:
    # spectral %>% 
    #   ggplot(aes(x = freq, y = power)) + geom_line() +
    #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")
    
    ## get estimate of spectral exponent over whole time series:
    model_output <- lm(spectral, formula = log10(power) ~ log10(freq)) %>%
      tidy(.) %>%
      filter(term == "log10(freq)") %>%
      mutate(lat = lat[y]) %>%
      mutate(long = long[x]) %>%
      mutate(time_interval = "whole_ts")
    
    ## store average spectral exponent in 2D array represting a static picture of projected variability:
    spec2D[y, x] <- model_output$estimate
    
    ## store:
    if(x == 1 & y == 1) {
      spec_exp <- model_output
    }
    else {
      spec_exp <- rbind(spec_exp, model_output)
    }
    
    y = y + 1
  }
  x = x + 1
}

dimnames(spec2D) <- NULL

## try plotting spectral exponents across space:
sp <-  expand.grid(lat, long)
colnames(sp) <- c("latitude", "longitude") 
sp$spec <- as.vector(spec2D)
sp <- sp %>%
  select(longitude, latitude, spec)

## plot in base R:
raster <- rasterize(sp[,1:2], r, sp[,3])
plot(raster, asp = 1)

## plot in ggplot:
gg3 <- sp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = spec)) + 
  ggtitle("Average spectral exponent over years 1850-2100") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "")


## Notes:
##  - dates each file spans is in title of file
##  - bnds means bounds, rows of lat and lon describe extent of each grid square 
##  - must check for gaps in time series and see how they deal with leap years!!!

## idea to find gaps in the time series:
## which(!(time[,2] - time[,1]) == 1)
## i think this will work?


ggsave(path = "figures/", filename = "air-temps-raw_Jan-01-2000.png", gg1, height = 6, width = 9)
ggsave(path = "figures/", filename = "air-temps-detrended_Jan-01-2000.png", gg2, height = 6, width = 9)
ggsave(path = "figures/", filename = "static-spec-exp.png", gg3, height = 6, width = 9)


## function to calculate spectral exponent over time series window
spectral_exponent_calculator <- function(ts_window) {
  l <- length(ts_window)
  
  # Fourier transform the time series window: 
  dft <- fft(ts_window)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
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

## function to calculate spectral exponent over sliding time series windows of varying widths (5-10 years) with a time step of one year 
## returns a list containing: matrix of change in spectral exponents using each window wdith, data frame of spectral exponents within each sliding window across all locations and for all window widths 
sliding_window_spec_exp <- function(detrended_ts) {
  spec3D <- array(dim = c(48, 96, 6)) ## 3D array to store spec exp in where third dimension represents years spectral exponent was calculated over (5,6,7,8,9,10)
  x = 1
  while (x < length(long)+1) {
    y = 1 
    while (y < length(lat)+1) {
      local_ts <- detrended_tas[x, y, ] ## get the local detrended time series 
      
      #########################################
      ##        SENSITIVITY ANALYSIS:        ##
      #########################################
      ## calculate spectral exponent using FFT over n year windows
      ## store spectral exponents and calculate slope
      n = 5
      while (n < 11) {
        ## find width of window in days (365.25*n)
        window_width <- round(365*n, digits = 0)
        
        ## calcualte spectral exponent in each window of n years, moving one window width at a time
        exp <- SlidingWindow(FUN = spectral_exponent_calculator, data = local_ts, step = 365, 
                             window = window_width)
        
        spec_exp <- exp %>%
          as.data.frame() %>%
          mutate(time_window_start = seq(from = 1, by = 365, 
                                         to = (floor(length(local_ts)-window_width)))) %>%
          rename(spec_exp = ".") %>%
          mutate(lat = lat[y]) %>%
          mutate(long = long[x]) %>%
          mutate(time_window_width = paste(n, "years")) %>%
          mutate(time_step = "1y")
        
        # ## plot:
        # spec_exp %>%
        #   ggplot(., aes (x = time_window_start, y = spec_exp)) + geom_point() +
        #   geom_smooth(method = "lm") + labs(x = "Time window start index",
        #                                     y = "Spectral exponent over window")
        
        ## regress spectral exponent and extract slope representing change in spectral exponent over time
        model_output <- lm(spec_exp, formula = spec_exp ~ time_window_start) %>%
          tidy(.) %>%
          filter(term == "time_window_start")
        
        spec_exp <- cbind(spec_exp, model_output)
        
        ## store in database 
        if(x == 1 & y == 1 & n == 5) {
          all_spec_exp <- spec_exp 
        }
        else {
          all_spec_exp <- rbind(all_spec_exp, spec_exp)
        }
        
        ## store in array
        spec3D[y, x, n-4] <- model_output$estimate
        
        ## advance to next time window width
        print(paste("On x = ",x, ", y = ", y, ", window width = ", n, " years.", sep = ""))
        n = n + 1
      }
      
      
      ## advance to next longitude
      y = y + 1
    }
    ## advance to next latitude
    x = x + 1
  }
  
  dimnames(spec3D) <- NULL
  
  return(list(spec3D, all_spec_exp))
}



## 4. calculate spectral exponent in each sliding window of varying widths for both linearly detrended data and seasonally detrended data:

l_spec_exp <- sliding_window_spec_exp(l_detrended_tas)
##saveRDS(l_spec_exp[[1]], "./data-processed/spec3d_year-time-step_l.rds")
write.csv(l_spec_exp[[2]], "./data-processed/cmcc-cesm_spectral-exponent-sliding-window_l.csv",
          row.names = FALSE)

s_spec_exp <- sliding_window_spec_exp(s_detrended_tas)
##saveRDS(s_spec_exp[[1]], "./data-processed/spec3d_year-time-step_s.rds")
write.csv(s_spec_exp[[2]], "./data-processed/cmcc-cesm_spectral-exponent-sliding-window_s.csv",
          row.names = FALSE)


spec3D <- s_spec_exp[[1]]

## try plotting:
n5 <- spec3D[,,1] ## extract change in spectral exponent over time windows with 5 year width

sp <-  expand.grid(lat, long)
colnames(sp) <- c("latitude", "longitude") 
sp$spec <- as.vector(n5)
sp <- sp %>%
  select(longitude, latitude, spec)

gg4 <- sp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = spec)) + 
  ggtitle("Change in spectral exponent over years 1850-2100 (5 year window, 1 year step)") + 
  coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Slope of spectral exponent")

## woohoo!!!!!
## before moving on to rest of the GCMs, lets try and calculate some velocities of variation! 
