## checking out sample GCM - "CMCC_CESM"
library(tidyverse)
library(ncdf4)
library(abind)
library(sp)
library(sf)
library(raster)
library(broom)
select <- dplyr::select

## make list of filenames to open: 
cmcc_cesm <- c("20960101-21001231.nc",
"20860101-20951231.nc",
"20760101-20851231.nc",
"20660101-20751231.nc",
"20560101-20651231.nc",
"20460101-20551231.nc",
"20360101-20451231.nc",
"20260101-20351231.nc",
"20160101-20251231.nc",
"20060101-20151231.nc",
"20050101-20051231.nc",
"20000101-20041231.nc")

i = 1
while (i < length(cmcc_cesm)+1) {
  filename <- paste("/Volumes/ADATA HV620/GCM-test/tas_day_CMCC-CESM_rcp85_r1i1p1_",
                    cmcc_cesm[i], sep = "")
  ncfile <- nc_open(filename)
  
  ## for first file, extract lattitude, longitude, air temp and time 
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

## make lat and long coords represent centre of grid cells:
lat <- (lat[1,] + lat[2,]) / 2
long <- (long[1,] + long[2,]) / 2

## try plotting air temps at one time point:
tp <-  expand.grid(long, lat)
colnames(tp) <- c("longitude", "latitude") 
tp$temp <- as.vector(tas[,,1])

## plot in base R:
r = raster(nrows = 48, ncols = 96, xmn=0, xmx=356.25, ymn=-87.65043, ymx= 87.65043, 
           crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

raster <- rasterize(tp[,1:2], r, tp[,3])
plot(raster, asp = 1)

## plot in base ggplot:
gg1 <- tp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = temp)) + 
  ggtitle("Mean air surface temperature on Jan 1, 2000") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Temperature (K)")




## for now, don't worry about standardizing grid size (will have to standardize across GCMs)
## just try FFT 

## 1. Linearly detrend temp time series for each spatial cell separately:
detrended_tas <- tas
x = 1 ## represents longitude index
while (x < length(long)+1) {
  y = 1 ## represents latitude index
  while (y < length(lat)+1) {
    local_ts <- detrended_tas[x, y, ] ## get the local time series 
    ts_df <- data.frame(time = 1:length(local_ts), temp = local_ts) ## add simple time column for now (days from 2000-01-01) 
    
    ## run linear regression for grid cell 
    output <- lm(ts_df, formula = temp ~ time)
    
    ## extract residuals and add to detrended_tas:
    detrended_tas[x, y, ] <- output$residuals
    
    y = y + 1
  }
  x = x + 1
}


## try plotting one time point of the new set of detrended temp time series:
detrended_tp <-  expand.grid(long, lat)
colnames(detrended_tp) <- c("longitude", "latitude") 
detrended_tp$temp <- as.vector(detrended_tas[,,1])

## plot in base R:
raster <- rasterize(detrended_tp[,1:2], r, detrended_tp[,3])
plot(raster, asp = 1)

## plot in ggplot:
gg2 <- detrended_tp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = temp)) + 
  ggtitle("Linearly detrended air surface temperature (Jan 1, 2000)") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Temperature residual (K)")


## 2. Compute periodogram using FFT across all years (to start) at each location 
spec <- array(dim = c(48, 96)) ## 2D array to store spec exp in 
x = 1
while (x < length(long)+1) {
  y = 1 
  while (y < length(lat)+1) {
    local_ts <- detrended_tas[x, y, ] ## get the local detrended time series 
    
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
    spec[y, x] <- model_output$estimate
    
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

dimnames(spec) <- list(lat, long)

## try plotting spectral exponents across space:
sp <-  expand.grid(lat, long)
colnames(sp) <- c("latitude", "longitude") 
sp$spec <- as.vector(spec)
sp <- sp %>%
  select(longitude, latitude, spec)

## plot in base R:
raster <- rasterize(sp[,1:2], r, sp[,3])
plot(raster, asp = 1)

## plot in ggplot:
gg3 <- sp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = spec)) + 
  ggtitle("Average spectral exponent over years 2000-2100") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "")


## Notes:
##  - dates each file spans is in title of file
##  - bnds means bounds, rows of lat and lon describe extent of each grid square 
## - muct check for gaps in time series and see how they deal with leap years!!!

## idea to find gaps in the time series:
## which(!(time[,2] - time[,1]) == 1)
## i think this will work?


ggsave(path = "figures/", filename = "air-temps-raw_Jan-01-2000.png", gg1, height = 6, width = 9)
ggsave(path = "figures/", filename = "air-temps-detrended_Jan-01-2000.png", gg2, height = 6, width = 9)
ggsave(path = "figures/", filename = "static-spec-exp.png", gg3, height = 6, width = 9)
