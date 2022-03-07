## script to investigate link between colour of noise in real climate data and population dynamics
library(tidyverse)
library(parallel)
library(MASS)
library(foreach)
library(doParallel)
library(ncdf4)
select <- dplyr::select

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
numCores

registerDoParallel(numCores)  # use multicore, set to the number of our cores

#################################################
###                 time window                ## 
#################################################
## write equation for temperature dependence of little r: 
d_j = 0.03 # juvenile mortality at reference temperature 
d_a = 0.05 # adult mortality rate at reference temperature 
Ad_j = 7500 # Arrhenius constant for juvenile mortality  
Ad_a = 10000 # Arrhenius constant for adult mortality  

Aa = -8000 # Arrhenius constant for development
alpha = 60
b_tr = 50 # average per capita fecundity at reference temperature 

Topt = 298 # optimal temperature
s = 4.8
Tr = 294 # reference temperature 

equation <- function(temp) {
  TD = ((1/Tr) - (1/temp))
  term1 = d_a*exp(Ad_a*TD)
  term2 = (1/(alpha*exp(Aa*TD))) * (pracma::lambertWp(b_tr*alpha*exp(Aa*TD - (temp-Topt)^2/(2*s^2) + alpha*exp(Aa*TD)*(d_a*exp(Ad_a*TD) - d_j*exp(Ad_j*TD)))))
  rT = -term1 + term2

  return(c(term1, term2, rT))
}

## recreate examples from figure 1
range <- seq(from = -5, to = 42, by = 0.1)
temps <- c(273.15 + range)
high_s <- sapply(FUN = equation, temps)
high_s = data.frame(term1 = high_s[1,], term2 = high_s[2,], rT = high_s[3,], temps = temps - 273.15)

## get rid of pts beyond tmin and tmax
high_s = high_s[-(1:first(which(high_s$rT > 0))),]
high_s = high_s[-(last(which(high_s$rT > 0)):nrow(high_s)),]
ggplot(high_s, aes(x = temps, y = term2)) + geom_smooth()

s = 2.5
low_s <- sapply(FUN = equation, temps)
plot(test[3,])
  
## thinking:
## want to vary growth rate from 0-2 (range of r producing stable dynamics in the Ricker model)
## want environmental temperatures to remain within Tmin and Tmax (so that r =/= 0)


## set model parameters
K0 = 100
N0 = 100
l = seq(0.1, 1, by = 0.1)

## read in a spatial chunk of climate data 
path = "/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/"
filepath = paste(path, "sp_files.rds", sep = "")
names = readRDS(filepath)

## read in the ocean coordinates
ocean <- read.csv("data-processed/masks/cmip5-ocean-coords.csv")
ocean_coords <- paste(ocean$lat, ocean$lon)

lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 

lon_index = 0
lat_index = 0
count = 1

## retrieve spatial chunk from nc file
s_filenames <-  str_replace_all(names, "spatial_temps", 's-detrended')
s_open = nc_open(s_filenames[count])
s_detrended_tas = ncvar_get(s_open, "var1_1")
nc_close(s_open)

all_results <- list()
spec_exp_list <- list()
element <- 1
count <- 1
x = 1 ## longitude
while (x < ncol(s_detrended_tas)+1) {
  y = 1 ## latitude
  while (y < nrow(s_detrended_tas)+1) {
    ## if not in the ocean
    if (!paste(lat[y + lat_index], lon[x + lon_index]) %in% ocean_coords) {
      local_ts <- s_detrended_tas[y,x,] ## get the local time series
      local_ts <- local_ts - mean(local_ts) ## centre around 0
      
      ## if no time series, skip to next latitude
      if (length(which(is.na(local_ts))) == 83950) {
        y = y + 1
      }
      else {
        ## run population model:
        print(paste("Running population model for lat = ", lat[y + lat_index],
                    ", lon = ", lon[x + lon_index], sep = ""))
        
        ## to see effects of demographic stochasticity, run each 10x
        ds <- 1
        while (ds <=  10) {
          ## parallelize:
          all <- foreach(icp=1:length(l), .combine=rbind) %dopar% {
            N = N0 ## set starting population size to 100
            K <- K0 + round(local_ts, digits = 0) ## generate new carrying capacity
            K[which(K <= 0)] <- 5 ## make sure none are negative
            
            i=1
            while(i <= length(local_ts)) {
              Nt = N*exp(1.5*(1 - (N/K[i])^l[icp]))
              Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
              
              ## stop if the pop goes extinct:
              if (is.na(Nt) | Nt == 0 | i == length(local_ts)) {
                extinction_time = i
                i = length(local_ts)+1
              }
              ## otherwise continue
              else {
                N = Nt   
                i = i + 1
              }
            }
            c(N = N0, K = K[extinction_time], l = l[icp], t = extinction_time, type = "local_ts",
              lat = lat[y + lat_index], lon = lon[x + lon_index])
          }
          all_results[[count]] <- all
          count = count + 1
          ds = ds + 1
        }
        y = y + 1
      }
    }
    else {
      ## if it was an ocean coordinate, move along
      print(paste("Skipping ocean (lat = ", lat[y + lat_index],
                  ", lon = ", lon[x + lon_index], ")", sep = ""))
      
      y = y + 1
    }
  }
  x = x + 1
}

#saveRDS(all_results, "data-processed/ricker_gcm_s_local_ts_samevar.rds")

## read in colours of noise for chunk:
file = "/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/spec-exp_long-0-60_lat-90-30.csv"
colours <- read.csv(file)

cols <- filter(colours, time_window_width == '10 years')

## calculate average noise colour across sliding windows in each loc
cols <- cols %>%
  group_by(lat, lon) %>%
  mutate(l_spec_exp_PSD_low = mean(l_spec_exp_PSD_low),
         s_spec_exp_PSD_low = mean(s_spec_exp_PSD_low),
         l_spec_exp_AWC = mean(l_spec_exp_AWC),
         s_spec_exp_AWC = mean(s_spec_exp_AWC),
         s_spec_exp_PSD_high = mean(s_spec_exp_PSD_high)) %>%
  ungroup() %>%
  select(lat, lon, l_spec_exp_PSD_low, s_spec_exp_PSD_low, l_spec_exp_AWC, s_spec_exp_AWC,
         s_spec_exp_PSD_high) %>%
  unique(.)

#all_results <- readRDS("data-processed/ricker_gcm_s_local_ts_samevar_reg-awc.rds")
all_results <- as.data.frame(do.call("rbind", all_results))
results <- filter(all_results, t > 1)
results$lat <- as.numeric(results$lat)
results$lon <- as.numeric(results$lon)
results$t = as.numeric(results$t)

results %>%
  ggplot(aes(x = as.numeric(K))) + geom_histogram()

results <- left_join(results, cols)

sims <- readRDS("data-processed/ricker_parallel_all-results_sandbox.rds")
sims <- as.data.frame(do.call("rbind", sims))
sims <- filter(sims, t > 1)

sims %>%
  filter(colour %in% c(0.3, 0.4, 0.5, 0.6, 0.7)) %>%
  ggplot(aes(x = K)) + geom_histogram()

sims <- select(sims, -true_colour) %>%
  mutate(type = "simulated noise")

both_datasets <- results %>%
  mutate("colour" = s_spec_exp_PSD_low) %>%
  select(colour, l, t, K, N) %>%
  mutate(type = "real noise") %>%
  rbind(sims, .)

both_datasets %>%
  ggplot(., aes(x = colour, y = t, colour = type, group = l)) + 
  geom_point() +
  scale_y_log10() +
  facet_wrap(~l)

results %>%
  select(s_spec_exp_PSD_low, l, t, K, N, s_spec_exp_PSD_high)%>%
  ggplot(., aes(x = s_spec_exp_PSD_low, y = t, colour = s_spec_exp_PSD_high)) + 
  geom_point() +
  scale_y_log10() +
  facet_wrap(~l)

results %>%
  filter(l==1) %>%
  ggplot(., aes(x = s_spec_exp_PSD_low, y = s_spec_exp_PSD_high, colour = t)) + 
  geom_point()

## as low spectral exponent increases, high spectral exponent decreases
results %>%
  mutate(diff = s_spec_exp_PSD_high - s_spec_exp_PSD_low) %>%
  select(s_spec_exp_PSD_low, l, t, K, N, s_spec_exp_PSD_high, diff) %>%
  ggplot(., aes(x = s_spec_exp_PSD_low, y = s_spec_exp_PSD_high, colour = diff)) + 
  geom_point() +
  scale_y_log10() +
  facet_wrap(~l)

  
tas_results <- results


## next: 
## add an sst chunk


#################################################
###                sst chunk                   ## 
#################################################
## set model parameters
K0 = 30
N0 = 30
l = seq(0.1, 1, by = 0.1)

## read in a spatial chunk of climate data 
path = "/Volumes/SundayLab/CMIP5-GCMs_tos/01_CMCC-CESM/"
filepath = paste(path, "sp_files_tos.rds", sep = "")
names = readRDS(filepath)

lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 

lon_index = 0
lat_index = 0
count = 1

## retrieve spatial chunk from nc file
s_filenames <-  str_replace_all(names, "spatial_temps", 's-detrended')
s_open = nc_open(s_filenames[count])
s_detrended_tas = ncvar_get(s_open, "var1_1")
nc_close(s_open)

all_results <- list()
spec_exp_list <- list()
element <- 1
count <- 1
x = 1 ## longitude
while (x < ncol(s_detrended_tas)+1) {
  y = 1 ## latitude
  while (y < nrow(s_detrended_tas)+1) {
    ## if not on land ocean
    if (!is.na(s_detrended_tas[y,x,1])) {
      local_ts <- s_detrended_tas[y,x,] - 273.15 ## get the local time series
      local_ts <- local_ts*1/sqrt(var(local_ts))*sqrt(4140)## add variance 
      local_ts <- local_ts - mean(local_ts) ## centre around 0
      
      ## if no time series, skip to next latitude
      if (length(which(is.na(local_ts))) == 83950) {
        y = y + 1
      }
      else {
        ## run population model:
        print(paste("Running population model for lat = ", lat[y + lat_index],
                    ", lon = ", lon[x + lon_index], sep = ""))
        
        ## to see effects of demographic stochasticity, run each 10x
        ds <- 1
        while (ds <=  10) {
          ## parallelize:
          all <- foreach(icp=1:length(l), .combine=rbind) %dopar% {
            N = N0 ## set starting population size to 100
            K <- K0 + round(local_ts, digits = 0) ## generate new carrying capacity
            K[which(K <= 0)] <- 5 ## make sure none are negative
            
            i=1
            while(i <= length(local_ts)) {
              Nt = N*exp(1.5*(1 - (N/K[i])^l[icp]))
              Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
              
              ## stop if the pop goes extinct:
              if (is.na(Nt) | Nt == 0 | i == length(local_ts)) {
                extinction_time = i
                i = length(local_ts)+1
              }
              ## otherwise continue
              else {
                N = Nt   
                i = i + 1
              }
            }
            c(N = N0, K = K[extinction_time], l = l[icp], t = extinction_time, type = "local_ts",
              lat = lat[y + lat_index], lon = lon[x + lon_index])
          }
          all_results[[count]] <- all
          count = count + 1
          ds = ds + 1
        }
        y = y + 1
      }
    }
    else {
      ## if it was a land coordinate, move along
      print(paste("Skipping land (lat = ", lat[y + lat_index],
                  ", lon = ", lon[x + lon_index], ")", sep = ""))
      
      y = y + 1
    }
  }
  x = x + 1
}

#saveRDS(all_results, "data-processed/ricker_gcm_s_local_ts_samevar_ocean.rds")

## read in colours of noise for chunk:
file = "/Volumes/SundayLab/CMIP5-GCMs_tos/01_CMCC-CESM/spec-exp_long-0-60_lat-90-30_new.csv"
colours <- read.csv(file)

cols <- filter(colours, time_window_width == '10 years')

## calculate average noise colour across sliding windows in each loc
cols <- cols %>%
  group_by(lat, lon) %>%
  mutate(l_spec_exp_PSD_low = mean(l_spec_exp_PSD_low),
         s_spec_exp_PSD_low = mean(s_spec_exp_PSD_low),
         l_spec_exp_AWC = mean(l_spec_exp_AWC),
         s_spec_exp_AWC = mean(s_spec_exp_AWC),
         s_spec_exp_PSD_high = mean(s_spec_exp_PSD_high)
         ) %>%
  ungroup() %>%
  select(lat, lon, l_spec_exp_PSD_low, s_spec_exp_PSD_low, l_spec_exp_AWC, s_spec_exp_AWC,
         s_spec_exp_PSD_high) %>%
  unique(.)

#all_results <- readRDS("data-processed/ricker_gcm_s_local_ts_samevar_ocean.rds")
all_results <- as.data.frame(do.call("rbind", all_results))
results <- filter(all_results, t > 1)
results$lat <- as.numeric(results$lat)
results$lon <- as.numeric(results$lon)
results$t = as.numeric(results$t)

results %>%
  ggplot(aes(x = as.numeric(K))) + geom_histogram()

results <- left_join(results, cols)

## join with tas:
tas_results <- select(tas_results, -colour_AWC, -colour_reg)
tas_results$realm = "land"
results$realm = "ocean"
results <- rbind(results, tas_results)

sims <- readRDS("data-processed/ricker_parallel_all-results_sandbox.rds")
sims <- as.data.frame(do.call("rbind", sims))
sims <- filter(sims, t > 1)

sims %>%
  filter(colour %in% c(0.3, 0.4, 0.5, 0.6, 0.7)) %>%
  ggplot(aes(x = K)) + geom_histogram()

sims <- select(sims, -true_colour) %>%
  mutate(type = "simulated noise")

both_datasets <- results %>%
  rename("colour" = s_spec_exp_PSD_low) %>%
  select(colour, l, t, K, N) %>%
  mutate(type = "real noise") %>%
  rbind(sims, .)

both_datasets %>%
  ggplot(., aes(x = colour, y = t, colour = type, group = type)) + 
  geom_point() +
  scale_y_log10() +
  facet_wrap(~l) +
  geom_smooth(method = "lm")

results %>%
  rename("colour" = s_spec_exp_PSD_low) %>%
  select(colour, l, t, K, N, s_spec_exp_PSD_high)%>%
  ggplot(., aes(x = colour, y = t, colour = s_spec_exp_PSD_high)) + 
  geom_point() +
  scale_y_log10() +
  facet_wrap(~l)

results %>%
  filter(l==0.9) %>%
  ggplot(., aes(x = s_spec_exp_PSD_low, y = s_spec_exp_PSD_high, colour = t)) + 
  geom_point() +
  geom_abline(slope = 1, intercept = 2)

## ones that have low scaling exponent across low and high frequencies have longer persistence time










#################################################
###                 garbage                    ## 
#################################################
##functions 
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
  
  ## plot spectrum:
  spectral %>%
    ggplot(aes(x = freq, y = power)) + geom_line() +
    scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")
  
  ## get estimate of spectral exponent over time series window:
  ## fit slope to low freqeuncies
  model_output_low <- spectral %>%
    filter(freq < 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  ## fit slope to high freqeuncies
  model_output_high <- spectral %>%
    filter(freq >= 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  ## fit slope for species with generation times < 1year
  model_output_reg <- spectral %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  return(list(-model_output_low$estimate, -model_output_high$estimate, -model_output_reg$estimate))
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

## set model parameters
K0 = 30
N0 = 30
l = seq(0.1, 1, by = 0.1)

## read in a spatial chunk of climate data 
path = "/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/"
filepath = paste(path, "sp_files.rds", sep = "")
names = readRDS(filepath)

## read in the ocean coordinates
ocean <- read.csv("data-processed/masks/cmip5-ocean-coords.csv")
ocean_coords <- paste(ocean$lat, ocean$lon)

lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 

lon_index = 0
lat_index = 0
count = 1

## retrieve spatial chunk from nc file
s_filenames <-  str_replace_all(names, "spatial_temps", 's-detrended')
s_open = nc_open(s_filenames[count])
s_detrended_tas = ncvar_get(s_open, "var1_1")
nc_close(s_open)

all_results <- list()
spec_exp_list <- list()
element <- 1
count <- 1
x = 1 ## longitude
while (x < ncol(s_detrended_tas)+1) {
  y = 1 ## latitude
  while (y < nrow(s_detrended_tas)+1) {
    ## if not in the ocean
    if (!paste(lat[y + lat_index], lon[x + lon_index]) %in% ocean_coords) {
      local_ts <- s_detrended_tas[y,x,] - 273.15 ## get the local time series
      local_ts <- local_ts*1/sqrt(var(local_ts))*sqrt(4140)## add variance 
      local_ts <- local_ts - mean(local_ts) ## centre around 0
      
      ## calculate spectral exponent
      ## get length
      L = length(local_ts)
      
      ## preprocess the time series:
      ## a. subtracting mean
      ts_s <- local_ts - mean(local_ts)
      
      ## b. windowing - multiply by a parabolic window 
      window_s <- parabolic_window(series = ts_s, N = L)
      ts_s <- ts_s*window_s
      
      ## c. bridge detrending (endmatching)
      ## ie. subtracting from the data the line connecting the first and last points of the series
      ts_s <- bridge_detrender(windowed_series = ts_s, N = L)
      
      ## calculate spectral exponent in window using PSD and AWC methods
      s_exp_PSD <- spectral_exponent_calculator_PSD(ts_s, l = L)
      
      s_exp_PSD_low <- s_exp_PSD[[1]]
      s_exp_PSD_reg <- s_exp_PSD[[3]]
      
      s_exp_AWC <- spectral_exponent_calculator_AWC(ts_s, N = L)
      
      
      ## if no time series, skip to next latitude
      if (length(which(is.na(local_ts))) == 83950) {
        y = y + 1
      }
      else {
        ## run population model:
        print(paste("Running population model for lat = ", lat[y + lat_index],
                    ", lon = ", lon[x + lon_index], sep = ""))
        
        ## to see effects of demographic stochasticity, run each 10x
        ds <- 1
        while (ds <=  10) {
          ## parallelize:
          all <- foreach(icp=1:length(l), .combine=rbind) %dopar% {
            N = N0 ## set starting population size to 100
            K <- K0 + round(local_ts, digits = 0) ## generate new carrying capacity
            K[which(K <= 0)] <- 5 ## make sure none are negative
            
            i=1
            while(i <= length(local_ts)) {
              Nt = N*exp(1.5*(1 - (N/K[i])^l[icp]))
              Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
              
              ## stop if the pop goes extinct:
              if (is.na(Nt) | Nt == 0 | i == length(local_ts)) {
                extinction_time = i
                i = length(local_ts)+1
              }
              ## otherwise continue
              else {
                N = Nt   
                i = i + 1
              }
            }
            c(N = N0, K = K[extinction_time], l = l[icp], t = extinction_time, type = "local_ts",
              lat = lat[y + lat_index], lon = lon[x + lon_index], colour_AWC = s_exp_AWC,
              colour_reg = s_exp_PSD_reg)
          }
          all_results[[count]] <- all
          count = count + 1
          ds = ds + 1
        }
        y = y + 1
      }
    }
    else {
      ## if it was an ocean coordinate, move along
      print(paste("Skipping ocean (lat = ", lat[y + lat_index],
                  ", lon = ", lon[x + lon_index], ")", sep = ""))
      
      y = y + 1
    }
  }
  x = x + 1
}

#saveRDS(all_results, "data-processed/ricker_gcm_s_local_ts_samevar.rds")
#saveRDS(all_results, "data-processed/ricker_gcm_s_local_ts_samevar_reg-awc.rds")
results <- as.data.frame(do.call("rbind", readRDS("data-processed/ricker_gcm_s_local_ts_samevar.rds"))) 
all_results = readRDS("data-processed/ricker_gcm_s_local_ts_samevar_reg-awc.rds")

## read in colours of noise for chunk:
file = "/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/spec-exp_long-0-60_lat-90-30.csv"
colours <- read.csv(file)

cols <- filter(colours, time_window_width == '10 years')

## calculate average noise colour across sliding windows in each loc
cols <- cols %>%
  group_by(lat, lon) %>%
  mutate(l_spec_exp_PSD_low = mean(l_spec_exp_PSD_low),
         s_spec_exp_PSD_low = mean(s_spec_exp_PSD_low),
         l_spec_exp_AWC = mean(l_spec_exp_AWC),
         s_spec_exp_AWC = mean(s_spec_exp_AWC),
         s_spec_exp_PSD_high = mean(s_spec_exp_PSD_high),) %>%
  ungroup() %>%
  select(lat, lon, l_spec_exp_PSD_low, s_spec_exp_PSD_low, l_spec_exp_AWC, s_spec_exp_AWC, s_spec_exp_PSD_high) %>%
  unique(.)

all_results <- as.data.frame(do.call("rbind", all_results))
all_results$colour_PSDlow = results$colour
results <- filter(all_results, t > 1)
results$lat <- as.numeric(results$lat)
results$lon <- as.numeric(results$lon)
results$t = as.numeric(results$t)
results$colour_PSDlow = as.numeric(results$colour_PSDlow)
results$colour_AWC = as.numeric(results$colour_AWC)
results$colour_reg = as.numeric(results$colour_reg)

results %>%
  ggplot(aes(x = as.numeric(K))) + geom_histogram()

results <- left_join(results, cols)

sims <- readRDS("data-processed/ricker_parallel_all-results_sandbox.rds")
sims <- as.data.frame(do.call("rbind", sims))
sims <- filter(sims, t > 1)

sims %>%
  filter(colour %in% c(0.3, 0.4, 0.5, 0.6, 0.7)) %>%
  ggplot(aes(x = K)) + geom_histogram()

sims <- select(sims, -true_colour) %>%
  filter(colour %in% c("0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1",
                       "1.1", "1.2", "1.3")) %>%
  mutate(type = "simulated noise") 

windowed <- results %>%
  mutate("colour" = s_spec_exp_PSD_low) %>%
  select(colour, l, t, K, N) %>%
  mutate(type = "windowed")

both_datasets <- results %>%
  mutate("colour" = colour_AWC) %>%
  select(colour, l, t, K, N) %>%
  mutate(type = "real noise") %>%
  rbind(sims,.) %>%
  rbind(windowed,. )
