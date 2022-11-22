### script that takes temperature time series from Berkeley Earth and:
## 1) computes the variance of each time series
## 2) preps each time series for use in population model 
##    a) adjusts the variance of each time series to equal the avg. variance 
##    b) converts temperatures to degrees Kelvin
##    c) centres the mean of each time series around 273 K 
##    d) records the min and max temperature in each time series 
library(tidyverse)
library(raster)
library(broom)
library(ncdf4)
library(easyNCDF)
library(foreach)
library(doParallel)

## detect cores
starts <- rep(100, 30)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
cl <- makePSOCKcluster(numCores)

registerDoParallel(numCores)

## get paths and filenames 
filenames <- readRDS("data-processed/BerkeleyEarth/be_sp_files.rds")

#################################################
###          1) compute variance               ## 
#################################################
lon_index = 0
lat_index = 0
count = 1
vars <- c()
maxs <- c()
mins <- c()
ranges <- c()
means <- c()
while(count < length(filenames)+1) {
  
  ## retrieve spatial chunk from nc file
  open = nc_open(filenames[count])
  tas = ncvar_get(open, "var1_1")
  nc_close(open)
  
  x = 1
  while(x <= 60) {
    y = 1
    while (y <= 60) {
      local_ts <- tas[y,x,]
      if(!is.na(max(local_ts))) {
      maxs <- append(maxs, max(local_ts))
      mins <- append(mins, min(local_ts))
      means <- append(means, mean(local_ts))
      vars <- append(vars, var(local_ts))
      ranges <- append(ranges, max(local_ts) - min(local_ts))
      
      }
      y = y+1
    }
    x = x+1
  }
  
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

list <- list(vars, means, maxs, mins, ranges)
saveRDS(list, "data-processed/be_stats.rds")
hist(vars)
avg_var = mean(vars) ## 109.086


#################################################
###      2) prep time series for forcing       ## 
#################################################
lon_index = 0
lat_index = 0
count = 1
avg_var = 109.086
ts_data <- c()
lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 
while(count < length(filenames)+1) {
  
  ## retrieve spatial chunk from nc file
  open = nc_open(filenames[count])
  tas = ncvar_get(open, "var1_1")
  nc_close(open)
  
  ## get lat and lon bounds to extract in between:
  lon_bound1 <- lon_index 
  lon_bound2 <- lon_index + 60
  lat_bound1 <- lat_index
  lat_bound2 <- lat_index - 60

  ## for each cell of the raster
  altered_tas <- array(dim = c(60,60,51830))
  x = 1 ## longitude
  while (x <= ncol(tas)) {
    y = 1 ## latitude
   
    while (y <= nrow(tas)) {
      ## get the local time series
      local_ts <- tas[y,x,]
      
      ## a) adjust the variance to equal the avg. variance
      local_ts = local_ts*1/sqrt(var(local_ts))*sqrt(avg_var)
      
      ## b) convert temperatures to degrees Kelvin
      local_ts = local_ts + 273.15
      
      ## c) centre the mean of the time series around 273.15 K
      diff = 273.15 - mean(local_ts)
      local_ts = local_ts + diff
      
      ## d) records the min and max temperature in each time series 
      min = min(local_ts, na.rm=T)
      max = max(local_ts, na.rm=T)
      
      ## save the new time series 
      altered_tas[y,x,] <- local_ts
      
      ts_data <- rbind(ts_data, 
                       c(min_temp = min, max_temp = max, lat =  lat[y + lat_index],
                         lon = lon[x + lon_index]))
      y = y + 1
    }
    x = x + 1
  }
  
  ## save spatial data 
  ArrayToNc(altered_tas, file_path = paste("data-processed/BerkeleyEarth/be_prepped-temps_lon-", 
                                           lon_bound1,"-", lon_bound2,
                                        "_lat-", lat_bound1 + 90, "-", lat_bound2 + 90,".nc", sep = ""))

  ## advance lat and long indecies to move to next chunk
  if (count == 6) {
    lat_index <- lat_index - 60
    lon_index <- 0
  }
  else if (count == 12) {
    lat_index <- lat_index - 60
    lon_index <- 0
  }
  else {
    lon_index <- lon_index + 60
  }
  
  count = count + 1
}

## save data
write.csv(ts_data, "data-processed/BerkeleyEarth/be_minmax-data.csv", row.names = FALSE)

ts_data <- data.frame(ts_data)
max(ts_data$max_temp)
min(ts_data$min_temp)

hist(ts_data$max_temp)
hist(ts_data$min_temp)
mins <- ts_data$min_temp[!is.infinite(ts_data$min_temp)]
quantile(mins, 0.1)
maxs <- ts_data$max_temp[!is.infinite(ts_data$max_temp)]
quantile(maxs, 0.9)


#################################################
###             3) run simulations             ## 
#################################################
## parameterise the r model using min max data 
## make centred from 227 - 311
fitness <- function(temp) {
  TD = ((1/Tr) - (1/temp)) ## 
  bT = b_tr*exp(-((temp-Topt)^2/(2*s^2))) 
  inv_alphaT = alpha*exp(Aa*TD)
  term1 = d_a*exp(Ad_a*TD)
  term2 = (1/(alpha*exp(Aa*TD)))*(LambertW::W(b_tr*alpha*exp(Aa*TD - ((temp-Topt)^2/(2*s^2)) + alpha*exp(Aa*TD)*(d_a*exp(Ad_a*TD) - d_j*exp(Ad_j*TD)))))
  rT = -term1 + term2
  
  return(c(term1, term2, rT, bT, inv_alphaT))
}

## vary parameters:
s = 15
d_j = 0.03 # juvenile mortality at reference temperature 
d_a = 0.07 # adult mortality rate at reference temperature 
Ad_j = 7500 # Arrhenius constant for juvenile mortality  
Ad_a = 16000 # Arrhenius constant for adult mortality  
Aa = -3000 # Arrhenius constant for development
alpha = 1 # age at maturity
b_tr = 50 # average per capita fecundity at reference temperature
Topt = 273 # temperature at which avg fecundity is maximal
Tr = 298 # reference temperature

# ## call function over range of temperatures
# range <- seq(from = -150, to = 50, by = 0.001)
# temps <- c(273.15 + range)
# data <- sapply(FUN = fitness, temps)
# 
# ## plot results
# data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps,
#                   bT = data[4,], inv_alphaT = data[5,]) %>%
#   select(-bT, inv_alphaT) %>%
#   gather(key = "term", value = "value", c(term1, term2, rT)) %>%
#   mutate(facet = ifelse(term == "rT", "rT", "terms"))
# 
# data %>%
#   ggplot(., aes(x = temps, y = value, colour = term)) + geom_point(size = 0.1) + theme_bw() +
#   labs(x = "Temperature (C)", y = "rm(T) or components of rm(T)") +
#   scale_y_continuous(limits = c(-1, 2)) +
#   facet_wrap(~facet, nrow = 2)  + 
#   geom_vline(xintercept = 223.15) + geom_vline(xintercept = 313.15)

## carrying capacity
M = 20 ## units: kg
Ea = 0.15  ## units: eV
k = 8.617e-5 ## units: eV K-1
Kt = M^(-3/4)*exp(Ea/(k*T))

# range <- seq(from = -50, to = 30, by = 0.001)
# temps <- c(273.15 + range)
# 
# Kt = M^(-3/4)*exp(Ea/(k*temps))
# 
# data = data.frame(Kt = Kt, temps = temps)
# 
# ggplot(data, aes(x = temps, y = Kt)) + geom_point() + theme_bw() +
#   labs(x = "Temperature (K)", y = "Carrying capacity") + 
#   geom_vline(xintercept = 223.15) + geom_vline(xintercept = 313.15)


## run population model simulations
## - force r, K and both r and K with temperature time series 
## - store population dynamics
## - 64800 x 100 x 4 

##### Ricker model #####
N0 = 100
Tmax = 200000 
K0 = 100
l = seq(0.1, 1, by = 0.1)
# l = seq(0.1, 1, by = 0.3)

lon_index = 0
lat_index = 0
count = 1
icp = 1
lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 
while(count < length(filenames)+1) {
  
  ## get lat and lon bounds to extract in between:
  lon_bound1 <- lon_index 
  lon_bound2 <- lon_index + 60
  lat_bound1 <- lat_index
  lat_bound2 <- lat_index - 60

  x = 1 ## longitude
  while (x <= 60) {
    y = 1 ## latitude
    
    ## retrieve local time series from chunk of prepped temps from nc file
    open = nc_open(paste("data-processed/BerkeleyEarth/be_prepped-temps_lon-",  lon_bound1,"-", 
                         lon_bound2, "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    tas = ncvar_get(open, "var1_1")
    tas <- tas[1:60,x,]
    nc_close(open)
    
    while (y <= nrow(tas)) {
      
      if(file.exists(paste("data-processed/BerkeleyEarth/popdynamics_lat-", lat[y + lat_index], "_lon-", lon[x + lon_index],
                           "_icp-", l[icp], ".csv", sep = ""))) {
        y = y + 1
        
      }
      else {
        local_ts <- tas[y,]
        
        ## if not na, force the population models 
        if(length(which(is.na(local_ts))) != length(local_ts)) {
          
          ## calculate new carrying capacity according to how it varies with temperature
          K <- M^(-3/4)*exp(Ea/(k*(local_ts))) 
          K[which(K <= 0)] <- 5 ## make sure none are negative
          K <- round(K, 0)
          
          ## C&Y version:
          K_cy <- 100 + (local_ts - 273.15)
          
          ## calculate new growth rate according to how it varies with temperature
          r <- sapply(FUN = fitness, local_ts)[3,]
          r[which(r<(-1))] <- -1
          ## make sure none are below -2
          #plot(ts(r))
          
          ## loop through ICPs
          #icp = 1 
          #while (icp < length(l)+1) {
          print(paste("On lat: ", y, " and lon: ", x, " and icp ", icp,  sep = ""))
          
          ## run 100 simulations per model:
          all <- foreach (z = 1:100, .combine=rbind)  %dopar% {
            
            N = N0 ## set starting population size to 100
            N_r = Nt = N_K = N_K_cy = N  # N_rK 
            Nts <- Nts_r <- Nts_K <- Nts_K_cy <- N #  Nts_rK
            i=1
            while(i < length(local_ts)) {
              
              ## vary r, K, neither and both according to environment
              Nt = N*exp(1.5*(1 - (N/K0)^l[icp])) ## constant r, constant K
              Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
              
              ## K:
              Nt_K = N_K*exp(1.5*(1 - (N_K/K[i])^l[icp])) # constant r, variable K
              Nt_K = sample(rpois(as.numeric(Nt_K), n = 1000), size = 1)
              
              ## C&Y K:
              Nt_K_cy = N_K_cy*exp(1.5*(1 - (N_K_cy/K_cy[i])^l[icp])) # constant r, variable K
              Nt_K_cy = sample(rpois(as.numeric(Nt_K_cy), n = 1000), size = 1)
              
              ## r:
              Nt_r = N_r*exp(r[i]*(1 - (N_r/K0)^l[icp])) # variable r, constant K
              Nt_r = sample(rpois(as.numeric(Nt_r), n = 1000), size = 1)
              
              # ## both:
              # Nt_rK = N_rK*exp(r[i]*(1 - (N_rK/K[i])^l[icp])) # variable r, variable K
              # Nt_rK = sample(rpois(as.numeric(N_rK), n = 1000), size = 1)
              
              N = Nt
              N_r = Nt_r
              N_K = Nt_K
              N_K_cy = Nt_K_cy
              # N_rK = Nt_rK
              Nts <- append(Nts, Nt)
              Nts_K <- append(Nts_K, Nt_K)
              Nts_K_cy <- append(Nts_K_cy, Nt_K_cy)
              Nts_r <- append(Nts_r, Nt_r)
              # Nts_rK <- append(Nts_rK, Nt_rK)
              
              i = i + 1
            }
            
            
            ## return population time series
            data.frame(lat = lat[y + lat_index], lon = lon[x + lon_index], 
                       sim = z, r = r, K = K,  K_cy = K_cy, Nts = Nts, Nts_r = Nts_r, 
                       Nts_K = Nts_K, Nts_K_cy = Nts_K_cy,
                       # Nts_rK = Nts_rK, 
                       l = l[icp], 
                       noise = local_ts, max_temp = max(local_ts), min_temp = min(local_ts))
            
          }
          
          ## save the data for location 
          write.csv(all, paste("data-processed/BerkeleyEarth/popdynamics_lat-", lat[y + lat_index], "_lon-", lon[x + lon_index],
                               "_icp-", l[icp], ".csv", sep = ""), 
                    row.names = F)
          
          #  icp = icp + 1
          #}
        }
        y = y + 1
      }
    }
    x = x + 1
  }
  
  
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

## plot one:
all <- read.csv("data-processed/BerkeleyEarth/popdynamics_lat-63.5_lon-32.5_icp-0.1.csv")

all %>%
    filter(sim == 20) %>%
    mutate(time = rep(1:51830, length(unique(.$sim)))) %>%
    ggplot(., aes(x = time, y = Nts_K_cy, group = sim)) + geom_line() + 
    geom_line(aes(y = r*100), colour = "red")

## how much variation is there between the 100 runs?
all %>%
  filter(sim %in% 1:100) %>%
  mutate(time = rep(1:51830, 100)) %>%
  group_by(time) %>%
  mutate(var = sd(Nts_r), 
         mean = mean(Nts_r)) %>%
  ggplot(., aes(x = time, y = var)) + geom_line() + 
  geom_line(aes(y = mean))


####################################################################
###         4) analyze changes in population dynamics             ## 
####################################################################
## analyze change in noise colour over time windows 
icp = 1
lon_index = 0
lat_index = 0
pop_exp_list <- list()
element = 1
l = seq(0.1, 1, by = 0.1)

lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 

## get lat and lon bounds to extract in between:
lon_bound1 <- lon_index 
lon_bound2 <- lon_index + 60
lat_bound1 <- lat_index
lat_bound2 <- lat_index - 60

## loop through grid cells
x = 1
while (x <= 360) {
  
  y = 1
  while (y <= 360) {
    
    filename = paste("data-processed/BerkeleyEarth/popdynamics_lat-", 
                     lat[y + lon_index], "_lon-", lon[x + lon_index],
                     "_icp-", l[icp], ".csv", sep = "")
    ## read in data:
    if (file.exists(filename)) {
      
      sims <- read.csv(filename)
      
      sims <- group_split(sims, sim)
      
      num <- 1
      while (num <= 100) {
        sim <- sims[[num]]
        
        ## for each time window, process time series and calculate colour
        n = 5
        while (n < 11) {
          year_start <- 0
          year_stop <- n*365 
          
          while (year_start <= (nrow(sim) - n*365)) {
            ## extract temps within time window
            ts_chunk_r <- sim$Nts_r[year_start:year_stop]
            ts_chunk_K <- sim$Nts_K[year_start:year_stop]
            ts_chunk_K_cy <- sim$Nts_K_cy[year_start:year_stop]
            noise <- sim$noise[year_start:year_stop]
            demsto <- sim$Nts[year_start:year_stop]
            
            ## if the window does not have NAs
            if (length(which(ts_chunk_r == 0)) == 0) {
              ## get length
              L = length(ts_chunk_r)
              
              ## preprocess the time series:
              ## a. subtracting mean
              ts_r <- ts_chunk_r - mean(ts_chunk_r)
              
              ## b. windowing - multiply by a parabolic window 
              window_r <- parabolic_window(series = ts_r, N = L)
              ts_r <- ts_r*window_r
              
              ## c. bridge detrending (endmatching)
              ## ie. subtracting from the data the line connecting the first and last points of the series
              ts_r <- bridge_detrender(windowed_series = ts_r, N = L)
              
              ## calculate spectral exponent in window using PSD and AWC methods
              r_exp_PSD <- spectral_exponent_calculator_PSD(ts_r, l = L)
              
              r_exp_PSD_low <- r_exp_PSD[[1]]
              r_exp_PSD_high <- r_exp_PSD[[2]]
              r_exp_PSD_all <- r_exp_PSD[[3]]
              # r_plot <- r_exp_PSD[[4]]
              # 
              # ggsave(r_plot,
              #        path = "figures/pop-spectral-change",
              #        filename = paste(" 0.9_r_", n, "_window-start-", year_start,
              #                         ".png", sep = ""),
              #        device = "png", width = 6, height = 4)
            }
            else {
              r_exp_PSD_low = r_exp_PSD_high = r_exp_PSD_all = NA
            }
            if (length(which(ts_chunk_K == 0)) == 0) {
              L = length(ts_chunk_K)
              ts_K <- ts_chunk_K - mean(ts_chunk_K)
              window_K <- parabolic_window(series = ts_K, N = L)
              ts_K <- ts_K*window_K
              ts_K <- bridge_detrender(windowed_series = ts_K, N = L)
              K_exp_PSD <- spectral_exponent_calculator_PSD(ts_K, l = L)
              
              K_exp_PSD_low <- K_exp_PSD[[1]]
              K_exp_PSD_high <- K_exp_PSD[[2]]
              K_exp_PSD_all <- K_exp_PSD[[3]]
              # K_plot <- K_exp_PSD[[4]]
              
              # ggsave(K_plot, 
              #        path = "figures/pop-spectral-change", 
              #        filename = paste(" 0.9_K_", n, "_window-start-", year_start, 
              #                         ".png", sep = ""), 
              #        device = "png", width = 6, height = 4)
            }
            else {
              K_exp_PSD_low = K_exp_PSD_high = K_exp_PSD_all = NA
            }
            if (length(which(ts_chunk_K_cy == 0)) == 0) {
              L = length(ts_chunk_K_cy)
              ts_K_cy <- ts_chunk_K_cy - mean(ts_chunk_K_cy)
              window_K_cy <- parabolic_window(series = ts_K_cy, N = L)
              ts_K_cy <- ts_K_cy*window_K_cy
              ts_K_cy <- bridge_detrender(windowed_series = ts_K_cy, N = L)
              K_cy_exp_PSD <- spectral_exponent_calculator_PSD(ts_K_cy, l = L)
              
              K_cy_exp_PSD_low <- K_cy_exp_PSD[[1]]
              K_cy_exp_PSD_high <- K_cy_exp_PSD[[2]]
              K_cy_exp_PSD_all <- K_cy_exp_PSD[[3]]
              # K_cy_plot <- K_cy_exp_PSD[[4]]
              
              # ggsave(K_cy_plot, 
              #        path = "figures/pop-spectral-change", 
              #        filename = paste(" 0.9_K_", n, "_window-start-", year_start, 
              #                         ".png", sep = ""), 
              #        device = "png", width = 6, height = 4)
            }
            else {
              K_cy_exp_PSD_low = K_cy_exp_PSD_high = K_cy_exp_PSD_all = NA
            }
            # if (length(which(ts_chunk_rK == 0)) == 0) {
            #   ## get length
            #   L = length(ts_chunk_K_cy)
            #   window_rK <- parabolic_window(series = ts_rK, N = L)
            #   ts_rK <- ts_chunk_rK - mean(ts_chunk_rK)
            #   ts_rK <- ts_rK*window_rK
            #   ts_rK <- bridge_detrender(windowed_series = ts_rK, N = L)
            #   rK_exp_PSD <- spectral_exponent_calculator_PSD(ts_rK, l = L)
            #   
            #   rK_exp_PSD_low <- rK_exp_PSD[[1]]
            #   rK_exp_PSD_high <- rK_exp_PSD[[2]]
            #   rK_exp_PSD_all <- rK_exp_PSD[[3]]
            # }
            # else {
            #   rK_exp_PSD_low = rK_exp_PSD_high = rK_exp_PSD_all = NA
            # }
            if (length(which(noise == 0)) == 0) {
              ## get length
              L = length(noise)
              ts_noise <- noise - mean(noise)
              window_noise <- parabolic_window(series = noise, N = L)
              ts_noise <- ts_noise*window_noise
              ts_noise <- bridge_detrender(windowed_series = ts_noise, N = L)
              noise_exp_PSD <- spectral_exponent_calculator_PSD(ts_noise, l = L)
              
              noise_exp_PSD_low <- noise_exp_PSD[[1]]
              noise_exp_PSD_high <- noise_exp_PSD[[2]]
              noise_exp_PSD_all <- noise_exp_PSD[[3]]
            }
            else {
              noise_exp_PSD_low = noise_exp_PSD_high = noise_exp_PSD_all = NA
            }
            if (length(which(demsto == 0)) == 0) {
              ## get length
              L = length(demsto)
              
              ## preprocess the time series:
              ## a. subtracting mean
              demsto_ts <- demsto - mean(demsto)
              
              ## b. windowing - multiply by a parabolic window 
              window_demsto <- parabolic_window(series = demsto_ts, N = L)
              demsto_ts <- demsto_ts*window_demsto
              
              ## c. bridge detrending (endmatching)
              ## ie. subtracting from the data the line connecting the first and last points of the series
              demsto_ts <- bridge_detrender(windowed_series = demsto_ts, N = L)
              
              ## calculate spectral exponent in window using PSD and AWC methods
              demsto_exp_PSD <- spectral_exponent_calculator_PSD(demsto_ts, l = L)
              
              demsto_exp_PSD_low <- demsto_exp_PSD[[1]]
              demsto_exp_PSD_high <- demsto_exp_PSD[[2]]
              demsto_exp_PSD_all <- demsto_exp_PSD[[3]]
            }
            else {
              demsto_exp_PSD_low = demsto_exp_PSD_high = demsto_exp_PSD_all = NA
            }
            
            ## store:
            pop_exp_list[[element]] <- c(r_exp_PSD_low, r_exp_PSD_high, r_exp_PSD_all,
                                         K_exp_PSD_low, K_exp_PSD_high, K_exp_PSD_all,
                                         K_cy_exp_PSD_low, K_cy_exp_PSD_high, K_cy_exp_PSD_all,
                                         noise_exp_PSD_low, noise_exp_PSD_high, noise_exp_PSD_all,
                                         demsto_exp_PSD_low, demsto_exp_PSD_high, demsto_exp_PSD_all,
                                         year_start, year_stop, paste(n*365, "days"),
                                         l[icp], lat[y + lon_index], lon[x + lon_index], num)
            
            
            
            
            ## move to next window
            year_start = year_stop + 1
            year_stop = year_stop + n*365 
            
            element = element + 1
          }
          
          print(paste("On lat: ", lat[y + lon_index], " and lon: ", lon[x + lon_index], 
                      " and window width ", n,  " number ", num, 
                      sep = ""))
          
          ## move to next window width
          n = n + 5
        }
        num = num + 1
      }
    }
    y = y + 1
  }
  x = x + 1
}


## bind rows in list into data frame
pop_exp_df <- data.frame(do.call(rbind, pop_exp_list), stringsAsFactors = FALSE)
colnames(pop_exp_df) <- c("r_spec_exp_PSD_low", "r_spec_exp_PSD_high", "r_spec_exp_PSD_all",
                          "K_spec_exp_PSD_low", "K_spec_exp_PSD_high", "K_spec_exp_PSD_all",
                          "K_cy_spec_exp_PSD_low", "K_cy_spec_exp_PSD_high", "K_cy_spec_exp_PSD_all",
                          "noise_spec_exp_PSD_low", "noise_spec_exp_PSD_high", "noise_spec_exp_PSD_all",
                          "demsto_spec_exp_PSD_low", "demsto_spec_exp_PSD_high","demsto_spec_exp_PSD_all",
                          "window_start_year",
                          "window_stop_year", "time_window_width", 
                          "icp", "lat", "lon", "sim")

#write.csv(pop_exp_df, "data-processed/BerkeleyEarth/forced-pop-dynams_chunk1_icp0.9_cy.csv", row.names = F)
pop_exp_df <- read.csv("data-processed/BerkeleyEarth/forced-pop-dynams_chunk1_icp0.1_cy.csv") %>%
 # rbind(., read.csv("data-processed/BerkeleyEarth/forced-pop-dynams_chunk1_icp0.9.csv")) %>%
  mutate(lat_lon = paste(lat, lon, sep = "_"),
         noise_spec_exp_PSD_low = as.numeric(as.character(noise_spec_exp_PSD_low)),
         r_spec_exp_PSD_low = as.numeric(as.character(r_spec_exp_PSD_low)),
         r_spec_exp_PSD_all = as.numeric(as.character(r_spec_exp_PSD_all)),
         noise_spec_exp_PSD_all = as.numeric(as.character(noise_spec_exp_PSD_all)))

pop_exp_df %>%
  ggplot(., aes(x = window_start_year, 
                y = noise_spec_exp_PSD_low)) + geom_point(colour = "black") +
  geom_smooth(method = "lm", se = F, colour = "black") +
  geom_point(aes(y = r_spec_exp_PSD_low, colour = sim), position = position_jitter(width = 0.05)) +
  geom_smooth(method = "lm", se = F, aes(y = r_spec_exp_PSD_low, colour = sim)) +
  facet_wrap(icp~time_window_width) 

pop_exp_df %>%
  ggplot(., aes(x = window_start_year, 
                y = noise_spec_exp_PSD_low)) + geom_point(colour = "black") +
  geom_smooth(method = "lm", se = F, colour = "black") +
  geom_point(aes(y = K_spec_exp_PSD_low, colour = sim), position = position_jitter(width = 0.05)) +
  geom_smooth(method = "lm", se = F, aes(y = K_spec_exp_PSD_low, colour = sim)) +
  facet_wrap(icp~time_window_width) 

## facet by unique grid cell location 
pop_exp_df %>%
  filter(lat_lon %in% unique(pop_exp_df$lat_lon)[1:20]) %>%
  filter(time_window_width == "3650 days") %>%
  ggplot(., aes(x = window_start_year, 
                y = noise_spec_exp_PSD_low)) + geom_point(colour = "black") +
  geom_smooth(method = "lm", se = F, colour = "black") +
  geom_point(aes(y = K_spec_exp_PSD_low, colour = sim), position = position_jitter(width = 0.05)) +
  geom_smooth(method = "lm", se = F, aes(y = K_spec_exp_PSD_low, colour = sim, group = sim)) +
  facet_wrap(icp~lat_lon) 

## regress and plot correlation between slopes
## for population dynamics without forcing:
pop_exp_df$window_start_year <- as.numeric(as.character(pop_exp_df$window_start_year))

noise <- pop_exp_df %>%
  mutate(lat_lon = paste(lat, lon, sep = "_")) %>%
  group_by(lat_lon, time_window_width, sim, icp) %>% # group by location, window width, and simulation
  do(tidy(lm(data = ., noise_spec_exp_PSD_all ~ window_start_year))) %>% # calculate slope 
  filter(term == "window_start_year") 

pops <- pop_exp_df %>%
  mutate(lat_lon = paste(lat, lon, sep = "_")) %>%
  filter(!is.na(demsto_spec_exp_PSD_all)) %>% 
  group_by(lat_lon, time_window_width, sim, icp) %>% # group by location, window width, and simulation
  do(tidy(lm(data = ., demsto_spec_exp_PSD_all ~ window_start_year))) %>% # calculate slope 
  filter(term == "window_start_year") 

## left join 
colnames(pops)[5:9] <- paste(colnames(pops)[5:9], "_popdynam", sep = "")
colnames(noise)[5:9] <- paste(colnames(noise)[5:9], "_noise", sep = "")

demsto <- left_join(pops, noise) 

## group by lat_lon, window width, icp and calculate mean slopes across all simulations
demsto <- demsto %>%
  group_by(lat_lon, time_window_width, icp) %>% 
  mutate(estimate_popdynam = mean(estimate_popdynam, na.rm=TRUE)) 

demsto$type = "demographic stochasticity"

## regress and plot correlation between slopes
## for r:
noise <- pop_exp_df %>%
  mutate(lat_lon = paste(lat, lon, sep = "_")) %>%
  group_by(lat_lon, time_window_width, sim, icp) %>% # group by location, window width, and simulation
  do(tidy(lm(data = ., noise_spec_exp_PSD_all ~ window_start_year))) %>% # calculate slope 
  filter(term == "window_start_year") 

pops <- pop_exp_df %>%
  mutate(lat_lon = paste(lat, lon, sep = "_")) %>%
  filter(!is.na(r_spec_exp_PSD_all)) %>% 
  group_by(lat_lon, time_window_width, sim, icp) %>% # group by location, window width, and simulation
  filter(n() >= 2) %>%
  do(tidy(lm(data = ., r_spec_exp_PSD_all ~ window_start_year))) %>% # calculate slope 
  filter(term == "window_start_year") 
  
## left join 
colnames(pops)[5:9] <- paste(colnames(pops)[5:9], "_popdynam", sep = "")
colnames(noise)[5:9] <- paste(colnames(noise)[5:9], "_noise", sep = "")
 
change_r <- left_join(pops, noise) 

## group by lat_lon, window width and calculate mean slopes across all simulations
change_r <- change_r %>%
  group_by(lat_lon, time_window_width, icp) %>% 
  mutate(estimate_popdynam = mean(estimate_popdynam, na.rm = TRUE)) 

change_r$type = "forced_r"

## regress and plot correlation between slopes
## for K:
noise <- pop_exp_df %>%
  mutate(lat_lon = paste(lat, lon, sep = "_")) %>%
  group_by(lat_lon, time_window_width, sim, icp) %>% # group by location, window width, and simulation
  do(tidy(lm(data = ., noise_spec_exp_PSD_all ~ window_start_year))) %>% # calculate slope 
  filter(term == "window_start_year") 

pops <- pop_exp_df %>%
  mutate(lat_lon = paste(lat, lon, sep = "_")) %>%
  filter(!is.na(K_spec_exp_PSD_all)) %>% 
  group_by(lat_lon, time_window_width, sim, icp) %>% # group by location, window width, and simulation
  filter(n() >= 2) %>%
  do(tidy(lm(data = ., K_spec_exp_PSD_all ~ window_start_year))) %>% # calculate slope 
  filter(term == "window_start_year") 

## left join 
colnames(pops)[5:9] <- paste(colnames(pops)[5:9], "_popdynam", sep = "")
colnames(noise)[5:9] <- paste(colnames(noise)[5:9], "_noise", sep = "")

change_K <- left_join(pops, noise) 

## group by lat_lon, window width and calculate mean slopes across all simulations
change_K <- change_K %>%
  group_by(lat_lon, time_window_width, icp) %>% 
  mutate(estimate_popdynam = mean(estimate_popdynam, na.rm = TRUE)) 

change_K$type = "forced_K"

## regress and plot correlation between slopes
## for K cy style:
noise <- pop_exp_df %>%
  mutate(lat_lon = paste(lat, lon, sep = "_")) %>%
  group_by(lat_lon, time_window_width, sim, icp) %>% # group by location, window width, and simulation
  do(tidy(lm(data = ., noise_spec_exp_PSD_all ~ window_start_year))) %>% # calculate slope 
  filter(term == "window_start_year") 

pops <- pop_exp_df %>%
  mutate(lat_lon = paste(lat, lon, sep = "_")) %>%
  filter(!is.na(K_cy_spec_exp_PSD_all)) %>% 
  group_by(lat_lon, time_window_width, sim, icp) %>% # group by location, window width, and simulation
  filter(n() >= 2) %>%
  do(tidy(lm(data = ., K_cy_spec_exp_PSD_all ~ window_start_year))) %>% # calculate slope 
  filter(term == "window_start_year") 

## left join 
colnames(pops)[5:9] <- paste(colnames(pops)[5:9], "_popdynam", sep = "")
colnames(noise)[5:9] <- paste(colnames(noise)[5:9], "_noise", sep = "")

change_K_cy <- left_join(pops, noise) 

## group by lat_lon, window width and calculate mean slopes across all simulations
change_K_cy <- change_K_cy %>%
  group_by(lat_lon, time_window_width, icp) %>% 
  mutate(estimate_popdynam = mean(estimate_popdynam, na.rm = TRUE)) 

change_K_cy$type = "forced_K_cy"

## combine:
change <- rbind(change_r, demsto) %>% rbind(., change_K)   %>% 
  rbind(., change_K_cy)

## plot correlation 
change %>%
  ggplot(aes(x = estimate_noise, y = estimate_popdynam, col = type)) + geom_point() +
  facet_wrap(icp~time_window_width) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) 

change %>%
  filter(type %in% c("forced_r", "forced_K", "forced_K_cy",
                                    "demographic stochasticity")) %>%
  ggplot(., aes(x = estimate_noise, y = estimate_popdynam, colour = icp)) + geom_point() +
  facet_wrap(~type) +
  labs(x = "Change in spectral exponent of temperature",
       y = "Change in spectral exponent of population size", colour = "ICP") +
  theme(panel.background = element_blank()) +
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) +
  scale_x_continuous(limits = c(-9e-8, 9e-6)) + 
  scale_y_continuous(limits = c(-5e-6, 5e-6)) 

## plot distribution of population size spectral exponents 
change %>%
  filter(type %in% c("forced_K_cy",
                     "demographic stochasticity"
                     #, "forced_r", "forced_K"
                     )) %>%
  ggplot(., aes(x = estimate_popdynam, fill = type)) + geom_histogram() +
  facet_wrap(~type) +
  geom_vline(xintercept = 0)


##########################################################
###             5) analyze extinction risk              ## 
##########################################################
## read in 100 population simulations for each location and calculate:
## 1) time to extinction (0, 10, 20, 30, 30, 50 individuals)
## 2) time to reach 50% population size
## 3) number of times 50% pop size is reached
## analyze change in noise colour over time windows 
icp = 1
lon_index = 0
lat_index = 0
l = seq(0.1, 1, by = 0.1)

lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
lon <-seq(from = 0.5, to = 359.5, length.out = 360) 

## get lat and lon bounds to extract in between:
lon_bound1 <- lon_index 
lon_bound2 <- lon_index + 60
lat_bound1 <- lat_index
lat_bound2 <- lat_index - 60

MTE_all <- data.frame()

## loop through grid cells
x = 1
while (x <= 360) {
  
  y = 1
  while (y <= 360) {
    
    filename = paste("data-processed/BerkeleyEarth/popdynamics_lat-", 
                     lat[y + lon_index], "_lon-", lon[x + lon_index],
                     "_icp-", l[icp], ".csv", sep = "")
    ## read in data:
    if (file.exists(filename)) {
      
      sims <- read.csv(filename)
    
      MTE_50 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 50)),
                  Nts_r = first(which(Nts_r <= 50)),
                  Nts_K = first(which(Nts_K <= 50)),
                  Nts_K_cy = first(which(Nts_K_cy <= 50)))
      
      MTE_50 <- gather(MTE_50, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_r, Nts_K_cy))
      MTE_50$extinction_metric <- "MTE_50"
      

      MTE_40 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 40)),
                  Nts_r = first(which(Nts_r <= 40)),
                  Nts_K = first(which(Nts_K <= 40)),
                  Nts_K_cy = first(which(Nts_K_cy <= 40)))
      
      MTE_40 <- gather(MTE_40, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_r, Nts_K_cy))
      MTE_40$extinction_metric <- "MTE_40"
      
      
      MTE_30 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 30)),
                  Nts_r = first(which(Nts_r <= 30)),
                  Nts_K = first(which(Nts_K <= 30)),
                  Nts_K_cy = first(which(Nts_K_cy <= 30)))
      
      MTE_30 <- gather(MTE_30, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_r, Nts_K_cy))
      MTE_30$extinction_metric <- "MTE_30"
      
      
      MTE_20 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 20)),
                  Nts_r = first(which(Nts_r <= 20)),
                  Nts_K = first(which(Nts_K <= 20)),
                  Nts_K_cy = first(which(Nts_K_cy <= 20)))
      
      MTE_20 <- gather(MTE_20, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_r, Nts_K_cy))
      MTE_20$extinction_metric <- "MTE_20"
      
      MTE_10 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 10)),
                  Nts_r = first(which(Nts_r <= 10)),
                  Nts_K = first(which(Nts_K <= 10)),
                  Nts_K_cy = first(which(Nts_K_cy <= 10)))
      
      MTE_10 <- gather(MTE_10, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_r, Nts_K_cy))
      MTE_10$extinction_metric <- "MTE_10"
      
      MTE_5 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 5)),
                  Nts_r = first(which(Nts_r <= 5)),
                  Nts_K = first(which(Nts_K <= 5)),
                  Nts_K_cy = first(which(Nts_K_cy <= 5)))
      
      MTE_5 <- gather(MTE_5, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_r, Nts_K_cy))
      MTE_5$extinction_metric <- "MTE_5"
      
      ## number of times 50 or fewer individuals is reached 
      n50 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = length(which(Nts <= 50)),
                  Nts_r = length(which(Nts_r <= 50)),
                  Nts_K = length(which(Nts_K <= 50)),
                  Nts_K_cy = length(which(Nts_K_cy <= 50)))
      
      n50 <- gather(n50, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_r, Nts_K_cy))
      n50$extinction_metric <- "n50"
      
      ## number of times 50 or fewer individuals is reached 
      n30 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = length(which(Nts <= 30)),
                  Nts_r = length(which(Nts_r <= 30)),
                  Nts_K = length(which(Nts_K <= 30)),
                  Nts_K_cy = length(which(Nts_K_cy <= 30)))
      
      n30 <- gather(n30, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_r, Nts_K_cy))
      n30$extinction_metric <- "n30"
  
      ## bind all:
      MTE <- rbind(MTE_10, MTE_20, MTE_30, MTE_40, MTE_50, MTE_5, n50, n30)
      
      ## add other info
      MTE$lat <- sims$lat[1]
      MTE$lon <- sims$lon[1]
      MTE$l <- sims$l[1]
      
      ## bind to rest of the data:
      MTE_all <- rbind(MTE_all, MTE)
      
      print(paste("On lat: ", lat[y + lon_index], " and lon: ", lon[x + lon_index],
                  sep = ""))
    }
    y = y + 1
  }
  x = x + 1
}


#write.csv(MTE_all, "data-processed/BerkeleyEarth/MTE_all_0.1.csv", row.names = F)
MTE_all <- read.csv("data-processed/BerkeleyEarth/MTE_all_0.1.csv")


## plot extinction risk
colnames(MTE_all)

MTE_all %>%
  ggplot(aes(x = extinction_time, fill = l)) + geom_histogram() + 
  facet_wrap(~extinction_metric)


MTE_all %>%
  filter(extinction_metric == "n50") %>%
  ggplot(aes(x = extinction_time)) + geom_histogram()  +
  facet_wrap(~pop_model)

## 30 individuals seems like good extinction metric for simplest forced model
MTE_all %>%
  filter(pop_model == "Nts_K_cy") %>%
  ggplot(aes(x = extinction_time)) + geom_histogram() +
  facet_wrap(~extinction_metric)

## next: plot MTE_30 (and other metrics) across baseline noise colour 
## bring in spectral analysis Berk Earth data 
spec <- read.csv("data-processed/BerkeleyEarth/forced-pop-dynams_chunk1_icp0.1_cy.csv") %>%
  filter(window_start_year == "0", time_window_width == "3650 days") %>%
  select(noise_spec_exp_PSD_low, noise_spec_exp_PSD_high, noise_spec_exp_PSD_all, lat, lon)

MTE_gr <- MTE_all %>%
  filter(pop_model == "Nts_K_cy") %>%
  select(-sim, -pop_model) %>%
  unique() %>%
  mutate(extinction_time = ifelse(is.na(extinction_time), 51830, extinction_time)) %>%
  group_by(extinction_metric, lat, lon) %>%
  mutate(extinction_time = mean(extinction_time, na.rm = T)) %>%
  unique()
  
ggplot(MTE_gr, aes(x = extinction_time)) + geom_histogram() +
  facet_wrap(~extinction_metric)

join <- left_join(MTE_gr, spec) %>%
  gather(key = "slope_type", value = "spectral_slope", c(noise_spec_exp_PSD_low, 
                                                         noise_spec_exp_PSD_high, noise_spec_exp_PSD_all))
join %>%
  filter(extinction_metric %in% c("n50")) %>%
  ggplot(aes(x = spectral_slope, y = extinction_time)) +
  geom_point() +
  facet_wrap(extinction_metric~slope_type)


## what about a change?

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
  
  high_freq <- filter(spectral, freq > 1/16)
  low_freq <- filter(spectral, freq < 1/16)
  
  # plot spectrum:
  # plot <- spectral %>%
  #   ggplot(aes(x = freq, y = power)) + geom_line() +
  #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm", colour = "red") +
  #   theme_light() +
  #   geom_smooth(data = high_freq, method = "lm", colour = "blue") +
  #   geom_smooth(data = low_freq, method = "lm", colour = "green")
  
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
  
  return(list(-model_output_low$estimate, -model_output_high$estimate, -model_output_all$estimate 
              ,plot
  ))
}

