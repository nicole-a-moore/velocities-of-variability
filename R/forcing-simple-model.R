## force population models with Berkeley Earth temperature time series 
## and with baseline-scenario time series
library(tidyverse)
library(foreach)
library(doParallel)
library(broom)
library(ncdf4)
theme_set(theme_bw())

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
cl <- makePSOCKcluster(numCores)

registerDoParallel(numCores)

## get paths and filenames 
filenames <- readRDS("data-processed/BerkeleyEarth/be_sp_files.rds")

##############################################
###             basic C & Y model           ## 
##############################################
## write population model:
popmod <- function(t, # number of time steps
                   N0, # initial population size
                   K # carrying capacity 
                   ){
  N <- numeric(t)
  N[1] <- N0
  for(i in 2:t){
    N[i] <- sample(rpois(as.numeric(N[i-1]*exp(1.5*(1 - (N[i-1]/K[i])^lambda[icp]))), n = 1000), 
                   size = 1)
  }
  return(N)
}

#### model parameters:
K0 = 100  #carrying capacity
N0 = 100  #starting population size
lambda = c(0.1, 0.5, 0.9) # intraspecific competition parameter 

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
      
      local_ts <- tas[y,]
    
      ## if time series is not na, force the population models 
      if(length(which(is.na(local_ts))) != length(local_ts)) {
        
        ## create fake local time series by repeating first 5 years of non-NA temperatures 
        fake_ts <- rep(local_ts[1:1825], 29)[1:51830]
        
        ## calculate new carrying capacity 
        K <- K0 + (local_ts - 273.15) 
        K[which(K <= 0)] <- 5 ## make sure none are negative
        K <- round(K, 0) ## round to whole individuals
        
        K_nc <- K0 + (fake_ts - 273.15) 
        K_nc[which(K_nc <= 0)] <- 5 ## make sure none are negative
        K_nc <- round(K_nc, 0) ## round to whole individuals
        
        ## loop through ICPs
        icp = 1 
        if (!file.exists(paste("data-processed/BerkeleyEarth/newpopdynamics_lat-", lat[y + lat_index], "_lon-", 
                               lon[x + lon_index],"_icp-", lambda[icp], ".csv", sep = ""))) {
          
          while (icp < length(lambda)+1) {
            print(paste("On lat: ", y, " and lon: ", x, " and icp ", icp,  sep = ""))
            
            ## run 100 simulations per model:
            all <- foreach (z = 1:100, .combine=rbind)  %dopar% {
              
              Nts = popmod(t = length(local_ts), N0 = N0, K = rep(100, length(local_ts))) # no forcing
              Nts_K = popmod(t = length(local_ts), N0 = N0, K = K) # C & Y model - change
              Nts_K_nc = popmod(t = length(local_ts), N0 = N0, K = K_nc) # C & Y model - no change
              
              ## return population time series
              data.frame(lat = lat[y + lat_index], lon = lon[x + lon_index], t = 1:51830,
                         sim = z, K = K, Nts = Nts, Nts_K = Nts_K, Nts_K_nc, lambda = lambda[icp], 
                         noise = local_ts - 273.15, fake_ts = fake_ts - 273.15)
              
              
            }
            
            ## save the data for location 
            write.csv(all, paste("data-processed/BerkeleyEarth/newpopdynamics_lat-", lat[y + lat_index], "_lon-", 
                                 lon[x + lon_index],"_icp-", lambda[icp], ".csv", sep = ""), 
                      row.names = F)
            
            icp = icp + 1
          }
        }
        else {
          icp = icp + 1
        }
      }
      y = y + 1
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


## plot some time series
all %>%
  filter(sim %in% 1:10) %>%
  gather(key = "data_type", value = "value", c(Nts, Nts_K, Nts_K_nc, K, noise)) %>%
  filter(data_type %in% c("Nts_K_nc", "noise")) %>%
  ggplot(aes(x = t, y = value, colour = data_type)) + geom_line() + 
  facet_wrap(~sim)
  


##########################################################
###             5) analyze extinction risk              ## 
##########################################################
## read in 100 population simulations for each location and calculate:
## 1) time to extinction (0, 10, 20, 30, 30, 50 individuals)
## 2) time to reach 50% population size
## analyze change in noise colour over time windows 
icp = 1
lon_index = 0
lat_index = 0
l = c(0.1, 0.5, 0.9)

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
    
    filename = paste("data-processed/BerkeleyEarth/newpopdynamics_lat-", 
                     lat[y + lon_index], "_lon-", lon[x + lon_index],
                     "_icp-", l[icp], ".csv", sep = "")
    ## read in data:
    if (file.exists(filename)) {
      
      sims <- read.csv(filename)
      
      MTE_50 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 50)),
                  Nts_K = first(which(Nts_K <= 50)),
                  Nts_K_nc = first(which(Nts_K_nc <= 50)))
      
      MTE_50 <- gather(MTE_50, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_K_nc))
      MTE_50$extinction_metric <- "MTE_50"
      
      
      MTE_40 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 40)),
                  Nts_K = first(which(Nts_K <= 40)),
                  Nts_K_nc = first(which(Nts_K_nc <= 40)))
      
      MTE_40 <- gather(MTE_40, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_K_nc))
      MTE_40$extinction_metric <- "MTE_40"
      
      
      MTE_30 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 30)),
                  Nts_K = first(which(Nts_K <= 30)),
                  Nts_K_nc = first(which(Nts_K_nc <= 30)))
      
      MTE_30 <- gather(MTE_30, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_K_nc))
      MTE_30$extinction_metric <- "MTE_30"
      
      
      MTE_20 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 20)),
                  Nts_K = first(which(Nts_K <= 20)),
                  Nts_K_nc = first(which(Nts_K_nc <= 20)))
      
      MTE_20 <- gather(MTE_20, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_K_nc))
      MTE_20$extinction_metric <- "MTE_20"
      
      MTE_10 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 10)),
                  Nts_K = first(which(Nts_K <= 10)),
                  Nts_K_nc = first(which(Nts_K_nc <= 10)))
      
      MTE_10 <- gather(MTE_10, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_K_nc))
      MTE_10$extinction_metric <- "MTE_10"
      
      MTE_5 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = first(which(Nts <= 10)),
                  Nts_K = first(which(Nts_K <= 10)),
                  Nts_K_nc = first(which(Nts_K_nc <= 10)))
      
      MTE_5 <- gather(MTE_5, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_K_nc))
      MTE_5$extinction_metric <- "MTE_5"
      
      ## number of times 50 or fewer individuals is reached 
      n50 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = length(which(Nts <= 50)),
                  Nts_K = length(which(Nts_K <= 50)),
                  Nts_K_nc = length(which(Nts_K_nc <= 50)))
      
      n50 <- gather(n50, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_K_nc))
      n50$extinction_metric <- "n50"
      
      ## number of times 50 or fewer individuals is reached 
      n30 <- sims %>%
        group_by(sim) %>%
        summarise(Nts = length(which(Nts <= 30)),
                  Nts_K = length(which(Nts_K <= 30)),
                  Nts_K_nc = length(which(Nts_K_nc <= 30)))
      
      n30 <- gather(n30, key = "pop_model", value = "extinction_time", c(Nts, Nts_K, Nts_K_nc))
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

## save:
# write.csv(MTE_all, "data-processed/BerkeleyEarth/simulations/newpopdynamics_extinction-metrics_0.9.csv", row.names = F)

MTE_all %>%
  ggplot(aes(x = extinction_time, fill = pop_model)) + geom_histogram() + 
  facet_wrap(~extinction_metric)

MTE_all %>%
  filter(extinction_metric == "n50") %>%
  ggplot(aes(x = extinction_time, fill = pop_model)) + geom_histogram(position = position_dodge()) 
## yay! distributions differ

## calculate average extinction risk across the 100 simulations for each lat x lon, metric, and pop model combo
## if population never reaches extinction for MTE, assign MTE = length of time series + 1
MTE_avg <- MTE_all %>%
  mutate(extinction_time = ifelse(is.na(extinction_time), 51831, extinction_time)) %>%
  group_by(lat, lon, pop_model, extinction_metric) %>%
  summarise(extinction_time = mean(extinction_time))

MTE_avg %>%
  filter(extinction_metric == "MTE_50") %>%
  ggplot(aes(x = extinction_time, fill = pop_model)) + geom_histogram(position = position_dodge()) 

## compare Nts_K results to Nts_K_cy results
## plot Nts_K vs. Nts_K_cy and see how different for each extinction metric
MTE_avg %>%
  filter(extinction_metric == "MTE_50") %>%
  spread(key = pop_model, value = extinction_time) %>%
  ggplot(aes(x = Nts_K, y = Nts_K_nc)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(~extinction_metric)

## plot colour of first time window x extinction time with shape = pop_model
## bring in colour of each lat x lon
col_data <- readRDS("data-processed/BerkeleyEarth/BE_noise-colour.rds")

MTE_col <- left_join(MTE_avg, col_data)

MTE_col %>%
  filter(pop_model != "Nts") %>%
  filter(extinction_metric == "MTE_50") %>%
  ggplot(aes(x = s_exp_PSD_low_all, y = extinction_time, col = pop_model)) +
  geom_point()

## if the population models do not have different results, dead end  

## see if magnitude of difference in extinction risk can be explained by change in spectral exponent 
MTE_col %>%
  filter(pop_model != "Nts") %>%
  spread(key = pop_model, value = extinction_time) %>%
  mutate(extinction_diff = Nts_K_nc - Nts_K) %>%
  filter(extinction_metric == "MTE_50") %>% 
  ggplot(aes(x = s_estimate_PSD_low, y = extinction_diff)) +
  geom_point()

