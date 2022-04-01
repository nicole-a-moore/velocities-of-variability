### script that takes temperature time series from Berkeley Earth and:
## 1) computes the variance of each time series
## 2) preps each time series for use in population model 
##    a) adjusts the variance of each time series to equal the avg. variance 
##    b) converts temperatures to degrees Kelvin
##    c) centres the mean of each time series around 273 K 
##    d) records the min and max temperature in each time series 
library(tidyverse)
library(raster)
library(ncdf4)

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
while(count < length(filenames)+1) {
  
  ## retrieve spatial chunk from nc file
  open = nc_open(filenames[count])
  tas = ncvar_get(open, "var1_1")
  nc_close(open)
  
  lat <- seq(from = 89.5, to = -89.5, length.out = 180) 
  lon <-seq(from = 0.5, to = 359.5, length.out = 360) 

  ## for each cell of the raster
  ts_list <- list()
  x = 1 ## longitude
  while (x < ncol(l_detrended_tas)+1) {
    y = 1 ## latitude
    ## parallelize:
    ts_list[[x]] <-  foreach(y=1:nrow(l_detrended_tas), combine=rbind)  %dopar% {
      
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
      
      c(min_temp = min, max_temp = max, lat =  lat[y + lat_index], lon = lon[x + lon_index])
    }
    
    x = x + 1
  }
  
  ## bind rows of data 
  ts_df <- data.frame(do.call(rbind, ts_list), stringsAsFactors = FALSE)
  
  ## add data to larger file:
  if (count == 1) {
    all_ts_df <- ts_df
  } 
  else {
    all_ts_df <- rbind(all_ts_df, ts_df)
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

## save data
write.csv(all_ts_df, "data-processed/BerkeleyEarth/be_minmax-data.csv", row.names = FALSE)



## parameterize the r model using minmax data 



## run population model simulations, forcing r, K and both r and K with temperature time series 






## make sure data looks good 
col <- seq(0, 3.2, by = 0.1)
names <- paste("data-processed/temperature-dependent-models_col-", 1:length(col), "_icp-0.1.rds", 
               sep = "")

data <- c()
i=1
while(i <= length(names)) {
  data <- rbind(data, read.csv(names[i]))
  print(i)
  i=i+1
}

## remake plot from c&y but for one ICP
sims <- filter(data, !sim %in% c(1:10)) %>%
  mutate(Nts = ifelse(is.na(Nts), 89000, Nts),
         Nts_r = ifelse(is.na(Nts_r), 89000, Nts_r),
         Nts_K = ifelse(is.na(Nts_K), 89000, Nts_K),
         Nts_rK = ifelse(is.na(Nts_rK), 89000, Nts_rK))

sims %>%
  ggplot(., aes(x = -true_colour, y = Nts_rK)) + geom_point()

sims %>%
  gather(key = "pop_model", value = "ext_time", c(Nts, Nts_r, Nts_K, Nts_rK)) %>%
  group_by(colour, pop_model) %>%
  mutate(mean_pers_time = mean(ext_time)) %>%
  select(-true_colour) %>%
  unique(.) %>%
  ggplot(., aes(y = mean_pers_time, x = colour, colour = pop_model)) + geom_point() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Noise colour", y = "Mean persistence time", colour = "Population model") + 
  facet_wrap(~pop_model)

sims %>%
  gather(key = "pop_model", value = "ext_time", c(Nts, Nts_r, Nts_K, Nts_rK)) %>%
  group_by(colour, pop_model) %>%
  mutate(cv = sd(ext_time)/mean(ext_time)) %>%
  ungroup() %>%
  select(-ext_time) %>%
  unique(.) %>%
  ggplot(., aes(y = cv, x = colour, colour = pop_model)) + geom_point()  +
  theme_light() +
  labs(x = "Noise colour", y = "CV persistence time", colour = "Population model") + 
  facet_wrap(~pop_model)

sims <- filter(data, sim == 2)

sims %>%
  filter(colour == 0.7) %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = Nts_K, group = sim)) + geom_line() +
  facet_wrap(~sim)





