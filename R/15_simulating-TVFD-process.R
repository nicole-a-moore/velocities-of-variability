## using fractal package to simulate noise with changing spectral exponent 
#install.packages("fractal", repos="http://R-Forge.R-project.org")
library(fractal)
library(tidyverse)
library(broom)

## function of interest: FDSimulate


## Example:
## create a time-varying FD parameter, linearly 
## varying from white to pink noise, then jump 
## to a red noise plateau 
delta <- c(seq(0, 0.5, by=0.01), rep(1,100))
delta <- c(seq(0.1, 1, by=0.001), rep(1,100))

delta <- c(seq(0, 0.49, by = 0.0005))
alpha <- c(rep(0.1,100), rep(0.2,100), rep(0.3,100), 
           rep(0.4,100), rep(0.5,100), rep(0.6,100),
           rep(0.7,100), rep(0.8,100), rep(0.9,100),
           rep(1,100), rep(1.1,100), rep(1.2,100), rep(1.3,100), 
           rep(1.4,100), rep(1.5,100), rep(1.6,100),
           rep(1.7,100), rep(1.8,100), rep(1.9,100),
           rep(2,100))
delta = -alpha/-2

alpha = rep(1:300, each = 100)/100
delta = -alpha/-2
## set the innovations variance to unity 
innovation <- rep(1, length(delta))

## simulate a time-varying FD process 
z <- FDSimulate(delta=delta, innovation=innovation, method = "ce")
print(z)
plot(z)




## try one with smooth change in delta
start = Sys.time()
alpha <- seq(-1, -3, by = -0.0002)
alpha <- rep(-2.9, 1001)
alpha <- c(rep(-1, 5000), rep(-3, 5000))
delta = alpha/-2

## set the innovations variance to unity 
innovation <- rep(1, length(delta))

## simulate a time-varying FD process 
noise <- FDSimulate(delta=delta, innovation=innovation)
plot(noise)

end = Sys.time()
start - end


## try another
start = Sys.time()
alpha <- seq(-0.1, -1, by = -0.001)
delta = alpha/-2

## set the innovations variance to unity 
innovation <- rep(1, length(delta))

## simulate a time-varying FD process 
noise <- FDSimulate(delta=delta, innovation=innovation)
plot(noise)

## remove mean and change variance to 4140:
noise <- noise*1/sqrt(var(noise))*sqrt(4140)
noise <- noise - mean(noise)

plot(ts(noise))

## things to figure out: 
## how much to change spectral exponent?
## what baseline colour?
## how long to simulate a 200,000 step time series?
## try measuring spectral change on one of these

## try simulating smaller time series and piecing together 
alpha <- seq(-0.0000015, -3, by = -1.5e-5)
delta = alpha/-2

breaks = seq(0, 200000, by = 500)

start = Sys.time()
for (i in 1:400) {
  alpha_sub = alpha[(breaks[i]+1):(breaks[i+1])]
  delta_sub = alpha_sub/-2
  
  if (i == 1) {
    ## set the innovations variance to unity 
    innovation <- rep(1, length(delta_sub))
    
    noise <- FDSimulate(delta = delta_sub, innovation=innovation) 
  }
  else {
    noise <- append(noise,  FDSimulate(delta = delta_sub, innovation = innovation))
  }
}
end = Sys.time()
end - start

plot(ts(noise))
length(which(is.na(noise))) # should be 0

noise <- noise*1/sqrt(var(noise, na.rm = TRUE))*sqrt(4140)
noise <- noise - mean(noise, na.rm = TRUE)

plot(ts(noise))

## measure spectral change using sliding window method
## also measure variance 

#########################################
##        SENSITIVITY ANALYSIS:        ##
#########################################
## calculate spectral exponent using FFT over n year windows
## store spectral exponents and calculate slope
local_ts = noise

start = Sys.time()
alpha <- seq(-2.00001, -2.000011, by = -0.0000000002)
delta = alpha/-2

## set the innovations variance to unity 
innovation <- rep(1, length(delta))

## simulate a time-varying FD process 
noise <- FDSimulate(delta=delta, innovation=innovation)
plot(noise)

local_ts = noise

## change NA to 0
local_ts[which(is.na(local_ts))] <- 0
element = 1
spec_exp_list <- list()

n = 5
while (n < 11) {
  year_start <- 1
  year_stop <- 1 + n*365
  
  while (year_start <= (150000 - n*365)) {
    ## extract temps within time window
    ts_chunk <- local_ts[year_start:year_stop]
    
    ## get length
    L = length(ts_chunk)
    
    ## preprocess the time series:
    ## a. subtracting mean
    ts <- ts_chunk - mean(ts_chunk)
    
    ## b. windowing - multiply by a parabolic window 
    window <- parabolic_window(series = ts, N = L)
    ts <- ts*window
    
    ## c. bridge detrending (endmatching)
    ## ie. subtracting from the data the line connecting the first and last points of the series
    ts <- bridge_detrender(windowed_series = ts, N = L)
    
    ## calculate spectral exponent in window using PSD and AWC methods
    exp_PSD <- spectral_exponent_calculator_PSD(ts, l = L)
    
    exp_PSD_low <- exp_PSD[[1]]
    exp_PSD_high <- exp_PSD[[2]]
    exp_PSD_all <- exp_PSD[[3]]
    
    spec_exp_list[[element]] <- c(exp_PSD_low, exp_PSD_high,exp_PSD_all,
                                    year_start, year_stop, 
                                  paste(n, "years"))
    
    
    ## move to next window
    year_start = year_stop + 1
    year_stop = year_stop + n*365 
    
    element = element + 1
  }
  
  ## move to next window width
  n = n + 1
}

## bind rows in list into data frame
spec_exp_df <- data.frame(do.call(rbind, spec_exp_list), stringsAsFactors = FALSE)
colnames(spec_exp_df) <- c("spec_exp_PSD_low", "spec_exp_PSD_high", "spec_exp_PSD_all",
                           "window_start_year","window_stop_year", "time_window_width")

## convert numbers to numeric
spec_exp_df[,1:4] <- sapply(spec_exp_df[,1:4], as.numeric)

spec_exp_df %>%
  ggplot(aes(x = window_start_year, y = spec_exp_PSD_low, colour = time_window_width)) +
  geom_point() +
  labs(x = "time", y = "spectral exponent") +
  facet_wrap(~time_window_width)














