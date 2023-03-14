## making a plot to show how we simulated and corrected noise with changing autocorrelation
library(tidyverse)
library(ggplot2)
library(foreach)
library(doParallel)
library(splus2R)
library(broom)
#library(ifultools)
source("R/population-simulations/fractal_functions.R")

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
  
  # ## plot spectrum:
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

## function to calculate spectral exponent over a time series window uisng average wavelet coefficient method
spectral_exponent_calculator_AWC <- function(ts_window, N) {
  ## a. compute wavelet transform
  wavelets <- biwavelet::wt(data.frame(time = 1:N, val = ts_window), do.sig = F)
  
  ## b. calculate arithmetic mean with respect to the translation coefficient (b)
  data <- data.frame(avg_wavelet_coeff = rowMeans(sqrt(wavelets$power), na.rm = TRUE), period = wavelets$period)
  
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

## write some parameters:
start_colour = 1.5
end_colour = 1.6
num_steps = 20
Tmax = 200000

by = (end_colour - start_colour)/num_steps
each = ceiling(Tmax/num_steps)

## simulate the time series
alpha <- seq(start_colour, end_colour, by = by)
alpha_inc = rep(alpha, each = each)[1:Tmax]  
alpha_stable = rep(c(start_colour, start_colour + 0.000001), num_steps)
alpha_stable = rep(alpha_stable, each = each)[1:200000]

## plot the different spectral exponents 
alpha_dec = rev(alpha_inc) - 0.095
plot(ts(alpha))
plot(ts(alpha_inc))
plot(ts(alpha_dec))

## calculate delta
delta_inc = -alpha_inc/-2
delta_dec = rev(delta_inc) - 0.095
delta_stable = -alpha_stable/-2

## set the innovations variance to unity
innovation <- rep(1, Tmax)

## simulate a time-varying FD process
noise_inc <- FDSimulate(delta = delta_inc, innovation = innovation)
noise_dec <- FDSimulate(delta = delta_dec, innovation = innovation)
noise_stable <- FDSimulate(delta = delta_stable, innovation = innovation)

## plot the different noises
plot(ts(noise_inc))
plot(ts(noise_dec))
plot(ts(noise_stable))

## if time series crosses from alpha < 1 to alpha > 1
if (any(alpha < 1) & !all(alpha < 1)) {
  
  ## define point at which alpha switches to greater/less than 1
  point_inc = first(which(alpha >= 1))
  point_dec = last(which(rev(alpha) >= 1))
  
  ## correct data so that each time step begins at last time step value IF alpha < 1:
  breaks = seq(1, 220000, by = each)
  new_noise_inc = c(noise_inc)
  new_noise_dec = c(noise_dec)
  new_noise_stable = c(noise_stable)
  
  i = 2
  while(i < length(breaks) - 1) {
    if(i > point_inc) {
      to_add = new_noise_inc[breaks[i] - 1] - new_noise_inc[breaks[i]]
      
      new_noise_inc[breaks[i]:(breaks[i+1] - 1)] = new_noise_inc[breaks[i]:(breaks[i+1] - 1)] + to_add
    } 
    if (i < point_dec)  {
      to_add = new_noise_dec[breaks[i] - 1] - new_noise_dec[breaks[i]]
      
      new_noise_dec[breaks[i]:(breaks[i+1] - 1)] = new_noise_dec[breaks[i]:(breaks[i+1] - 1)] + to_add
    }
    i = i + 1
  }
} else if(any(alpha >= 1)) {
  
  ## fix it so that each time step begins at last time step value:
  breaks = seq(1, 220000, by = each)
  new_noise_inc = c(noise_inc)
  new_noise_dec = c(noise_dec)
  new_noise_stable = c(noise_stable)
  
  i = 2
  while(i < length(breaks)) {
    ## increasing
    to_add = new_noise_inc[breaks[i] - 1] - new_noise_inc[breaks[i]]
    
    new_noise_inc[breaks[i]:(breaks[i+1] - 1)] = new_noise_inc[breaks[i]:(breaks[i+1] - 1)] + to_add
    
    ## decreasing
    to_add = new_noise_dec[breaks[i] - 1] - new_noise_dec[breaks[i]]
    
    new_noise_dec[breaks[i]:(breaks[i+1] - 1)] = new_noise_dec[breaks[i]:(breaks[i+1] - 1)] + to_add
    
    ## stable
    to_add = new_noise_stable[breaks[i] - 1] - new_noise_stable[breaks[i]]
    
    new_noise_stable[breaks[i]:(breaks[i+1] - 1)] = new_noise_stable[breaks[i]:(breaks[i+1] - 1)] + to_add
    
    i = i + 1
  }
} else {
  new_noise_inc = c(noise_inc)
  new_noise_dec = c(noise_dec)
  new_noise_stable = c(noise_stable)
}

list = list(new_noise_inc[1:200000], new_noise_dec[1:200000], new_noise_stable[1:200000])

## remove mean, standardize variance
new_list = lapply(list, 
                  FUN = function(x) {
  ## remove mean and change variance to 4140:
  n <- x*1/sqrt(var(x, na.rm = TRUE))*sqrt(4140)
  n <- n - mean(n,  na.rm = TRUE)
  return(n)
})
  

## plot the new noise 
plot(ts(new_list[[1]])) # inc 
plot(ts(new_list[[2]])) # dec
plot(ts(new_list[[3]])) # stable

## now calculate spectral exponent over time after the fact for each time series 
spec_exp_df_all <- c()
types = c("Increasing autocorrelation", "Decreasing autocorrelation", "Stable autocorrelation")
for(i in 1:3) {
  local_ts = new_list[[i]]
  
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
  spec_exp_df$time_series_type = types[i]
  
  ## convert numbers to numeric
  spec_exp_df[,1:4] <- sapply(spec_exp_df[,1:4], as.numeric)
  
  spec_exp_df_all <- rbind(spec_exp_df, spec_exp_df_all)
}





## now plot them all
## one panel: showing stepwise increase / decrease 
df <- data.frame(time = 1:200000, 
                 stable_noise = alpha_stable, 
                 alpha_inc = alpha_inc,
                 alpha_dec = alpha_dec) %>%
  gather(key = "ts_type", value = "spectral_exponent", 
         c(stable_noise, alpha_inc, alpha_dec))

panel_a <- df %>%
  ggplot(., aes(x = time, y = spectral_exponent, colour = ts_type)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("#BA3B46", "#3FA34D", "black"),
                     labels = c("Decreasing autocorrelation", 
                                "Increasing autocorrelation",
                                "Stable autocorrelation")) +
  labs(colour = "", x = "Time", y = '\nSpectral exponent') 

ggsave(panel_a, path = "figures/noise-simulation", filename = "panel_a.png", 
        "png", width = 10, height = 2)

## now plot simulated noise 
df <- data.frame(time = 1:200000, 
                 noise_stable = c(noise_stable), 
                 noise_inc = c(noise_inc),
                 noise_dec = c(noise_dec), 
                 new_noise_stable = c(new_list[[3]]), 
                 new_noise_inc = c(new_list[[1]]),
                 new_noise_dec = c(new_list[[2]])) %>%
  gather(key = "ts_type", value = "ts_value", 
         c(noise_stable, noise_inc, noise_dec,
           new_noise_stable, new_noise_inc, new_noise_dec)) %>%
  mutate(corrected_or_uncorrected = ifelse(str_detect(ts_type, "new"), "Corrected", 
                                                      "Uncorrected"))
## panel showing uncorrected noise
panel_b <- df %>%
  filter(corrected_or_uncorrected == "Uncorrected") %>%
  ggplot(., aes(x = time, y = ts_value, colour = ts_type)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("#BA3B46", "#3FA34D", "black"),
                     labels = c("Decreasing autocorrelation", 
                                "Increasing autocorrelation",
                                "Stable autocorrelation")) +
  labs(colour = "", x = "Time", y = '\nSimulated value') 

ggsave(panel_b, path = "figures/noise-simulation", filename = "panel_b.png", 
       "png", width = 10, height = 2)

## another showing corrected noise with adjusted mean and variance 
panel_c <- df %>%
  filter(corrected_or_uncorrected == "Corrected") %>%
  ggplot(., aes(x = time, y = ts_value, colour = ts_type)) +
  geom_line() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("#BA3B46", "#3FA34D", "black"),
                     labels = c("Decreasing autocorrelation", 
                                "Increasing autocorrelation",
                                "Stable autocorrelation")) +
  labs(colour = "", x = "Time", y = 'Corrected\nsimulated value') 

ggsave(panel_c, path = "figures/noise-simulation", filename = "panel_c.png", 
       "png", width = 10, height = 2)


## another showing spectral exponent of corrected versus uncorrected noise 
df <- data.frame(time = 1:200000, 
                 stable_noise = alpha_stable, 
                 alpha_inc = alpha_inc,
                 alpha_dec = alpha_dec) %>%
  gather(key = "ts_type", value = "spectral_exponent", 
         c(stable_noise, alpha_inc, alpha_dec)) %>%
  mutate(time_series_type = ifelse(ts_type == "stable_noise", "Stable autocorrelation", 
                                   ifelse(ts_type == "alpha_inc", "Increasing autocorrelation", 
                                          "Decreasing autocorrelation")))
 
panel_d <- spec_exp_df_all %>%
  filter(time_window_width == "5 years") %>%
  ggplot(aes(x = window_start_year, y = spec_exp_PSD_low, colour = time_series_type)) +
  geom_point() +
  labs(x = "time", y = "spectral exponent") +
  facet_wrap(~time_series_type, nrow = 3) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  scale_color_manual(values = c("#BA3B46", "#3FA34D", "black"),
                     labels = c("Decreasing autocorrelation", 
                                "Increasing autocorrelation",
                                "Stable autocorrelation")) +
  labs(colour = "", x = "Time", y = 'Spectral exponent') +
  geom_line(data = df, aes(x = time,  ## add lines showing what change in spectral exponent *should* look like
                           y = spectral_exponent,
                           colour = time_series_type),
            inherit.aes = FALSE, lwd = 0.5) + 
  geom_smooth(method = "lm", se = FALSE)
  
ggsave(panel_d, path = "figures/noise-simulation", filename = "panel_d.png", 
       "png", width = 10, height = 2)
 
## simulate 500 time series
## run spectral analysis on corrected versus uncorrected version 
## record slope estimates and plot distribution 
## write some parameters:
start_colour = 1.5
end_colour = 1.6
num_steps = 20
Tmax = 200000

by = (end_colour - start_colour)/num_steps
each = ceiling(Tmax/num_steps)

spec_exp_df_all <- c()
for(x in 1:100) {
  ## simulate the time series
  alpha <- seq(start_colour, end_colour, by = by)
  alpha_inc = rep(alpha, each = each)[1:Tmax]  
  alpha_stable = rep(c(start_colour, start_colour + 0.000001), num_steps)
  alpha_stable = rep(alpha_stable, each = each)[1:200000]
  
  ## calculate delta
  delta_inc = -alpha_inc/-2
  delta_dec = rev(delta_inc) - 0.095
  delta_stable = -alpha_stable/-2
  
  ## set the innovations variance to unity
  innovation <- rep(1, Tmax)
  
  ## simulate a time-varying FD process
  noise_inc <- FDSimulate(delta = delta_inc, innovation = innovation)
  noise_dec <- FDSimulate(delta = delta_dec, innovation = innovation)
  noise_stable <- FDSimulate(delta = delta_stable, innovation = innovation)
  
  ## if time series crosses from alpha < 1 to alpha > 1
  if (any(alpha < 1) & !all(alpha < 1)) {
    
    ## define point at which alpha switches to greater/less than 1
    point_inc = first(which(alpha >= 1))
    point_dec = last(which(rev(alpha) >= 1))
    
    ## correct data so that each time step begins at last time step value IF alpha < 1:
    breaks = seq(1, 220000, by = each)
    new_noise_inc = c(noise_inc)
    new_noise_dec = c(noise_dec)
    new_noise_stable = c(noise_stable)
    
    i = 2
    while(i < length(breaks) - 1) {
      if(i > point_inc) {
        to_add = new_noise_inc[breaks[i] - 1] - new_noise_inc[breaks[i]]
        
        new_noise_inc[breaks[i]:(breaks[i+1] - 1)] = new_noise_inc[breaks[i]:(breaks[i+1] - 1)] + to_add
      } 
      if (i < point_dec)  {
        to_add = new_noise_dec[breaks[i] - 1] - new_noise_dec[breaks[i]]
        
        new_noise_dec[breaks[i]:(breaks[i+1] - 1)] = new_noise_dec[breaks[i]:(breaks[i+1] - 1)] + to_add
      }
      i = i + 1
    }
  } else if(any(alpha >= 1)) {
    
    ## fix it so that each time step begins at last time step value:
    breaks = seq(1, 220000, by = each)
    new_noise_inc = c(noise_inc)
    new_noise_dec = c(noise_dec)
    new_noise_stable = c(noise_stable)
    
    i = 2
    while(i < length(breaks)) {
      ## increasing
      to_add = new_noise_inc[breaks[i] - 1] - new_noise_inc[breaks[i]]
      
      new_noise_inc[breaks[i]:(breaks[i+1] - 1)] = new_noise_inc[breaks[i]:(breaks[i+1] - 1)] + to_add
      
      ## decreasing
      to_add = new_noise_dec[breaks[i] - 1] - new_noise_dec[breaks[i]]
      
      new_noise_dec[breaks[i]:(breaks[i+1] - 1)] = new_noise_dec[breaks[i]:(breaks[i+1] - 1)] + to_add
      
      ## stable
      to_add = new_noise_stable[breaks[i] - 1] - new_noise_stable[breaks[i]]
      
      new_noise_stable[breaks[i]:(breaks[i+1] - 1)] = new_noise_stable[breaks[i]:(breaks[i+1] - 1)] + to_add
      
      i = i + 1
    }
  } else {
    new_noise_inc = c(noise_inc)
    new_noise_dec = c(noise_dec)
    new_noise_stable = c(noise_stable)
  }
  
  list = list(new_noise_inc[1:200000], new_noise_dec[1:200000], new_noise_stable[1:200000])
  
  ## remove mean, standardize variance
  new_list = lapply(list, 
                    FUN = function(x) {
                      ## remove mean and change variance to 4140:
                      n <- x*1/sqrt(var(x, na.rm = TRUE))*sqrt(4140)
                      n <- n - mean(n,  na.rm = TRUE)
                      return(n)
                    })
  
  ## now calculate spectral exponent over time after the fact for each time series 
  
  ## make dataframe of noise
  data = data.frame(noise_stable = noise_stable[1:200000], 
                    noise_inc = noise_inc[1:200000],
                    noise_dec = noise_dec[1:200000],
                    corr_noise_stable = new_list[[3]],
                    corr_noise_inc = new_list[[1]],
                    corr_noise_dec = new_list[[2]],
                    time = 1:200000)
  
  ## change NA to 0
  data[which(is.na(data))] <- 0
  element = 1
  spec_exp_list <- list()
  
  n = 5
  while (n < 11) {
    year_start <- 1
    year_stop <- 1 + n*365
    
    while (year_start <= (150000 - n*365)) {
      ## extract temps within time window
      ts_chunk <- data[year_start:year_stop,c(1,2,3,4,5,6)]
      
      ## get length
      L = nrow(ts_chunk)
      
      ## preprocess the time series:
      ## a. subtracting mean
      ts_ns <- ts_chunk$noise_stable - mean(ts_chunk$noise_stable)
      ts_ni <- ts_chunk$noise_inc - mean(ts_chunk$noise_inc)
      ts_nd <- ts_chunk$noise_dec - mean(ts_chunk$noise_dec)
      ts_cns <- ts_chunk$corr_noise_stable - mean(ts_chunk$corr_noise_stable)
      ts_cni <- ts_chunk$corr_noise_inc - mean(ts_chunk$corr_noise_inc)
      ts_cnd <- ts_chunk$corr_noise_dec - mean(ts_chunk$corr_noise_dec)
      
      
      ## b. windowing - multiply by a parabolic window 
      window <- parabolic_window(series = ts_ns, N = L)
      ts_ns <- ts_ns*window
      window <- parabolic_window(series = ts_ni, N = L)
      ts_ni <- ts_ni*window
      window <- parabolic_window(series = ts_nd, N = L)
      ts_nd <- ts_nd*window
      window <- parabolic_window(series = ts_cns, N = L)
      ts_cns <- ts_cns*window
      window <- parabolic_window(series = ts_cni, N = L)
      ts_cni <- ts_cni*window
      window <- parabolic_window(series = ts_cnd, N = L)
      ts_cnd <- ts_cnd*window
      
      ## c. bridge detrending (endmatching)
      ## ie. subtracting from the data the line connecting the first and last points of the series
      ts_ns <- bridge_detrender(windowed_series = ts_ns, N = L)
      ts_ni <- bridge_detrender(windowed_series = ts_ni, N = L)
      ts_nd <- bridge_detrender(windowed_series = ts_nd, N = L)
      ts_cni <- bridge_detrender(windowed_series = ts_cni, N = L)
      ts_cns <- bridge_detrender(windowed_series = ts_cns, N = L)
      ts_cnd <- bridge_detrender(windowed_series = ts_cnd, N = L)
      
      ## calculate spectral exponent in window using PSD and AWC methods
      exp_PSD_ns <- spectral_exponent_calculator_PSD(ts_ns, l = L)
      exp_PSD_ni <- spectral_exponent_calculator_PSD(ts_ni, l = L)
      exp_PSD_nd <- spectral_exponent_calculator_PSD(ts_nd, l = L)
      exp_PSD_cns <- spectral_exponent_calculator_PSD(ts_cns, l = L)
      exp_PSD_cni <- spectral_exponent_calculator_PSD(ts_cni, l = L)
      exp_PSD_cnd <- spectral_exponent_calculator_PSD(ts_cnd, l = L)
      
      exp_PSD_low_ns <- exp_PSD_ns[[1]]
      exp_PSD_low_ni <- exp_PSD_ni[[1]]
      exp_PSD_low_nd <- exp_PSD_nd[[1]]
      exp_PSD_low_cns <- exp_PSD_cns[[1]]
      exp_PSD_low_cni <- exp_PSD_cni[[1]]
      exp_PSD_low_cnd <- exp_PSD_cnd[[1]]
      
      spec_exp_list[[element]] <- c(exp_PSD_low_ns, exp_PSD_low_ni, exp_PSD_low_nd,
                                    exp_PSD_low_cns, exp_PSD_low_cni, exp_PSD_low_cnd,
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
  colnames(spec_exp_df) <- c("exp_PSD_low_ns", "exp_PSD_low_ni", "exp_PSD_low_nd",
                             "exp_PSD_low_cns", "exp_PSD_low_cni", "exp_PSD_low_cnd",
                             "window_start_year","window_stop_year", "time_window_width")
  
  spec_exp_df <- spec_exp_df %>%
    gather(key = "Time series type", value = "Measured spectral exponent",
           c(exp_PSD_low_ns, exp_PSD_low_ni, exp_PSD_low_nd, 
             exp_PSD_low_cns, exp_PSD_low_cni,exp_PSD_low_cnd)) %>%
    mutate(time_series_type = ifelse(`Time series type` %in% c("exp_PSD_low_cnd", "exp_PSD_low_nd"), 
                                                               "Decreasing autocorrelation",
                                                               ifelse(`Time series type` %in% c("exp_PSD_low_cni", "exp_PSD_low_ni"),
                                                                                                "Increasing autocorrelation", 
                                                                                                "Stable autocorrelation"))) %>%
    mutate(ts_type = ifelse(`Time series type` %in% c("exp_PSD_low_ns", "exp_PSD_low_ni", "exp_PSD_low_nd"),
                            "Original simulation",
                            ifelse(`Time series type` %in% c("exp_PSD_low_cns", "exp_PSD_low_cni", "exp_PSD_low_cnd"),
                                   "Corrected simulation",  NA)))
  

  ## now calculate change in spectral exponent and store the estimates 
  estimates <- spec_exp_df %>%
    mutate(window_start_year = as.numeric(as.character(window_start_year))) %>%
    group_by(time_window_width, `Time series type`, time_series_type, ts_type) %>%
    do(tidy(lm(`Measured spectral exponent` ~ window_start_year, data = .), conf.int = TRUE)) %>%
    filter(term == "window_start_year") %>%
    ungroup() 
  
  estimates$simulation_num = x
  
  spec_exp_df_all <- rbind(spec_exp_df_all, estimates)
}

#write.csv(spec_exp_df_all, "data-processed/after-the-fact-spec-exp.csv", row.names = FALSE)
spec_exp_df_all <- read.csv("data-processed/after-the-fact-spec-exp.csv")

panel_e <- spec_exp_df_all %>%
  mutate(time_window_width = factor(.$time_window_width, 
                                    levels = c("5 years", "6 years", "7 years", 
                                               "8 years", "9 years", "10 years"), 
                                    ordered = TRUE)) %>%
  mutate(ts_type = factor(.$ts_type, levels = c("Original simulation",
                                                "Corrected simulation"), 
                                    ordered = TRUE)) %>%
  mutate(expected_slope = ifelse(time_series_type == "Increasing autocorrelation", 0.1/200000,
                               ifelse(time_series_type == "Decreasing autocorrelation", -0.1/200000,
                                      0))) %>%
  filter(time_window_width == "5 years") %>%
  ggplot(., aes(x = ts_type, colour = time_series_type,
                y = expected_slope)) + 
  #facet_grid(time_window_width~time_series_type) +
  facet_wrap(~time_series_type, nrow = 1) + 
  geom_point(aes(y = estimate, x = ts_type, colour = time_series_type), inherit.aes = FALSE) +
  theme_bw() +
  geom_point(size = 4, shape = 3, colour = "black") +
  theme(panel.grid = element_blank(), legend.position = "right") +
  labs(x = "", colour = "", y = "Estimated change\nin spectral exponent") +
  scale_color_manual(values = c("#BA3B46", "#3FA34D", "black")) +
  scale_x_discrete(labels = c("Original\nsimulation", "Corrected\nsimulation"))

ggsave(panel_e, path = "figures/noise-simulation", filename = "panel_e.png", 
       "png", width = 10, height = 2)

  

