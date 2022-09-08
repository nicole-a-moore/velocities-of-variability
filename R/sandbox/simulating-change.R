## experimenting with wave alterations


## goal:
## take time series of simulated noise a la C & Y
## change the autocorrelation
## measure change 
beta = 1.5

# library(Rmpfr)
# vec <- 1 + 5*mpfr(1/c(1:51830), 50) # set arbitrary precision that's greater than R default
# vec <- append(vec, rep(1, 524288-51830))

## create 1/fB noise as described in Cuddington and Yodzis
n = rnorm(51830, mean = 0, sd = 1) ## random numbers with 0 mean and unit variance 
phases <- runif(51830, 0, 2*pi) ## random phases
f = 1:51830 

a <- n*1/(f^(beta/2)) ### amplitudes = random normally distributed numbers * 1/f^beta/2

complex <- a*cos(phases) + a*sin(phases) ## complex coefficients

dft <- fft(as.numeric(complex), inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
noise = as.numeric(dft[1:51830]) ## crop the noise series to first 200,000 points

## remove mean and change variance to 109.086:
noise <- noise*1/sqrt(var(noise))*sqrt(109.086)
noise <- noise - mean(noise)

## estimate noise colour from a linear regression of power spectrum:
l <- length(noise)
dft <- fft(noise)/l
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

true_colour <- lm(data = spectral, log(power) ~ log(freq))


## problem: performing inverse fourier transform removes the signal we added

## create data frame of each wave (frequency f, amplitude a) at discrete time points from 1:524288
## y=a*sin⁡(bx)
df <- data.frame(nrow = 524288, ncol = 51830)
freq = 1
t <- 1:51830 

f <- 1:(51830)/51830
while(freq <= 51830) {
  y = a[freq]*sin((f[freq])*t) + phases[freq]
  
  ## save in matrix where rows = different waves, columns = time 
  #df[freq, 1:51830] <- y
  
  ## add to vector:
  if(freq == 1) {
    sums <- y
  }
  else {
    sums <- sums + y
  }
  
  print(freq)
  freq = freq + 1
}

plot(ts(sums))
noise <- sums

## check colour: should be 1.5 
noise <- noise*1/sqrt(var(noise))*sqrt(109.086)

## subtract  mean
noise <- noise - mean(noise)

## estimate noise colour from a linear regression of power spectrum:
l <- length(noise)
dft <- fft(noise)/l
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)%>%
  filter(freq < 1/16)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

lm(data = spectral, log(power) ~ log(freq))


## check change in colour
## (this one should be stationary)
## for each time window, process time series and calculate colour
element <- 1
exp_list <- list()
n = 5
while (n < 11) {
  year_start <- 0
  year_stop <- n*365 
  
  while (year_start <= (length(noise) - n*365)) {
    ## extract temps within time window
    ts_chunk <- noise[year_start:year_stop]
    
    ## subtract mean
    ts <- ts_chunk - mean(ts_chunk)
    
    ## estimate noise colour from a linear regression of power spectrum:
    l <- length(ts)
    dft <- fft(ts)/l
    amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
    amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
    freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
    
    ## create periodogram data by squaring amplitude of FFT output
    spectral <- data.frame(freq = freq, power = amp^2)%>%
      filter(freq < 1/16)
    
    spectral %>%
      ggplot(aes(x = freq, y = power)) + geom_line() +
      scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
      theme_minimal()
    
    exp_PSD_low <- spectral %>%
      lm(., formula = log10(power) ~ log10(freq)) %>%
      tidy(.) %>%
      filter(term == "log10(freq)")
    
    ## move to next window
    year_start = year_stop + 1
    year_stop = year_stop + n*365 
    
    exp_list[[element]] <- c(-exp_PSD_low$estimate, year_start, year_stop, paste(n*365, "days"))
    element = element + 1
  }
  
  ## move to next window width
  n = n + 5
}

exp_df <- data.frame(do.call(rbind, exp_list), stringsAsFactors = FALSE)
colnames(exp_df) <- c("spec_exp_PSD_low", "window_start_year","window_stop_year", "time_window_width")

ggplot(data = exp_df, aes(x = as.numeric(as.character(window_start_year)),
                          y = as.numeric(as.character(spec_exp_PSD_low)))) + geom_point() +
  facet_wrap(~time_window_width) +
  geom_smooth(method = "lm")


## next: try changing amplitude over time
df <- data.frame(nrow = 51830, ncol = 51830)
freq = 1
t <- 1:51830 

f <- 1:(51830)/51830

## for first 100 frequencies, increase amplitude over time
while(freq <= 51830) {

  if(freq %in% 1:10000) {
    new_a <- a[freq]*rev(1/10:100000)[1:51830]
    y = new_a*sin((f[freq])*t) + phases[freq]
  }
  else {
    y = a[freq]*sin((f[freq])*t) + phases[freq]
  }
  
  ## save in matrix where rows = different waves, columns = time 
  #df[freq, 1:51830] <- y
  
  ## add to vector:
  if(freq == 1) {
    sums_change <- y
  }
  else {
    sums_change <- sums_change + y
  }
  
  print(freq)
  freq = freq + 1
}

plot(ts(sums_change))
noise_change <- sums_change

## check colour: should be 1.5 
noise_change <- noise_change*1/sqrt(var(noise_change))*sqrt(109.086)

## subtract  mean
noise_change <- noise_change - mean(noise_change)

## estimate noise colour from a linear regression of power spectrum:
l <- length(noise_change)
dft <- fft(noise_change)/l
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral_change <- data.frame(freq = freq, power = amp^2) %>%
  filter(freq < 1/16)

spectral_change %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

lm(data = spectral_change, log(power) ~ log(freq))



element <- 1
exp_list_change <- list()
n = 5
while (n < 11) {
  year_start <- 0
  year_stop <- n*365 
  
  while (year_start <= (length(noise_change) - n*365)) {
    ## extract temps within time window
    ts_chunk <- noise_change[year_start:year_stop]
    
    ## subtract mean
    ts <- ts_chunk - mean(ts_chunk)
    
    ## estimate noise colour from a linear regression of power spectrum:
    l <- length(ts)
    dft <- fft(ts)/l
    amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
    amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
    freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
    
    ## create periodogram data by squaring amplitude of FFT output
    spectral <- data.frame(freq = freq, power = amp^2)%>%
      filter(freq < 1/16)
    
    spectral %>%
      ggplot(aes(x = freq, y = power)) + geom_line() +
      scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
      theme_minimal()
    
    exp_PSD_low <- spectral %>%
      lm(., formula = log10(power) ~ log10(freq)) %>%
      tidy(.) %>%
      filter(term == "log10(freq)")
    
    ## move to next window
    year_start = year_stop + 1
    year_stop = year_stop + n*365 
    
    exp_list_change[[element]] <- c(-exp_PSD_low$estimate, year_start, year_stop, paste(n*365, "days"))
    element = element + 1
  }
  
  ## move to next window width
  n = n + 5
}

exp_df_change <- data.frame(do.call(rbind, exp_list_change), stringsAsFactors = FALSE)
colnames(exp_df_change) <- c("spec_exp_PSD_low", "window_start_year","window_stop_year", "time_window_width")

ggplot(data = exp_df_change, aes(x = as.numeric(as.character(window_start_year)),
                          y = as.numeric(as.character(spec_exp_PSD_low)))) + geom_point() +
  facet_wrap(~time_window_width) +
  geom_smooth(method = "lm")



## wavelet approach 
library(WaveletComp)
x1 <- periodic.series(start.period = 80, length = 1000)
x2 <- periodic.series(start.period = 30, length = 1000)
x <- x1 + x2 + 0.2*rnorm(1000)

my.data <- data.frame(x = x)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels") )

reconstruct(my.w, sel.period = 80, plot.waves = TRUE, lwd = c(1,2),
            legend.coords = "bottomleft")


## try analyzing Berkeley Earth TS with a lot of change 
## lat = 69.5, lon = 34.5
ts <- read.csv("data-processed/BerkeleyEarth/popdynamics_lat-69.5_lon-34.5_icp-0.1.csv") %>%
  filter(sim == 1)

ts <- ts$noise
plot(ts(ts))

my.data <- data.frame(x = ts)
my.w <- analyze.wavelet(my.data, "x",
                        loess.span = 0,
                        dt = 1, dj = 1/20,
                        lowerPeriod = 16,
                        upperPeriod = 128,
                        make.pval = TRUE, n.sim = 10)
wt.image(my.w, n.levels = 250,
         legend.params = list(lab = "wavelet power levels") )




