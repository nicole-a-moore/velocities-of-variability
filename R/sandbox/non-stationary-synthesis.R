### trying to make non-stationary time series using wavelets 
library(longmemo)
library(wsyn)
library(biwavelet)
library(tidyverse)

## generate a time series of coloured noise:
one_over_f <- function(beta){
  ## create 1/fB noise as described in Cuddington and Yodzis
  n = rnorm(524288, mean = 0, sd = 1) ## random numbers with 0 mean and unit variance 
  phases <- runif(524288, 0, 2*pi) ## random phases
  f = 1:524288 ## f
  
  a <- n*1/(f^(beta/2)) ### amplitudes = random normally distributed numbers * 1/f^beta/2
  
  complex <- a*cos(phases) + a*sin(phases) ## complex coefficients
  
  dft <- fft(complex, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
  noise = as.numeric(dft[1:2000]) ## crop the noise series to first 200,000 points
  
  ## remove mean and change variance to 4140:
  noise <- noise*1/sqrt(var(noise))*sqrt(4140)
  noise <- noise - mean(noise)
  
  ## estimate noise colour from a linear regression of pwoer spectrum:
  l <- length(noise)
  dft <- fft(noise)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral <- data.frame(freq = freq, power = amp^2)
  
  # spectral %>%
  #   ggplot(aes(x = freq, y = power)) + geom_line() +
  #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  #   theme_minimal()
  
  true_colour <- lm(data = spectral, log(power) ~ log(freq))
  
  return(list(noise, as.numeric(true_colour$coefficients[2])))
}

noise <- one_over_f(1)[[1]]
plot(ts(noise))

## estimate noise colour from a linear regression of pwoer spectrum:
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

lm(data = spectral, log(power) ~ log(freq))

## perform wavelet transform
wavelets <- biwavelet::wt(data.frame(time = 1:length(noise), val = noise), do.sig = F)
plot(wavelets, type = "wave")

## pull out wavelet coefficients
coeff <- wavelets$wave

## mess with it: '
## increase amplitude at long periods, decrease amplitude at short periods 
coeff_new <- coeff[,]
i=1
while (i <= nrow(coeff_new)) {
  coeff_new[i,] <- coeff_new[i,]*(1:ncol(coeff_new))*1.5
  i = i+1
}

# plot(ts(coeff[51,]))
# plot(ts(coeff_new[51,]))

wavelets$wave <- coeff_new
plot(wavelets, type = "wave")

ts <- colSums(coeff_new)
ts <- ts*1/sqrt(var(ts))*sqrt(4140)

plot(ts(ts))
plot(ts(noise))

## estimate noise colour from a linear regression of pwoer spectrum:
l <- length(ts)
dft <- fft(ts)/l
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral2 <- data.frame(freq = freq, power = amp^2)

spectral2 %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

lm(data = spectral2, log(power) ~ log(freq))





### try discrete wavelet transform
library(wavelets)
dwt = dwt(noise, boundary = "reflection")

w = dwt@W
plot(ts(w[[9]]))

ts = sum(w[[1:9]])
idwt <- idwt(dwt)

plot(ts(idwt))
plot(ts(noise))



