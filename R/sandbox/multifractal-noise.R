## mutltifractal spectral synthesis

one_over_f_x2 <- function(beta1, beta2){
# 524288
  n = rnorm(89000, mean = 0, sd = 1) ## random numbers with 0 mean and unit variance 
  n1 = n[1:100]
  n2 = n[1:89000]
  phases <- runif(89000, 0, 2*pi) ## random phases
  f1 = 1:100 ## f
  f2 = 1:89000 ## f
  
  a1 <- n1*1/(f1^(beta1/2)) ### amplitudes = random normally distributed numbers * 1/f^beta/2
  a2 <- n2*1/(f2^(beta2/2))
  
  a2 <- a2*1/sqrt(var(a2))*sqrt(var(a1))
  
  ## keep same overall power but shift 

  a = append(a1[1:100], a2[101:88900])
    
  complex <- a*cos(phases) + a*sin(phases) ## complex coefficients
  
  dft <- fft(complex, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
  noise = as.numeric(dft[1:89000]) ## crop the noise series to first 200,000 points
  
  ## remove mean and change variance to 4140:
  noise <- noise*1/sqrt(var(noise))*sqrt(4140)
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
    theme_minimal() +
    geom_vline(xintercept = 1/(100*10))
  
  true_colour1 <- filter(spectral, freq < 1/1000) %>%
    lm(data = ., log(power) ~ log(freq))
  true_colour2 <- filter(spectral, freq >= 1/1000) %>%
    lm(data = ., log(power) ~ log(freq))
  
  return(list(noise, as.numeric(true_colour$coefficients[2])))
}