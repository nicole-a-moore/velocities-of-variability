### replicating Cuddington and Yodzis population model results 
## using multiple cores!!
library(tidyverse)
library(parallel)
library(MASS)
library(foreach)
library(doParallel)

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
numCores

registerDoParallel(numCores)  # use multicore, set to the number of our cores

## demographic stochasticity: 
## draw the population size at a given time (Nt) from a Poisson distribution, Z, with the mean given by the expectation from the population model
# Nt= Z[Nt-1 exp(r(1- (Nt-1/Kt))^l)]
## N0 = K0 = 100
## r = 1.5

## environmental noise
## Kt = K0 + ot
## where ot was drawn from noise of certain colour
## where ot generated negative Kt value, carrying capacity was set to 5

## l = degree of overcompensatory/undercompensatory dynamics
## ranges 0-1, 1 indicates overcompensation

## extinction occurs when Nt < 1, and time to extinction t was recorded

## 500 simulations, run for a maximum duration of 200 000 generations

#### model parameters:
K0 = 100  #carrying capacity
Tmax = 200000 #length of time to run model
N0 = 100 #starting conditions

one_over_f <- function(beta){
  ## create 1/fB noise as described in Cuddington and Yodzis
  n = rnorm(524288, mean = 0, sd = 1) ## random numbers with 0 mean and unit variance 
  phases <- runif(524288, 0, 2*pi) ## random phases
  f = 1:524288 ## f
  
  a <- n*1/(f^(beta/2)) ### amplitudes = random normally distributed numbers * 1/f^beta/2
  
  complex <- a*cos(phases) + a*sin(phases) ## complex coefficients
  
  dft <- fft(complex, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
  noise = as.numeric(dft[1:200000]) ## crop the noise series to first 200,000 points
  
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

## parameters:
s = 5
d_j = 0.03 # juvenile mortality at reference temperature 
d_a = 0.05 # adult mortality rate at reference temperature 
Ad_j = 9000 # Arrhenius constant for juvenile mortality  
Ad_a = 10000 # Arrhenius constant for adult mortality  
Aa = -17000 # Arrhenius constant for development
alpha = 25 # height of term 2
b_tr = 50 # average per capita fecundity at reference temperature, height and shape of term 2
Topt = 275.15 # temperature at which avg fecundity is maximal, shifts term 2 
Tr = 268 # reference temperature, changes height of term curves

fitness <- function(temp) {
  TD = ((1/Tr) - (1/temp))
  bT = b_tr*exp(-((temp-Topt)^2/(2*(s^2))))
  inv_alphaT = alpha*exp(Aa*TD)
  term1 = d_a*exp(Ad_a*TD)
  term2 = (1/(alpha*exp(Aa*TD)))*(LambertW::W(b_tr*alpha*exp(Aa*TD - ((temp-Topt)^2/(2*(s^2))) + alpha*exp(Aa*TD)*(d_a*exp(Ad_a*TD) - d_j*exp(Ad_j*TD)))))
  rT = -term1 + term2
  
  return(c(term1, term2, rT, bT, inv_alphaT))
}

## plot:
range <- seq(from = -30, to = 40, by = 0.001)
temps <- c(273.15 + range)

## plot results 
data <- sapply(FUN = fitness, temps)
data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps - 273.15,
                  bT = data[4,], inv_alphaT = data[5,])

basic_plot <- function(data)  {
  bp <- data %>%
    select(-bT, inv_alphaT) %>%
    gather(key = "term", value = "value", c(term1, term2, rT)) %>%
    ggplot(., aes(x = temps, y = value, colour = term)) + geom_point(size = 0.1) + 
    labs(x = "Temperature (C)", y = "rm(T) or components of rm(T)") +
    scale_y_continuous(limits = c(0, 2))
  
  return(bp)
}

basic_plot(data) +
  geom_vline(xintercept = 13) +
  geom_vline(xintercept = -13)



colours = seq(0, 3.2, by = 0.1)
all_results <- list()
l = seq(0.1, 1, by = 0.1)


## approach: let maximum per capita rate of increase (r) vary as a function of T

icp = 1
while (icp <= length(l)) {
 
  col = 1
  all_col <- list()
  while (col <=  length(colours)) {
    print(paste("On colour ", col, " and icp ", icp, sep = ""))
    
    all <- foreach (1:500, .combine=rbind)  %dopar% {
      N = N0 ## set starting population size to 100
      
      ## generate coloured noise: 
      output <- one_over_f(beta = colours[col])
      noise = output[[1]]
      
      noise_K <- noise*1/sqrt(var(noise))*sqrt(8)
      K <- K0 + round(noise, digits = 0) ## generate new carrying capacity
      K[which(K <= 0)] <- 5 ## make sure none are negative
      K_scal <- K0 + round(noise_K, digits = 0) ## generate new carrying capacity
      K_scal[which(K_scal <= 0)] <- 5 ## make sure none are negative
      
      ## generation new growth rate r
      noise_temp <- noise*1/sqrt(var(noise))*sqrt(8) ## make noise less variable 
      rs <- sapply(FUN = fitness, noise_temp + 273.15)
      rmt = rs[3,]
      
      N_K = N_r = N_rK = N
      Nts <- Nts_rK <- Nts_r <- Nts_K <- N
      i=1
      while(i <= Tmax) {
       
        ## vary r, K, niether and both according to environment
        Nt = N*exp(1.5*(1 - (N/K0)^l[icp])) ## constant r, constant K
        Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
        
        Nt_r = N_r*exp(rmt[i]*(1 - (N_r/K0)^l[icp])) # variable r, constant K
        Nt_r = sample(rpois(as.numeric(Nt_r), n = 1000), size = 1)
        
        Nt_K = N_K*exp(1.5*(1 - (N_K/K[i])^l[icp])) # constant r, variable K
        Nt_K = sample(rpois(as.numeric(Nt_K), n = 1000), size = 1)
        
        Nt_rK = N_rK*exp(rmt[i]*(1 - (N_rK/K_scal[i])^l[icp])) # variable r, variable K
        Nt_rK = sample(rpois(as.numeric(Nt_rK), n = 1000), size = 1)
        
        N = Nt 
        N_r = Nt_r
        N_K = Nt_K
        N_rK = Nt_rK
        Nts <- append(Nts, Nt)
        Nts_rK <- append(Nts_rK, Nt_rK)
        Nts_K <- append(Nts_K, Nt_K)
        Nts_r <- append(Nts_r, Nt_r)
        i = i + 1
      }
      c(N = N0, K = K[extinction_time], l = l[icp], t = extinction_time, colour = colours[col],
        true_colour = output[[2]])
    }
    
    all_col[[col]] <- all
    results <- c()
    col = col + 1
  }
  
  ## bind results from icp together
  all_results[[icp]] <- do.call("rbind", all_col) 
  
  icp = icp + 1
}


plot(ts(Nts))
plot(ts(Nts_r))
plot(ts(Nts_K))
plot(ts(Nts_rK))

## estimate noise colour from a linear regression of pwoer spectrum:
L <- length(Nts)
dft <- fft(Nts)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()


## estimate noise colour from a linear regression of pwoer spectrum:
ext <- first(which(Nts_r == 0))
#L <- length(Nts_r[1:ext])
L <- length(Nts_r)
#dft <- fft(Nts_r[1:ext])/L
dft <- fft(Nts_r)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

## estimate noise colour from a linear regression of pwoer spectrum:
ext <- first(which(Nts_K == 0))
#L <- length(Nts_K[1:ext])
L <- length(Nts_K)
#dft <- fft(Nts_K[1:ext])/L
dft <- fft(Nts_K)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

## estimate noise colour from a linear regression of pwoer spectrum:
ext <- first(which(Nts_rK == 0))
#L <- length(Nts_rK[1:ext])
L <- length(Nts_rK)
#dft <- fft(Nts_rK[1:ext])/L
dft <- fft(Nts_rK)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()





