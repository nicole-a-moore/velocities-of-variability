### replicating Cuddington and Yodzis population model results 
library(tidyverse)

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

colours = runif(50, 0, 3.2)

all_results <- list()
col = 1
l = 0
while (col <=  length(colours)) {
  print(paste("On spectral exponent:", col, sep = " "))
  x = 1 
  while (x <= 500) {
    N = N0 ## set starting population size to 100
    
    ## generate coloured noise: 
    output <- one_over_f(beta = colours[col])
    noise = output[[1]]
    
    K <- K0 + round(noise, digits = 0) ## generate new carrying capacity
    K[which(K <= 0)] <- 5 ## make sure none are negative
    
    
    i=1
    while(i <= Tmax) {
        Nt = N*exp(1.5*(1 - (N/K[i])^l))
        Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
        
        ## save if the pop goes extinct:
        if (is.na(Nt) | Nt == 0 | i == Tmax) {
          if (i == 1 & x == 1) {
            results <- data.frame(run = x, N = N0, K = K[i], l = l, t = i, colour = colours[col],
                                  true_colour = output[[2]])
            i = i + 1
            
          }
          else {
            results <- rbind(results, data.frame(run = x, N = N0, K = K[i], l = l, t = i,
                                                 colour = colours[col], true_colour = output[[2]]))
          }
          
          i = Tmax+1
        }
        ## otherwise continue
        else {
          
          N = Nt   
          
          i = i + 1
        }
      
    }
    
    x = x+1
  }
  
  all_results[[col]] <- results
  results <- c()
  col = col + 1
}

saveRDS(all_results, "data-processed/ricker_all-results_l0.rds")

all <- do.call("rbind", all_results) 

all %>%
  ggplot(., aes(x = -true_colour, y = log(t), colour = l)) + geom_point()

## next l
l = 1
all_results <- list()
col = 1
while (col <=  length(colours)) {
  print(paste("On spectral exponent:", col, sep = " "))
  x = 1 
  while (x <= 500) {
    N = N0 ## set starting population size to 100
    
    ## generate coloured noise: 
    output <- one_over_f(beta = colours[col])
    noise = output[[1]]
    
    K <- K0 + round(noise, digits = 0) ## generate new carrying capacity
    K[which(K <= 0)] <- 5 ## make sure none are negative
    
    
    i=1
    while(i <= Tmax) {
      Nt = N*exp(1.5*(1 - (N/K[i])^l))
      Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
      
      ## save if the pop goes extinct:
      if (is.na(Nt) | Nt == 0 | i == Tmax) {
        if (i == 1 & x == 1) {
          results <- data.frame(run = x, N = N0, K = K[i], l = l, t = i, colour = colours[col],
                                true_colour = output[[2]])
          i = i + 1
          
        }
        else {
          results <- rbind(results, data.frame(run = x, N = Nt, K = K[i], l = l, t = i,
                                               colour = colours[col], true_colour = output[[2]]))
        }
        
        i = Tmax+1
      }
      ## otherwise continue
      else {
       
        N = Nt   
        
        i = i + 1
      }
      
    }
    
    x = x+1
  }
  
  all_results[[col]] <- results
  results <- c()
  col = col + 1
}

saveRDS(all_results, "data-processed/ricker_all-results_l1.rds")

l1 <- readRDS("data-processed/ricker_all-results_l1.rds")
l0 <- readRDS("data-processed/ricker_all-results_l0.rds")

all_l1 <- do.call("rbind", l1) 
all_l0 <- do.call("rbind", l0) 

all <- rbind(all_l1, all_l0)
  
all <- filter(all, t > 1)

all %>%
  ggplot(., aes(x = -true_colour, y = log(t), colour = l)) + geom_point()

all %>%
  group_by(colour, l) %>%
  mutate(mean_pers_time = mean(t))%>%
  select(-true_colour) %>%
  unique(.) %>%
  ggplot(., aes(y = mean_pers_time, x = colour, colour = l)) + geom_point() +
  scale_y_log10()

all %>%
  group_by(colour, l) %>%
  mutate(extinction_risk  = length(which(t <= 100))/length(t)*100) %>%
  select(-true_colour) %>%
  unique(.) %>%
  ggplot(., aes(y = extinction_risk, x = colour, colour = l)) + geom_point() 

all %>%
  group_by(colour, l) %>%
  mutate(cv = sd(t)/mean(t)) %>%
  ungroup() %>%
  select(-t) %>%
  unique(.) %>%
  ggplot(., aes(y = cv, x = colour, colour = l)) + geom_point() 



### garbage:
## check that exponent matches
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
