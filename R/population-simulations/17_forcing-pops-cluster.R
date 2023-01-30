## script to run simulations on cluster 
## simulating population dynamics under changes in autocorrelation
.libPaths(c("~/projects/def-jsunday/nikkim/VoV/packages", .libPaths()))
library(tidyverse)
library(fractal)
library(foreach)
library(doParallel)

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- 10

registerDoParallel(numCores)  # use multicore, set to the number of our cores


# Plan:
#   Create 200k long time series, 100 of each type: 
#     Increase spectral exponent by 0.1 for each interval between 0 - 3 
#     Decrease spectral exponent by 0.1 for each interval between 3 - 1 
#     Stable spectral exponent for each baseline noise colour (0.1, 0.2... 3.0)
# 
# For each simulated time series:
#   Run simple C&Y population model for ICP = 0.1, 0.5, 0.9
#   Track time to extinction 
#   Store both simulated time series and population time series

#################################################
###                   FUNCTIONS                ## 
#################################################
## write population model:
popmod <- function(t, # number of time steps
                   N0, # initial population size
                   K){ # carrying capacity 
  N <- numeric(t)
  N[1] <- N0
  for(i in 2:t){
    N[i] <- sample(rpois(as.numeric(N[i-1]*exp(1.5*(1 - (N[i-1]/K[i])^lambda[icp]))), n = 1000), 
                   size = 1)
  }
  return(N)
}

## write model to simulate and correct tvfd noise:
tvfd <- function(start_colour, end_colour, num_steps, Tmax) {
  
  by = (end_colour - start_colour)/num_steps
  each = ceiling(Tmax/num_steps)
  
  ## simulate the time series
  alpha <- seq(start_colour, end_colour, by = by)
  alpha_inc = rep(alpha, each = each)[1:Tmax]
  alpha_stable = rep(start_colour, Tmax)
  
  delta_inc = -alpha_inc/-2
  delta_dec = rev(delta_inc)
  delta_stable = -alpha_stable/-2
  
  ## set the innovations variance to unity
  innovation <- rep(1, Tmax)
  
  ## simulate a time-varying FD process
  noise_inc <- FDSimulate(delta = delta_inc, innovation = innovation)
  noise_dec <- FDSimulate(delta = delta_dec, innovation = innovation)
  noise_stable <- FDSimulate(delta = delta_stable, innovation = innovation)
  
  ## if time series crosses from alpha < 1 to alpha > 1
  if(any(alpha < 1) & !all(alpha < 1)) {
    
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
  }
  
  else if (any(alpha > 1)) {
    
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
      
      i = i + 1
    }
  }
  else {
    new_noise_inc = c(noise_inc)
    new_noise_dec = c(noise_dec)
    new_noise_stable = c(noise_stable)
  }
  
  list = list(new_noise_inc[1:200000], new_noise_dec[1:200000], new_noise_stable[1:200000])
  
  ## remove mean, standardize variance
  new_list = lapply(list, FUN = function(x) {
    ## remove mean and change variance to 4140:
    n <- x*1/sqrt(var(x, na.rm = TRUE))*sqrt(4140)
    n <- n - mean(n,  na.rm = TRUE)
    return(n)
  })
  
  return(new_list)
}

#################################################
###                  SIMULATIONS               ## 
#################################################
base_col = seq(0, 3, by = 0.1)
all_results <- list()
lambda = c(0.1, 0.5, 0.9)

K0 = 100
N0 = 100
Tmax = 200000

col = 1
while (col <  length(base_col)) {
  
  start_colour = base_col[col]
  end_colour = base_col[col+1]
  
  ## run 100 simulations per model:
  icp = 1
  while (icp <= length(lambda)) {
    print(paste("On base colour ", col, " and icp ", icp, sep = ""))
    
    ## run 100 simulations per colour:
    all <- foreach (z = 1:100, .combine=rbind)  %dopar% {
      
      noise_stable = noise_dec = noise_inc = rep(NA, Tmax)
      
      ## create variable is_empty to keep track of whether one time series is empty or not
      is_empty = length(which(is.na(noise_inc))) == Tmax |
        length(which(is.na(noise_dec))) == 
        Tmax | length(which(is.na(noise_stable))) == Tmax 
      
      ## while at least one time series is empty 
      while(is_empty) {
        
        ## generate increasing, decreasing, and stable coloured noise: 
        output <- tvfd(start_colour = start_colour, 
                       end_colour = end_colour,
                       num_steps = 20,
                       Tmax = Tmax)
        
        noise_inc = output[[1]]
        noise_dec = output[[2]]
        noise_stable = output[[3]]
        
        ## check whether any are still empty (improperly simulated)
        is_empty = length(which(is.na(noise_inc))) == Tmax |
          length(which(is.na(noise_dec))) == Tmax | 
          length(which(is.na(noise_stable))) == Tmax 
      }
      
      ## force carrying capacity based on environment 
      K_stable <- K0 + round(noise_stable, digits = 0) 
      K_inc <- K0 + round(noise_inc, digits = 0) 
      K_dec <- K0 + round(noise_dec, digits = 0) 
      
      ## make sure none are negative
      K_stable[which(K_stable <= 0)] <- 5
      K_inc[which(K_inc <= 0)] <- 5
      K_dec[which(K_dec <= 0)] <- 5
      
      ## run pop models
      Nts_stable = popmod(t = Tmax, N0 = N0, K = K_stable) # C & Y model 
      Nts_inc = popmod(t = Tmax, N0 = N0, K = K_inc)
      Nts_dec = popmod(t = Tmax, N0 = N0, K = K_dec)
      
      data.frame(t = 1:Tmax, sim = z, 
                 N_stable = Nts_stable,
                 N_inc = Nts_inc,
                 N_dec = Nts_dec,
                 lambda = lambda[icp], 
                 noise_stable = noise_stable, 
                 noise_inc = noise_inc, 
                 noise_dec = noise_dec, 
                 base_col = start_colour)
    }
    
    ## write results of 100 sims for base colour x icp combo
    write.csv(all, paste("/Volumes/NIKKI/pop-sims/tvfdpopdynam_icp-", lambda[icp], "_base-col-", 
                         start_colour, "_20-steps.csv", sep = ""), row.names = F)
    
    icp = icp + 1
  }
  ## move to next colour
  col = col + 1
}
  