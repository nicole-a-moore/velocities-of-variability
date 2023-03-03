## script to run marine simulations on cluster 
## simulating population dynamics under changes in autocorrelation
.libPaths(c("~/projects/def-jsunday/nikkim/VoV/packages", .libPaths()))
library(tidyverse)
library(foreach)
library(doParallel)
library(splus2R)
library(ifultools)
source("R/population-simulations/fractal_functions.R")


## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- 10

registerDoParallel(numCores)  # use multicore, set to the number of our cores

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

## run population models but only save time to extinction data, not entire series
run_100_sims_ER <- function(inc, start_colour, end_colour, Tmax, ICP) {
  ## run 100 simulations per colour:
  all <- foreach (z = 1:100, .combine=rbind)  %dopar% {
    
    noise = rep(NA, Tmax)
    
    ## create variable is_empty to keep track of whether time series is empty or not
    is_empty = length(which(is.na(noise))) == Tmax
    
    ## while time series is empty 
    while(is_empty) {
      
      ## generate increasing, decreasing, and stable coloured noise: 
      noise <- tvfd_single(start_colour = start_colour, 
                           end_colour = end_colour,
                           num_steps = 20,
                           Tmax = Tmax, 
                           inc = inc)
      
      ## check whether still empty
      is_empty = length(which(is.na(noise))) == Tmax 
    }
    
    ## force carrying capacity based on environment 
    K <- K0 + round(noise, digits = 0) 
    
    ## make sure none are negative
    K[which(K <= 0)] <- 5
    
    ## run pop models
    Nts = popmod(t = Tmax, N0 = N0, K = K) # C & Y model 
    
    ## calculate time to extinction
    t_extinction = first(which(Nts == 0))
    t_extinction = ifelse(is.na(t_extinction), Tmax, t_extinction)
    
    data.frame(t_extinction = t_extinction, 
               sim = z, 
               start_colour = start_colour,
               end_colour = end_colour,
               lambda = ICP)
  }
  return(all)
}

## write model to simulate and correct tvfd noise, no matter whether noise is increasing, decreasing, or stable:
tvfd_single <- function(start_colour, end_colour, num_steps, Tmax, inc) {
  
  by = (end_colour - start_colour)/num_steps
  each = ceiling(Tmax/num_steps)
  
  ## simulate the time series
  ## if no change in noise colour
  if(by == 0) {
    # make alpha constant
    alpha = rep(start_colour, Tmax)
  }
  ## other wise, create changing alpha
  else {
    alpha = seq(start_colour, end_colour, by = by)
    alpha = rep(alpha, each = each)[1:Tmax]
  }
  
  ## calculate delta
  delta = -alpha/-2
  
  ## set the innovations variance to unity
  innovation <- rep(1, Tmax)
  
  ## simulate a time-varying FD process
  noise <- FDSimulate(delta = delta, innovation = innovation)
  
  ## if time series is stable, do not correct for stepwise jumps
  if(by == 0) {
    new_noise = c(noise)
  }
  else {
    ## if time series crosses from alpha < 1 to alpha > 1
    if(any(alpha < 1) & !all(alpha < 1)) {
      
      ## define point at which alpha switches to greater/less than 1
      if(inc == TRUE) {
        point = first(which(alpha >= 1))
      }
      else {
        point = last(which(alpha >= 1)) + 1
      }
      
      ## correct data so that each time step begins at last time step value when alpha < 1:
      breaks = seq(1, 220000, by = each)
      new_noise = c(noise)
      
      i = 2
      while(i < length(breaks) - 1) {
        if(inc == TRUE & breaks[i] >= point) {
          to_add = new_noise[breaks[i] - 1] - new_noise[breaks[i]]
          
          new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
        } 
        else if (inc == FALSE & breaks[i] <= point)  {
          to_add = new_noise[breaks[i] - 1] - new_noise[breaks[i]]
          
          new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
        }
        i = i + 1
      }
    }
    else if (any(alpha > 1)) {
      
      ## fix it so that each time step begins at last time step value:
      breaks = seq(1, 220000, by = each)
      new_noise = c(noise)
      
      i = 2
      while(i < length(breaks)) {
        ## increasing
        if(inc == TRUE) {
          to_add = new_noise[breaks[i] - 1] - new_noise[breaks[i]]
          
          new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
        }
        else if(inc == FALSE){
          ## decreasing
          to_add = new_noise[breaks[i] - 1] - new_noise[breaks[i]]
          
          new_noise[breaks[i]:(breaks[i+1] - 1)] = new_noise[breaks[i]:(breaks[i+1] - 1)] + to_add
        }
        
        i = i + 1
      }
    }
  }
  
  ## remove mean and change variance to 4140:
  new_noise <- new_noise*1/sqrt(var(new_noise, na.rm = TRUE))*sqrt(4140)
  new_noise <- new_noise - mean(new_noise,  na.rm = TRUE)
  
  return(new_noise)
}



#################################################
###                  SIMULATIONS               ## 
#################################################
min_colour = 1.149833
max_colour = 1.997999
min_change = -0.1895386
max_change = 0.2035578
start_cols = seq(from = min_colour, to = max_colour, by = 0.1)
lambda = c(0.1, 0.5, 0.9)

incs = seq(from = 0, to =  max_change+0.01, by = 0.01)
decs = seq(from = 0-0.01, to =  min_change-0.01, by = -0.01)

K0 = 100
N0 = 100
Tmax = 200000

all_results <- list()

## loop through each amount of increase/decrease
i = 1
while (i <= length(incs) | i <= length(decs)) {
  
  ## loop through each start colour
  col = 1
  while(col <= length(start_cols)) {
    start_colour = start_cols[col]
    
    ## loop through each icp
    icp = 1
    while (icp <= length(lambda)) {
      print(paste("On inc/dec number ", i, " and icp ", icp, " and start colour ", start_colour, sep = ""))
      
      ## if there are still increases left
      if(i <= length(incs)) {
        end_colour = start_colour + incs[i]
        
        inc_sims = run_100_sims_ER(inc = TRUE, 
                                start_colour = start_colour,
                                end_colour = end_colour, 
                                Tmax = Tmax,
                                ICP = lambda[icp])
        
        inc_sims$spec_change = incs[i]
        
        ## write results of 100 sims for inc x icp combo
        write.csv(inc_sims, paste("pop-sims/marine_icp-", lambda[icp], "_change-", 
                                  incs[i], "_start-col-", start_colour, ".csv", sep = ""), row.names = F)
      }
      
      if(i <= length(decs)) {
        end_colour = start_colour + decs[i]
        
        dec_sims = run_100_sims_ER(inc = FALSE, 
                                   start_colour = start_colour,
                                   end_colour = end_colour, 
                                   Tmax = Tmax,
                                   ICP = lambda[icp])
        
        dec_sims$spec_change = decs[i]
        
        ## write results of 100 sims for dec x icp combo
        write.csv(dec_sims, paste("pop-sims/marine_icp-", lambda[icp], "_change-", 
                                  decs[i], "_start-col-", start_colour, ".csv", sep = ""), row.names = F)
      }
      
      ## move to next icp
      icp = icp + 1
    }
    ## move to next start colour 
    col = col + 1
  }
  
  ## move to next amount of change:
  i = i + 1
}

