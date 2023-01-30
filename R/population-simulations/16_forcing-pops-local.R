## simulating population dynamics under changes in autocorrelation
library(tidyverse)
library(ggplot2)
library(fractal)
library(foreach)
library(doParallel)

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
numCores

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
  alpha_stable = rep(c(start_colour, start_colour + 0.000001), num_steps)
  alpha_stable = rep(alpha_stable, each = each)[1:200000]
  
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
  
  else if (any(alpha >= 1)) {
   
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
  
  ## run 1000 simulations per model:
  icp = 1
  while (icp <= length(lambda)) {
    print(paste("On base colour ", col, " and icp ", icp, sep = ""))
    
    ## run 100 simulations per colour:
    all <- foreach (z = 1:1000, .combine=rbind)  %dopar% {
      
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
  
      ## save extinction risk metrics 
      data.frame(sim = z,
                 N_stable = first(which(Nts_stable == 0)),
                 N_inc = first(which(Nts_inc == 0)),
                 N_dec = first(which(Nts_dec == 0)),
                 lambda = lambda[icp],
                 base_col = start_colour)
      
      ## save whole time series 
      # data.frame(t = 1:Tmax, sim = z, 
      #            N_stable = Nts_stable,
      #            N_inc = Nts_inc,
      #            N_dec = Nts_dec,
      #            lambda = lambda[icp], 
      #            noise_stable = noise_stable, 
      #            noise_inc = noise_inc, 
      #            noise_dec = noise_dec, 
      #            base_col = start_colour)
    }
    
    ## write results of 1000 sims for base colour x icp combo
    write.csv(all, paste("data-processed/pop-sims_semi-stable/tvfdpopdynam_icp-", lambda[icp], "_base-col-", 
                         start_colour, "_20-steps_1000.csv", sep = ""), row.names = F)
    
    icp = icp + 1
  }
  ## move to next colour
  col = col + 1
}

## plot a few:
all %>%
  filter(sim %in% 1) %>%
  ggplot(., aes(x = t, y = N_inc)) + geom_line() + facet_wrap(~sim) +
  geom_line(aes(y = noise_stable + 100), colour = "red")



########################################################
###            calculate extinction risk              ## 
########################################################
## read in 100 simulations for each noise colour and calculate time to extinction 
## analyze change in noise colour over time windows 
base_col = seq(0, 3, by = 0.1)
lambda = c(0.1, 0.5, 0.9)

ER_all <- data.frame()

## loop through icps
icp = 1
while (icp <= length(lambda)) {
  
  ## loop through colours
  col = 1
  while (col < length(base_col)) {
    start_colour = base_col[col]
    
    filename = paste("data-processed/tvfd-pop-sims/tvfdpopdynam_icp-", lambda[icp], "_base-col-", 
                     start_colour, "_20-steps.csv", sep = "")
    
    if(file.exists(filename)) {
      
      ## read in data:
      sims <- read.csv(filename)
      
      ## calculate time to extinction for each time series 
      ER <- sims %>%
        group_by(sim) %>%
        summarise(N_stable = first(which(N_stable == 0)),
                  N_inc = first(which(N_inc == 0)),
                  N_dec = first(which(N_dec == 0))) %>%
        mutate(lambda = sims$lambda[1],
               start_colour = start_colour)
      
      ## bind to rest of the data:
      ER_all <- rbind(ER_all, ER)
      
      print(paste("On icp: ", lambda[icp], " and colour: ", base_col[col],
                  sep = ""))
    }
    
    col = col + 1
  }
  icp = icp + 1
}

#write.csv(ER_all, "data-processed/tvfd-pop-sims/ER-metrics_20-steps.csv", row.names = F)
ER <- read.csv("data-processed/tvfd-pop-sims/ER-metrics_20-steps.csv")

length(which(!is.na(ER$N_stable)))
length(which(!is.na(ER$N_inc)))
length(which(!is.na(ER$N_dec)))

## reformat data for plotting 
ER <- ER %>%
  gather(key = "noise_type", value = "time_to_extinction", c(N_stable, N_inc, N_dec)) 

## replace NAs (cases where population doesn't go extinct) with 200000
ER <- ER %>%
  mutate(time_to_extinction = ifelse(is.na(time_to_extinction), 
                                     200000,
                                     time_to_extinction)) 

## make start noise colour column accurate for decreasing noise
ER <- ER %>%
  mutate(start_colour = ifelse(noise_type == "N_dec", 
                                     start_colour + 0.1,
                                     start_colour)) 

## calculate mean time to extinction, coefficient of variation in time to extinction 
ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction))%>%
  unique(.) %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, colour = noise_type)) + geom_point() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "Mean persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) 

ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(cv = sd(time_to_extinction)/mean(time_to_extinction)) %>%
  unique(.) %>%
  ggplot(., aes(y = cv, x = start_colour, colour = noise_type)) + geom_point() +
  theme_light() +
  labs(x = "Starting noise colour", y = "CV persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda)

## try a new plot:
## arrow from N_stable to point (N_inc, start_colour + 0.1 ) to show how stable changes under 
pts = ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction))%>%
  filter(noise_type %in% c("N_stable"))

dec = ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction))%>%
  filter(noise_type %in% c("N_stable", "N_dec")) 

ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction))%>%
  filter(noise_type %in% c("N_stable", "N_inc")) %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, group = start_colour)) + 
  geom_line(colour = "red", aes(x = start_colour - 0.02)) +
  geom_line(data = dec, colour = "blue", aes(x = start_colour + 0.02)) +
  geom_point(data = pts) +
  scale_y_log10() +
  theme_light() +
  labs(x = "Noise colour", y = "Mean persistence time") +
  facet_wrap(~lambda) 


### next:
## plot some 
## se what happens when the num_steps parameter is changed - same results?


## try making triple axis plot for terrestrial environment, marine environment 

### TERRESTRIAL ###
# min_colour = 0.07916834
# max_colour = 1.507507 
# min_change = -0.07106262
# max_change = 0.05052732
# start_cols = seq(from = min_colour, to = max_colour, by = 0.1)
# lambda = c(0.1, 0.5, 0.9)
# 
# incs = seq(from = 0, to =  max_change+0.005, by = 0.005)
# decs = seq(from = 0-0.005, to =  min_change-0.005, by = -0.005)

# ### MARINE ###
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
        write.csv(inc_sims, paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[icp], "_change-", 
                                  incs[i], "_start-col-", start_colour, "_corrected.csv", sep = ""), row.names = F)
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
        write.csv(dec_sims, paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[icp], "_change-", 
                                  decs[i], "_start-col-", start_colour, "_corrected.csv", sep = ""), row.names = F)
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

## write function to simulate 100 time series with the same change in noise colour, length, icp
run_100_sims <- function(inc, start_colour, end_colour, Tmax, ICP) {
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
    
    data.frame(t = 1:Tmax, 
               sim = z, 
               N = Nts,
               start_colour = start_colour,
               end_colour = end_colour,
               lambda = ICP, 
               noise = noise)
  }
  return(all)
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
  ## if time series crosses from alpha < 1 to alpha > 1
  else if(any(alpha < 1) & !all(alpha < 1)) {
    
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
  ## remove mean and change variance to 4140:
  new_noise <- new_noise*1/sqrt(var(new_noise, na.rm = TRUE))*sqrt(4140)
  new_noise <- new_noise - mean(new_noise,  na.rm = TRUE)
  
  return(new_noise)
}



## check data
## read in and combine
filenames = paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[1:3], "_change-", 
                  decs[1], "_start-col-", start_cols[1:length(start_cols)],  "_corrected.csv",sep = "") %>%
  append(., paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[1:3], "_change-", 
                  incs[1], "_start-col-", start_cols[1:length(start_cols)],  "_corrected.csv", sep = "")) %>%
  append(., paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[1:3], "_change-", 
                incs[2], "_start-col-", start_cols[1:length(start_cols)],  "_corrected.csv", sep = ""))%>%
  append(., paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[1:3], "_change-", 
                  decs[2], "_start-col-", start_cols[1:length(start_cols)],  "_corrected.csv", sep = ""))%>%
  append(., paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[1:3], "_change-", 
                  incs[3], "_start-col-", start_cols[1:length(start_cols)], "_corrected.csv", sep = "")) %>%
  append(., paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[1:3], "_change-", 
                  decs[3], "_start-col-", start_cols[1:length(start_cols)],  "_corrected.csv",sep = ""))%>%
  append(., paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[1:3], "_change-", 
                  incs[4], "_start-col-", start_cols[1:length(start_cols)],  "_corrected.csv",sep = "")) %>%
  append(., paste("data-processed/tvfd-pop-sims/marine_icp-", lambda[1:3], "_change-", 
                  decs[4], "_start-col-", start_cols[1:length(start_cols)],  "_corrected.csv", sep = ""))

data = c()
for (i in 1:length(filenames)) {
  if(file.exists(filenames[i])) {
    data = rbind(data, read.csv(filenames[i]))
  }
}

data %>%
  group_by(spec_change, lambda, start_colour) %>%
  mutate(mean_pers_time = mean(t_extinction))%>%
  unique(.) %>%
  mutate(inc_or_dec = ifelse(spec_change < 0, "decrease", 
                             ifelse(spec_change > 0, "increase", 
                                    ifelse(spec_change == 0, 
                                    "stable", 
                                    NA)))) %>%
  filter(spec_change %in% c(0.03, 0.01, 0, -0.01, -0.04)) %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, colour = inc_or_dec,
                shape = as.factor(spec_change))) + 
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "Mean persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) 
  


ER$source = "general"
ER$spec_change = 0.1
data$source = "marine"

ER <- ER %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction))%>%
  unique(.) %>%
  mutate(noise_type = ifelse(noise_type == "N_inc", 
                             "increase", 
                             ifelse(noise_type == "N_dec", 
                                    "decrease",
                                    "stable")))

data %>%
  group_by(spec_change, lambda, start_colour) %>%
  mutate(mean_pers_time = mean(t_extinction)) %>%
  unique(.) %>%
  mutate(noise_type = ifelse(spec_change < 0, "decrease", 
                             ifelse(spec_change > 0, "increase", 
                                    ifelse(spec_change == 0, 
                                           "stable", 
                                           NA)))) %>%
  filter(spec_change %in% c(0.03, 0.01, 0, -0.01, -0.04)) %>%
  select(-end_colour) %>%
  rbind(., ER) %>%
  filter(start_colour > 1, start_colour < 2) %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, colour = noise_type,
                shape = as.factor(spec_change))) + 
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "Mean persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) 

data %>%
  group_by(spec_change, lambda, start_colour) %>%
  mutate(mean_pers_time = mean(t_extinction)) %>%
  unique(.) %>%
  mutate(noise_type = ifelse(spec_change < 0, "decrease", 
                             ifelse(spec_change > 0, "increase", 
                                    ifelse(spec_change == 0, 
                                           "stable", 
                                           NA)))) %>%
  filter(spec_change %in% c( 0, 0.01, 0.02, 0.03)) %>%
  select(-end_colour) %>%
  rbind(., ER) %>%
  filter(noise_type != "decrease") %>%
  filter(start_colour > 1, start_colour < 2) %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, colour = noise_type,
                shape = as.factor(spec_change))) + 
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "Mean persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) 

data %>%
  group_by(spec_change, lambda, start_colour) %>%
  mutate(mean_pers_time = mean(t_extinction)) %>%
  unique(.) %>%
  mutate(noise_type = ifelse(spec_change < 0, "decrease", 
                             ifelse(spec_change > 0, "increase", 
                                    ifelse(spec_change == 0, 
                                           "stable", 
                                           NA)))) %>%
  filter(spec_change %in% c( 0, -0.01, -0.02,-0.03, -0.04)) %>%
  select(-end_colour) %>%
  rbind(., ER) %>%
  filter(noise_type != "increase") %>%
#  filter(start_colour > 1, start_colour < 2) %>%
  ggplot(., aes(y = mean_pers_time, x = start_colour, colour = noise_type,
                shape = as.factor(spec_change))) + 
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "Mean persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) 

ER <- ER %>%
  mutate(cv = sd(time_to_extinction)/mean(time_to_extinction)) %>%
  group_by(start_colour, noise_type, lambda) %>%
  mutate(mean_pers_time = mean(time_to_extinction))%>%
  unique(.) 
  
data %>%
  group_by(spec_change, lambda, start_colour) %>%
  mutate(mean_pers_time = mean(t_extinction)) %>%
  mutate(cv = sd(t_extinction)/mean(t_extinction)) %>%
  unique(.) %>%
  mutate(noise_type = ifelse(spec_change < 0, "decrease", 
                             ifelse(spec_change > 0, "increase", 
                                    ifelse(spec_change == 0, 
                                           "stable", 
                                           NA)))) %>%
  filter(spec_change %in% c(0.03, 0.01, 0, -0.01, -0.04)) %>%
  select(-end_colour) %>%
  rbind(., ER) %>%
  filter(start_colour > 1, start_colour < 2) %>%
  ggplot(., aes(y = cv, x = start_colour, colour = noise_type,
                shape = as.factor(spec_change))) + 
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "CV persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) 

data %>%
  group_by(spec_change, lambda, start_colour) %>%
  mutate(mean_pers_time = mean(t_extinction)) %>%
  mutate(cv = sd(t_extinction)/mean(t_extinction)) %>%
  unique(.) %>%
  mutate(noise_type = ifelse(spec_change < 0, "decrease", 
                             ifelse(spec_change > 0, "increase", 
                                    ifelse(spec_change == 0, 
                                           "stable", 
                                           NA)))) %>%
  filter(spec_change %in% c( 0, 0.01, 0.02, 0.03)) %>%
  select(-end_colour) %>%
  rbind(., ER) %>%
  filter(noise_type != "decrease") %>%
  filter(start_colour > 1, start_colour < 2) %>%
  ggplot(., aes(y = cv, x = start_colour, colour = noise_type,
                shape = as.factor(spec_change))) + 
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "CV persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) 

data %>%
  group_by(spec_change, lambda, start_colour) %>%
  mutate(mean_pers_time = mean(t_extinction)) %>%
  mutate(cv = sd(t_extinction)/mean(t_extinction)) %>%
  unique(.) %>%
  mutate(noise_type = ifelse(spec_change < 0, "decrease", 
                             ifelse(spec_change > 0, "increase", 
                                    ifelse(spec_change == 0, 
                                           "stable", 
                                           NA)))) %>%
  filter(spec_change %in% c( 0, -0.01, -0.02,-0.03, -0.04)) %>%
  select(-end_colour) %>%
  rbind(., ER) %>%
  filter(noise_type != "increase") %>%
  filter(start_colour > 1, start_colour < 2) %>%
  ggplot(., aes(y = cv, x = start_colour, colour = noise_type,
                shape = as.factor(spec_change))) + 
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Starting noise colour", y = "CV persistence time", 
       colour = "Change in autocorrelation:") +
  facet_wrap(~lambda) 

