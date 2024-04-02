## building spatially explicit coupled map lattice 
#.libPaths( c( "~/projects/def-jsunday/nikkim/packages" , .libPaths() ) )
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)

## set path
path = "data-processed/metapop-models/climate-shift/"

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
numCores

registerDoParallel(numCores)  # use multicore, set to the number of our cores

## function to create n_ts series of 1/f noise that have p degree of synchrony
## method adapted from Lögdberg & Wennergren
## when p = 1, all are synchronous 
## when p = 0, all are asynchronous
one_over_f_synchro <- function(beta, p, n_ts, L){
  ## set shared random element of phase:
  randf <- runif(524288, 0, 2*pi) 
  
  ## for each requested time series:
  list_ts <- list()
  list_colours <- list()
  x = 1
  while(x <= n_ts) {
    ## create a series of 1/fB noise as described in Cuddington and Yodzis
    n = rnorm(524288, mean = 0, sd = 0.01) ## random numbers with 0 mean and 0.01 sd 
    f = 1:524288 ## f
    a <- 1/(f^(beta/2))*exp(1)^n ### amplitudes = 1/f^beta/2 * tiny random component
    
    ## apply phase shift 
    phases <- randf + (1-p)*runif(1, 0, 2*pi)
    
    ## calculate wave coeffs and inverse dft
    complex <- a*cos(phases) + a*sin(phases) ## complex coefficients
    
    dft <- fft(complex, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
    noise = as.numeric(dft[1:L]) ## crop the noise series to first L points
    
    ## remove mean and change variance to 1:
    noise <- noise*1/sqrt(var(noise))*sqrt(1^2)
    noise <- noise - mean(noise)
    
    ## estimate noise colour from a linear regression of power spectrum:
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
    
    ## save 
    list_ts[[x]] <- noise
    list_colours[[x]] <- as.numeric(true_colour$coefficients[2])
    
    x = x + 1
  }
  
  return(list(list_ts, list_colours))
}

r = 1.2 ## set maximum growth rate 
K = 100 ## set mean carrying capacity
d = 0.1 ## set prop of offspring dispersing 
L = 100 ## set length of time series 
nrow = 100
ncol = 100

## work in latitudinal variation in r
h = 10 ## half-saturation constant that defines distance at which Eit = 0.5
s = -3 ## shape parameter defining direction (negative = negative slope) and shape (s > 1 gives signmoid)
Emax = 1 ## set Emax to 1
n_pops = nrow*ncol # number of grid cells
lambda = c(0.1, 1) # intraspecifc competition parameter 
icp = 1 ## set lambda
reps = 10 ## set number of replicates 
p = c(0, 0.75, 1)

## read the command line arguments 
command_args <- commandArgs(trailingOnly = TRUE)
i = as.numeric(command_args[1])
p_curr = i

## loop through noise colours
betas = seq(0, 2, 0.5)
beta_curr = 1
all_reps_pops <- c()
while(beta_curr <= length(betas)) {
  
  print("Script started!")
  
  ## loop through replicates 
  all <- foreach (rep = 1:reps, .combine=rbind)  %dopar% {
    
    ## cellular lattice of 100 x 600 habitat patches
    lattice_N_it = array(dim = c(nrow,ncol,L))
    lattice_r = matrix(ncol = ncol, nrow = nrow)
    lattice_E_it = matrix(ncol = ncol, nrow = nrow)
    
    ## higher proportion dispersing = less pronounced effect of suitability gradient
    lattice_r[1:nrow,1:ncol] <- r ## start with growth rate = max growth rate 
    lattice_N_it[1:nrow,1:ncol,1] <- K/2 ## start with population size = carrying capacity / 2
    
    ## position optimum climatic conditions as a row 50 on the lattice (Emax)
    lattice_E_it[100/2,] = Emax
    
    ## assume that conditions decline sigmoidally away from this optimum in both directions
    ## for nrow/10 cells
    lattice_E_it[1:(100/2-1),] = Emax*rev((1:(100/2-1))^s/ ((1:(100/2-1))^s + h^s))
    lattice_E_it[(100/2+1):nrow,] = Emax*((1:nrow)^s / ((1:nrow)^s + h^s))[1:(nrow-100/2)]
    # plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
    
    ## replicate latitudinal gradient 1000 times 
    lattice_E_it_array <- replicate(L, lattice_E_it)
    
    ## create autocorrelated time series for each cell, 2000 time steps long
    lattice_ac_it = array(dim = c(nrow,ncol,L))
    
    ## start with complete synchrony, high autocorrelation
    noise <- one_over_f_synchro(beta = betas[beta_curr],
                                p = p[p_curr],
                                n_ts = nrow*ncol, ## number of cells
                                L = L)
    
    ## calculate mean measured beta
    beta_star =  mean(sapply(noise[[2]], cbind))
    
    ## measure cross correlation:
    ## first-difference the time series values
    ts_all <- data.frame(sapply(noise[[1]], cbind))
    diffs <- ts_all[2:(nrow(ts_all)-1),] - ts_all[1:(nrow(ts_all)-2),]
    
    cors <- c()
    i=1
    while(i <= ncol(diffs)) {
      for(n in ncol(diffs)) {
        if(n != i) {
          cors <- append(cors, cor(diffs[,i], diffs[,n], method = 'pearson'))
        }
      }
      i = i + 1
    }
    p_star = mean(cors)
    
    
    ## assign noise to each cell in the lattice
    df <- data.frame(sapply(noise[[1]], cbind))
    n = 1
    for(a in 1:nrow) {
      for(b in 1:ncol) {
        lattice_ac_it[a,b,] <- df[,n]
        n = n + 1
      }
    }
    
    ## set carring capacity to max carrying capacity for all cells
    lattice_K_it <- K*lattice_ac_it ## multiplicative
    lattice_K_it[,,] <- K
    
    ## prevent extinction from catastrophe (only allow extinction from demographic stochasticity)
    lattice_K_it[lattice_K_it < 1] = 0.1 ## was 0
    #plot(lattice_K_it[1,1,])
    
    ## let noise affect r for each cell in the lattice
    lattice_r_array = (lattice_E_it_array + lattice_ac_it)*r
    # plot(x = 1:nrow, y = lattice_r_array[1:nrow,1,1])
    
    ## run at stable climate conditions for 500 time steps
    t = 1 
    N_global <- c()
    N_ext_local <- 0
    while(t < L) {
      
      ## get grid of growth rates at time t
      lattice_r_curr <- lattice_r_array[,,t]
      
      ## get grid of carrying capacity at time t
      lattice_K_curr <- lattice_K_it[,,t]
      
      ## get grid of current population size at time t
      lattice_N_curr <- lattice_N_it[,,t]
      
      ##### DISPERSE #####
      ## fixed proportion of individuals have same probability of dispersing into any patch
      curr_allpops <- lattice_N_curr
      curr_allpops_new <- curr_allpops
      ## loop through grid cells
      x = 1
      while (x <= ncol) {
        y = 1
        while(y <= nrow) {
          ## if all pops are dead, stop
          if(length(which(curr_allpops == 0)) == n_pops) {
            y = n_pops + 1
            i = L
          }
          else {
            ## get number of dispersers 
            dispersers = curr_allpops[y,x]*d
            ## figure out where each one disperses 
            ## individuals can disperse into neighbouring 8 cells
            othercells <- expand.grid(x = c(x-1,x,x+1), y = c(y-1,y,y+1)) 
            othercells = othercells[!(x == othercells$x & y == othercells$y),]
            othercells$num = 1:8
            
            ## randomly sample from 1:8
            sample <- data.frame(num = sample(1:8, dispersers, replace = TRUE))
            
            othercells <- left_join(sample, othercells, by = "num") %>%
              filter(!(y == 0 | x == 0 | y > nrow | x > ncol))
            
            ## subtract dispersers
            curr_allpops_new[y,x] <- curr_allpops_new[y,x] - dispersers
            
            if(length(othercells) != 0) {
              ## subtract dispersers from population size, add them to population size at their new home   
              v = 1
              while(v <= nrow(othercells)) {
                curr_allpops_new[othercells$y[v], othercells$x[v]] = 
                  curr_allpops_new[othercells$y[v], othercells$x[v]] + 1
                v = v + 1
              }
            }
            y = y + 1
          }
        }
        x = x + 1
      }
      
      
      ##### REPRODUCE #####
      # for each subpopulation, calculate the new population size after reproduction:
      new_sizes <- matrix(ncol = ncol, nrow = nrow)
      K_curr = lattice_K_it[,,t]
      r_curr = lattice_r_array[,,t]
      x = 1
      while (x <= ncol) {
        y = 1
        while(y <= nrow) {
          N = curr_allpops_new[y,x]
          
          new_sizes[y,x] = round(as.numeric(N*exp(r_curr[y,x]*(1 - as.complex(N/(K_curr[y,x]))^lambda[icp]))))
        
          ## add demographic stochasticity
          new_sizes[y,x] = sample(rpois(new_sizes[y,x], n = 1000), size = 1)
          
          y = y + 1
        } 
        x = x + 1
      }
      ## if pop size < 0, set to 0
      new_sizes[which(new_sizes < 0)] <- 0
      new_sizes[which(is.na(new_sizes))] <- 0
      
      ## update population matrix:
      lattice_N_it[,,t+1] = new_sizes
      
      ## save stats:
      ## calculate global population size (sum of all pop sizes)
      N_global <- append(N_global, sum(lattice_N_it[,,t+1]))
      ## calculate local population extinction frequency (how many pops go extinct per time step)
      N_ext_local = N_ext_local + length(which(new_sizes == 0))
    
      ## on last step, save grid of population sizes
      if((t+1) == L) {
        # ras <- raster(nrow = 100, ncols = 100, xmn = 0, xmx = 100, ymn = 0, ymx = 100)
        # values(ras) <- c(t(lattice_N_it[,,t]))
        # plot(ras)
        # values(ras) <- c(t(lattice_r_array[,,t]))
        # plot(ras)
       
        saveRDS(lattice_N_it[,,t+1], paste0(path, "lattice-N-it_beta_",
                                          betas[beta_curr], "_p", p[p_curr], "_icp_", lambda[icp], 
                                          "_population_", rep, "_t100.rds"))
       
      }
      
      t = t + 1
      
    }
    
    
    ## save run:
    pops_all <- data.frame(loc_Emax = rep((nrow/10)/2,L-1),
                           time = rep(1:length(N_global)),
                           replicate = rep,
                           p = p[p_curr],
                           p_star = p_star,
                           beta = betas[beta_curr],
                           beta_star = beta_star,
                           N_global = N_global,
                           N_ext_local = rep(N_ext_local,L-1)) 
    pops_all
   }
  
  ## save
  write.csv(all, paste0(path, "all_beta", betas[beta_curr], 
                        "_p", p[p_curr], "_icp_", lambda[icp], "_r_t100.csv"), 
            row.names = FALSE)

  print(paste0("Finished p=", p[p_curr], " and beta=", betas[beta_curr]))

  beta_curr = beta_curr + 1
}

