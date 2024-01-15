## creating simple metapopulation model 
library(tidyverse)
library(foreach)
library(doParallel)

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- 10

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
    
    ## remove mean and change variance to 0.625^2:
    noise <- noise*1/sqrt(var(noise))*sqrt(0.625^2)
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

L = 1000 # number of time steps
reps = 100 # number of replicates; they used 500
n_pops = 10 # number of subpopulaions
Km = 100 # mean carrying capacity 
lambda = c(0.1, 1) # intraspecifc competition parameter 
r = c(0.4, 0.8, 1.2) # growth rate values 
beta = 0.5

## dispersal parameters:
d = 0.1 # fixed proportion of individuals that emigate the patch

## set lambda and r
icp = 1
r_n = 2


## loop through different amounts of synchronization 
p = c(0, 0.5, 1)
betas = seq(0, 2, 0.1)
all_reps_pops <- c()
all_reps_noise <- c()
p_curr = 1
while(p_curr <= length(p)) {
  
  ## loop through different amounts of autocorrelation
  beta_curr <- 1
  while(beta_curr <= length(betas)) {
    
    ## loop through replicates 
    all <- foreach (rep = 1:100, .combine=rbind)  %dopar% {
      ##### NOISE #####
      ## generate the synchronous noise
      noise <- one_over_f_synchro(beta = betas[beta_curr], p = p[p_curr], n_ts = n_pops, L = L)
      
      ## first-difference the time series values
      ts_all <- data.frame(sapply(noise[[1]], cbind))
      diffs <- ts_all[2:(nrow(ts_all)-1),] - ts_all[1:(nrow(ts_all)-2),]
      
      ## measure cross correlation:
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
      
      ## save the noise 
      ts_all$time = rep(1:nrow(ts_all))
      colnames(ts_all)[1:(ncol(ts_all)-1)] <- paste0("tsnum_", 1:(ncol(ts_all)-1))
      ts_all$replicate = rep
      ts_all$p = p[p_curr]
      ts_all$p_star = p_star
      ts_all$beta = betas[beta_curr]
      ts_all$beta_star = mean(sapply(noise[[2]], cbind))
      
      all_reps_noise <- rbind(all_reps_noise, ts_all)
      
      ##### POPULATION DYNAMICS #####
      ## let noise affect carrying capacity 
      K = list()
      for(i in 1:length(noise[[1]])) {
        noiselist = noise[[1]]
        K_temp = Km*(1+noiselist[[i]])
        K_temp[K_temp < 1] = 0
        K[[i]] = K_temp
        ## set sizes below 1 to 0
      }
      
      ## set initial population size of subpopulations 
      N_j = as.list(rep(round(Km/2), n_pops)) # initial population size of subpopulations, set to Km/2
      ## for each time step t:
      N_global <- sum(sapply(N_j, sum))
      i = 1
      while(i < L){
        
        curr_allpops_new <- sapply(N_j, function(x) x[i])
      
        ##### reproduce #####
        # for each subpopulation, calculate the new population size after reproduction:
        new_sizes <- c()
        y=1
        while(y <= n_pops) {
          N = curr_allpops_new[y]
          K_curr = K[[y]]
          
          new_sizes = append(new_sizes, 
                             round(as.numeric(N*exp(r[r_n]*(1 - as.complex(N/(K_curr[i]))^lambda[icp])))))
          
          y = y+1
        }
        ## if pop size < 0, set to 0
        new_sizes[which(new_sizes < 0)] <- 0
        new_sizes[which(is.na(new_sizes))] <- 0
        
        N_j <- mapply(c, N_j, new_sizes, SIMPLIFY=FALSE)
        
        ##### census ####
        ## calculated global population size (sum of all pop sizes)
        N_global <- append(N_global, sum(sapply(new_sizes, sum)))
        
        i = i + 1
      }
      
      ## save run:
      pops_all <- data.frame(sapply(N_j, cbind))
      pops_all$time = rep(1:nrow(pops_all))
      colnames(pops_all)[1:(ncol(pops_all)-1)] <- paste0("popnum_", 1:(ncol(pops_all)-1))
      pops_all$replicate = rep
      pops_all$p = p[p_curr]
      pops_all$p_star = p_star
      pops_all$beta = betas[beta_curr]
      pops_all$beta_star = mean(sapply(noise[[2]], cbind))
      pops_all$N_global = N_global
      
      pops_all
    }
    
    all_reps_pops <- rbind(all_reps_pops, all)
    
    print(paste0("Finished p=", p[p_curr], " and beta=", betas[beta_curr]))
    
    beta_curr = beta_curr + 1
  }
   
  p_curr = p_curr + 1
}

write.csv(all_reps_pops, "data-processed/metapop-models/simple-model_no-dispersal_r0.8_lambda0.1_popdynam.csv",
          row.names = FALSE)
write.csv(all_reps_noise, "data-processed/metapop-models/simple-model_no-dispersal_r0.8_lambda0.1_noise.csv",
          row.names = FALSE)

all_reps_pops <- read.csv("data-processed/metapop-models/simple-model_no-dispersal_r0.8_lambda0.1_popdynam.csv")
  

## calculate extinction risk:
## proportion of extinct populations / all replicates 
ext_summary <- all_reps_pops %>%
  gather(., key = "popnum", value = "size", 
         paste0("popnum_", 1:10)) %>%
  group_by(replicate, p, beta) %>%
  mutate(goes_extinct = any(N_global == 0)) %>%
  select(replicate, p, beta, goes_extinct) %>%
  unique() %>%
  ungroup() %>%
  group_by(p, beta) %>%
  mutate(extinction_risk = length(which(goes_extinct == "TRUE")) /reps) %>%
  ungroup() %>%
  select(beta, p, extinction_risk) %>%
  unique()

## there should be 1 for each beta x p 

## plot noise colour vs. extinction risk 
ext_summary %>%
  ggplot(aes(x = beta, y = extinction_risk, colour = p, group = p)) +
  geom_line() + 
  theme_bw()

## what about local extinction risk?
loc_extinction <- all_reps_pops %>%
  gather(., key = "popnum", value = "size", 
         paste0("popnum_", 1:10)) %>%
  group_by(p, beta) %>%
  mutate(num_local_extinctions = length(which(size == 0))) %>%
  select(p, beta, num_local_extinctions) %>%
  unique() %>%
  ungroup() 

loc_extinction %>%
  ggplot(aes(x = beta, y = num_local_extinctions, colour = p, group = p)) +
  geom_line() + 
  theme_bw()

## next steps: 
## check these results 
## run for different noise colours across same icp x r 
## (compare to figures of paper)





##### GARBAGE ######

split_pops <- group_split(all_reps_pops, replicate)


pops <- bind_rows(split_pops[1:10])


## plot the pops:
pops_all <- gather(pops, key = "popnum", value = "size", 
                 paste0("popnum_", 1:10))

pops_all %>%
  ggplot(aes(x = time, y = size, colour = popnum, group = popnum)) + 
  geom_path(linewidth = 0.5) +
  theme_bw() +
  theme(legend.position = "none") +
  facet_wrap(~replicate)

min(all_reps_pops$N_global)


## plot the noise
noise <- one_over_f_synchro(beta = beta, p = p, n_ts = n_pops, L = L)

## first-difference the time series values
ts_all <- data.frame(sapply(noise[[1]], cbind))
ts_all$time = rep(1:nrow(ts_all))
colnames(ts_all)[1:10] <- paste0("tsnum_", 1:10)

ts_all <- gather(ts_all, key = "tsnum", value = "noise",
                 paste0("tsnum_", 1:10))

ts_all %>%
  ggplot(aes(x = time, y = noise, colour = tsnum, group = tsnum)) +
  geom_path(linewidth = 0.5) +
  theme_bw() +
  theme(legend.position = "none")


one = filter(ts_all, tsnum == "tsnum_1")
sd(one$noise)
mean(one$noise)
