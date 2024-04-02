## building spatially explicit coupled map lattice 
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
theme_set(theme_bw())

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

r = 1.2 ## set maximum growth rate 
K = 100 ## set mean carrying capacity
d = 0.1 ## set prop of offspring dispersing 
L = 1000 ## set length of time series 
nrow = 30
ncol = 10

## work in latitudinal variation in r
h = 10 ## half-saturation constant that defines distance at which Eit = 0.5
s = -3 ## shape parameter defining direction (negative = negative slope) and shape (s > 1 gives signmoid)
Emax = 1 ## set Emax to 1
n_pops = nrow*ncol # number of grid cells
lambda = c(0.1, 1) # intraspecifc competition parameter 
icp = 1 ## set lambda
reps = 100 ## set number of replicates 
p = c(0, 0.75, 1)
p_curr = 1


## loop through noise colours
betas = seq(0, 2, 0.1)
beta_curr = 1
all_reps_pops <- c()
while(beta_curr <= length(betas)) {
  
  ## loop through replicates 
  all <- foreach (rep = 1:reps, .combine=rbind)  %dopar% {
    
    ## cellular lattice of 100 x 300 (later 2700) habitat patches
    lattice_N_it = array(dim = c(nrow,ncol,L))
    lattice_r = matrix(ncol = ncol, nrow = nrow)
    lattice_E_it = matrix(ncol = ncol, nrow = nrow)
    
    ## higher proportion dispersing = less pronounced effect of suitability gradient
    lattice_r[1:nrow,1:ncol] <- r ## start with growth rate = max growth rate 
    lattice_N_it[1:nrow,1:ncol,1] <- K/2 ## start with population size = carrying capacity / 2
    
    ## position optimum climatic conditions as a row on the lattice (Emax)
    lattice_E_it[30/2,] = Emax
    
    ## assume that conditions decline sigmoidally away from this optimum in both directions
    ## for nrow/10 cells
    lattice_E_it[1:(30/2-1),] = Emax*rev((1:(30/2-1))^s/ ((1:(30/2-1))^s + h^s))
    lattice_E_it[(30/2+1):nrow,] = Emax*((1:nrow)^s / ((1:nrow)^s + h^s))[1:(nrow-30/2)]
    # plot(x = 1:nrow, y = lattice_E_it[1:nrow,1])
    
    ## replicate latitudinal gradient 1000 times 
    lattice_E_it_array <- replicate(L, lattice_E_it)
    ## set growth rate r to equal max growth rate*environmental suitability
    lattice_r_array = lattice_E_it_array*r
    # plot(x = 1:nrow, y = lattice_r_array[1:nrow,1,1])
    
    ## create autocorrelated time series for each cell, 1000 time steps long
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
    
    ## let noise affect K for each cell in the lattice
    ## create latitudinal gradient in mean K
    lattice_K_it <- K*lattice_E_it[,1] + 100*lattice_ac_it ## multiplicative
    #plot(lattice_K_it[1,1,])
    
    ## older versions:
    #lattice_K_it <- K*lattice_ac_it ## multiplicative
    #lattice_K_it <- K + lattice_ac_it ## additive 
    
    ## prevent extinction from catastrophe (only allow extinction from demographic stochasticity)
    lattice_K_it[lattice_K_it < 1] = 0.1 ## was 0
    #plot(lattice_K_it[1,1,])
    
    ## run at stable climate conditions for 1000 time steps
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
            ## start with global dispersal:
            ## randomly sample x and y pairs:
            othercells <- expand.grid(1:nrow, 1:ncol)
            othercells <- filter(othercells, !(Var1 == y & Var2 == x))
            
            ## subtract dispersers
            curr_allpops_new[y,x] <- curr_allpops_new[y,x] - dispersers
            
            if(length(othercells) != 0) {
              ## if some pops can still be dispersed to
              new_homes = sample(1:nrow(othercells), dispersers, replace = TRUE)
              ## subtract dispersers from population size, add them to population size at their new home   
              for(v in new_homes) {
                curr_allpops_new[othercells$Var1[v], othercells$Var2[v]] = 
                  curr_allpops_new[othercells$Var1[v], othercells$Var2[v]] + 1
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
    
      t = t + 1
      
      ## on last step, save grid of population sizes
      if(t == 1000) {
        # ras <- raster(nrow = 300, ncols = 50, xmn = 0, xmx = 10, ymn = 0, ymx = 300)
        # values(ras) <- c(t(lattice_N_it[,,t]))
        # plot(ras)
       
        saveRDS(lattice_N_it[,,t], paste0("data-processed/metapop-models/grid_beta", betas[beta_curr], 
                                          "_p", p[p_curr], "_icp_", lambda[icp], 
                                          "_population_", rep, ".rds"))
       
      }
    }
    
    ## save run:
    pops_all <- data.frame(loc_Emax = rep((nrow/10)/2,999),
                           time = rep(1:length(N_global)),
                           replicate = rep,
                           p = p[p_curr],
                           p_star = p_star,
                           beta = betas[beta_curr],
                           beta_star = beta_star,
                           N_global = N_global,
                           N_ext_local = rep(N_ext_local, 999)) 
    pops_all
   }
  
  ## save
  write.csv(all, paste0("data-processed/metapop-models/all_beta", betas[beta_curr], 
                        "_p", p[p_curr], "_icp_", lambda[icp], "_no-catastrophe_5_additive_gradK.csv"), 
            row.names = FALSE)

  print(paste0("Finished p=", p[p_curr], " and beta=", betas[beta_curr]))

  beta_curr = beta_curr + 1
}

## results for complete synchrony
all <- read.csv("data-processed/metapop-models/all_beta0_p0.csv")

all %>%
  mutate(unique_id = paste(replicate, beta)) %>%
  ggplot(aes(x = time, y = N_global, colour = beta, group = unique_id)) + geom_path()



## calculate extinction risk for each beta:
ext_summary <- c()
l=1
while(l <= length(betas)) {
  p_curr = 1
  while(p_curr <= length(p)) {
    filename = paste0("data-processed/metapop-models/all_beta", betas[l],
                      "_p", p[p_curr], "_icp_", lambda[icp], "_no-catastrophe_5_additive.csv")
    # filename = paste0("data-processed/metapop-models/all_beta", betas[l], "_p", p[p_curr],
    #                  ".csv")
    if(file.exists(filename)) {
      all <- read.csv(filename)
      
      ext_summary <- all %>%
        group_by(p, beta) %>%
        mutate(mean_local_freq_ext = mean(N_ext_local)) %>%
        ungroup() %>%
        group_by(replicate, p, beta) %>%
        mutate(goes_extinct = any(N_global == 0)) %>%
        select(replicate, p, beta, goes_extinct, mean_local_freq_ext) %>%
        unique() %>%
        ungroup() %>%
        group_by(p, beta) %>%
        mutate(extinction_risk = length(which(goes_extinct == "TRUE")) / reps) %>%
        ungroup() %>%
        select(beta, p, extinction_risk, mean_local_freq_ext) %>%
        unique() %>%
        rbind(ext_summary, .)
      
      # all %>%
      #   mutate(unique_id = paste(replicate, beta)) %>%
      #   ggplot(aes(x = time, y = N_global, colour = beta, group = unique_id)) + geom_path()
    }
  
    p_curr = p_curr + 1
  }

  l = l+1
}


## plot noise colour vs. global extinction risk 
ext_summary %>%
  ggplot(aes(x = beta, y = extinction_risk, colour = p, group = p)) +
  geom_line() + 
  theme_bw() +
  labs(x = "Spectral exponent", y = "Proportion of replicates where\nglobal population undergoes extinction",
       colour = "Spatial\nsynchrony") +
  theme(panel.grid = element_blank())

ggsave(filename = "metapop-global-extinction_no-shift.png", path = "figures/proposal", 
       width = 5, height = 3)
## same results as paper suggested


## plot noise colour versus local extinction risk
ext_summary %>%
  ggplot(aes(x = beta, y = mean_local_freq_ext, colour = p, group = p)) +
  geom_line() + 
  theme_bw() +
  labs(x = "Spectral exponent", y = "Mean number of times\na local population\nexperiences extinction",
       colour = "Spatial\nsynchrony") +
  theme(panel.grid = element_blank())

ggsave(filename = "metapop-local-extinction_no-shift.png", path = "figures/proposal", 
       width = 5, height = 3)
## red noise decreases extinction risk only when conditions are spatially asynchronous 

## I thought white noise would have lower local extinction risk than red noise (based on idea that more autocorrelation = higher chance of long bad period leading to extinction), but no
# - maybe amplitude of white noise is too high to see this effect (C&Y ran loooong simulations with smaller temporal fluctuations)

# interpret as: 
# - when a population is not close to extinction, over long time scales red noise increases extinction risk because it increases the chances a population will experience a long bout of bad conditions 
# - when a population is close to extinction, red noise can actually decrease extinction risk compared to white noise  


## now:
## - try adding more populations / changing steepness of suitability gradient 
## - visualize what's happening in some of the sims using raster brick



## things to check: 
## is extinction risk higher at the edges, where suitability is lower?

## rasterize and plot the lattice 
rast = brick(lattice_N_it)
plot(rast[[9]])

animate(rast)



xyz = reshape2::melt(lattice_N_it)
colnames(xyz) <- c("y", "x", "z", "n")

xyz %>%
  filter(z == 3) %>%
  ggplot(aes(x = x, y = y, fill = n)) +
  geom_raster() + 
  coord_fixed()

## look for evidence of latitudinal gradient in suitability 
xyz %>%
  group_by(y) %>%
  mutate(mean = mean(n)) %>%
  distinct() %>%
  ggplot(aes(x = y, y = mean)) +
  geom_point()


## FOR LATER ##
## calculate range shift parameters
## create dataframe of points 
xyz = reshape2::melt(lattice_N_it)
colnames(xyz) <- c("y", "x", "z")

## turn into presence/absence using suitablity threshold of 3
prab <- xyz %>%
  mutate(PrAb = ifelse(z > 3, "1", "0")) 

## find mean lat
prs <- prab %>% 
  filter(PrAb == "1") 

centroid <- mean(prs$y)

## find mean min lat
mins <- prs %>%
  arrange(y) %>%
  .[1:10,] 

trailing <- mean(mins$y)

## find mean max lat
maxs <- prs %>%
  arrange(y) %>%
  arrange(-y) %>%
  .[1:10,] 

leading <- mean(maxs$y)


## define the range edges and centroid
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## play with sensitivity to definitions later (see Jake's spreadsheet)
## for now use:
## centroid = mean latitude of occupied cells 
## trailing edge = mean latitude of 10 southernmost occupied cells
## leading edge = mean latitude of 10 northernmost occupied cells
## defining "occupied cell" as cell in which suitability > 3

## turn into presence/absence using suitablity threshold of 3
prab <- xyz %>%
  mutate(PrAb = ifelse(z > 3, "1", "0")) 

prab %>%
  ggplot(aes(x = x, y = y, fill = PrAb)) +
  geom_raster() + 
  coord_fixed()

## find mean lat
prs <- prab %>% 
  filter(PrAb == "1") 

centroid <- mean(prs$y)

## find mean min lat
mins <- prs %>%
  arrange(y) %>%
  .[1:10,] 

trailing <- mean(mins$y)

## find mean max lat
maxs <- prs %>%
  arrange(y) %>%
  arrange(-y) %>%
  .[1:10,] 

leading <- mean(maxs$y)

## plot in versus out of range 
prab %>%
  ggplot(aes(x = x, y = y, fill = z)) +
  geom_raster() + 
  geom_hline(yintercept = leading, colour = "red") +
  geom_hline(yintercept = trailing, colour = "red") +
  geom_hline(yintercept = centroid, colour = "red") +
  coord_fixed()

## are range limits shifting?
stats %>%
  gather(key = "range_edge", value = "latitude", c(centroid, trailing_edge, leading_edge)) %>%
  ggplot(aes(x = time_step, y = latitude, colour = range_edge)) +
  geom_line() +
  geom_point(aes(y = loc_Emax), colour = "black")

## then shift climate conditions for 1500 time steps
## rapid = 0.33 rows / time step
## slow = 0.25 rows / time step
## assessed extinction risk over 30, 50, 100 years (time steps)
## and measure what they did in the paper + the position of the leading edge / trailing edge / range centroid
## variable to keep track of climate optimum
loc_Emax = 600
shift_by = 1 
stats = data.frame()
x = 1
while(x <= 1500) {
  ## calculate climatic suitability  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~
  ## shift the climate optimum
  
  ## position optimum climatic conditions as a row on the lattice (Emax)
  lattice_r <- matrix(ncol = 100, nrow = 2700)
  lattice_r[loc_Emax,] = Emax
  
  ## assume that conditions decline sigmoidally away from this optimum in both directions
  lattice_r[(1+loc_Emax-600):(loc_Emax-1),] = Emax*(1:599 / (1:599 + h^s))
  lattice_r[(loc_Emax+1):(loc_Emax + 600),] = Emax*rev(1:600 / (1:600 + h^s))
  lattice_r[is.na(lattice_r)] <- Emax*(1 / (1 + h^s))
  # plot(x = c(1:2700), y = c(lattice_r[1:2700,1]))
  
  ## add the noise:
  lattice_r <- r*(lattice_r + noise_array[,,x+500])
  
  ## calculate number of offspring produced per individual  
  all_offspring = foreach (i=1:length(lattice_N_it), .combine=rbind)  %dopar% {
    
    ## get number of individuals in cell
    num_ind = lattice_N_it[i]
    
    ## get current growth rate in cell
    curr_r = lattice_r[i]
    
    # for each individual in the cell, add up their offspring
    num_off <- 0
    for(j in 1:lattice_N_it[i]) {
      num_off <- num_off + sample_distribution(N_it = num_ind, r = curr_r, K = K)
    }
    ## save number of offspring 
    num_off
  }
  
  ## individual dies after reproducing
  ## (so population size beocmes number of offspring)
  lattice_N_it[,] <- all_offspring
  
  ## fixed proportion of offspring disperse into any unoccupied patch from nearest 8 neighbouring patches
  lattice_disp <- matrix(ncol = 100, nrow = 2700)
  lattice_disp[1:2700,1:100] <- 0
  lattice_non_disp <- lattice_N_it
  for(i in 1:nrow(lattice_N_it)) {
    for(j in 1:ncol(lattice_N_it)) {
      
      ## if cell isn't empty
      if(lattice_N_it[i,j] != 0) {
        ## find unoccupied cells within 8 neighbours 
        ## must index corners and edges differently 
        ## top left 
        if(i == 1 & j == 1) {
          ind <- rbind(c(0, 1), c(0, 1))
          ncell = 4
        }
        ## top right 
        else if(i == 1 & j == ncol(lattice_N_it)) {
          ind <- rbind(c(0, 1), c(-1, 0))
          ncell = 4
        }
        ## bottom left
        else if(i == nrow(lattice_N_it) & j == 1) {
          ind <- rbind(c(-1, 0), c(0, 1))
          ncell = 4
        }
        ## bottom right 
        else if(i == nrow(lattice_N_it) & j == ncol(lattice_N_it)) {
          ind <- rbind(c(-1, 0), c(-1, 0))
          ncell = 4
        }
        ## left-most column
        else if(j == 1) {
          ind <- rbind(c(-1, 1), c(0, 1))
          ncell = 6
        }
        ## right-most column
        else if(j == ncol(lattice_N_it)) {
          ind <- rbind(c(-1, 1), c(-1, 0))
          ncell = 6
        }
        ## top row
        else if(i == 1) {
          ind <- rbind(c(0, 1), c(-1, 1))
          ncell = 6
        }
        ## bottom row
        else if(i == nrow(lattice_N_it)) {
          ind <- rbind(c(-1, 0), c(-1, 1))
          ncell = 6
        }
        ## rest of the middle cells 
        else {
          ind <- rbind(c(-1, 1), c(-1, 1))
          ncell = 9
        }
        
        unocc <- which(lattice_N_it[(i+ind[1,1]):(i+ind[1,2]), (j+ind[2,1]):(j+ind[2,2])] == 0)
        
        ## if there are unoccupied cells 
        if(length(unocc) != 0) {
          ## get number of dispersers 
          disp = floor(lattice_N_it[i,j]*p)
          
          # ## distribute evenly between unoccupied cells 
          # remainder <- disp%%length(unocc)
          # per_cell <- rep(floor(disp / length(unocc)), length(unocc))
          # 
          # ## add remainder of dispersers to random cell 
          # per_cell[sample(1:length(unocc), 1)] <- per_cell[sample(1:length(unocc), 1)] + remainder
          # 
          # ## check that is adds up 
          # disp == sum(per_cell)
          # 
          # ## create new vector of # of dispersers per surrounding cell
          # new <- rep(0, ncell)
          # new[unocc] <- per_cell
          
          ## send all dispersers to one unoccupied cell at random
          position = sample(unocc, size = 1)
          
          ## create new vector of # of dispersers per surrounding cell
          new <- rep(0, ncell)
          new[position] <- disp
          
          ## add dispersers to disperser lattice
          lattice_disp[(i+ind[1,1]):(i+ind[1,2]), (j+ind[2,1]):(j+ind[2,2])] <- new
          
          ## update population size in the cell
          lattice_non_disp[i,j] <- lattice_N_it[i,j] - disp
        }
        else {
          ## update population size in the cell
          lattice_non_disp[i,j] <- lattice_N_it[i,j] 
        }
      }
    }
  }
  
  ## add disperser count and non-disperser count to get population sizes at next time step:
  lattice_N_it <- lattice_disp + lattice_non_disp
  
  ## calculate range shift parameters
  ## create dataframe of points 
  xyz = reshape2::melt(lattice_N_it)
  colnames(xyz) <- c("y", "x", "z")
 
  ## turn into presence/absence using suitablity threshold of 3
  prab <- xyz %>%
    mutate(PrAb = ifelse(z > 3, "1", "0")) 
  
  ## find mean lat
  prs <- prab %>% 
    filter(PrAb == "1") 
  
  centroid <- mean(prs$y)
  
  ## find mean min lat
  mins <- prs %>%
    arrange(y) %>%
    .[1:10,] 
  
  trailing <- mean(mins$y)
  
  ## find mean max lat
  maxs <- prs %>%
    arrange(y) %>%
    arrange(-y) %>%
    .[1:10,] 
  
  leading <- mean(maxs$y)
  
  ## save stats:
  stats <- rbind(stats, data.frame(time_step = x+500,
             loc_Emax = loc_Emax,
             num_occupied = length(which(lattice_N_it != 0)),
             num_extinct = length(which(lattice_N_it == 0)),
             centroid = centroid,
             leading_edge = leading, 
             trailing_edge = trailing))
  
  ## shift climatic optimum 
  loc_Emax = loc_Emax + 1
  
  print(paste("On loop number: ", x, sep = ""))
  x = x + 1
}

## save stats 
# write.csv(stats, "data-processed/range-shift-sims/stats_lattice_N_it_100x2700_p0.1_large_K0.99_amp0.05_disp.csv", 
#         row.names = FALSE)

stats = read.csv("data-processed/range-shift-sims/stats_lattice_N_it_100x2700_p0.1_large_K0.99_amp0.05_disp.csv")

## are range limits shifting?
stats %>%
  gather(key = "range_edge", value = "latitude", c(centroid, trailing_edge, leading_edge)) %>%
  ggplot(aes(x = time_step, y = latitude, colour = range_edge)) +
  geom_line() +
  geom_point(aes(y = loc_Emax), colour = "black")

mod <- lm(centroid ~ time_step,
   data = stats)
summary(mod)
# yes

## is range size changing?
stats %>%
  ggplot(aes(x = time_step, y = num_occupied)) +
  geom_point() +
  geom_smooth(method = "lm")

mod <- lm(num_occupied ~ time_step,
          data = stats)
summary(mod)

## with no noise added, no change in range size and near perfect climate tracking
## next: add white noise (K=0) and should see no extinction/change


## interpretation of first runs:
## - no extinction with white noise (matches expectation)
## - can't really interpret range size; need average of multiple runs 
## - under red noise, range limits shift less variably (consistent increase/decrease in position)
## - under white noise, range limits shift sporadically (highly variable increase/decrease)
## - found stable rear edges (expansion at trailing edge)

## next: 
## - try allowing dispersal to only 1 unoccupied neighbouring patch 

## run 100 times
## - generate 100 time series and save in a file 
## - for K = 0, 0.5, 0.99



