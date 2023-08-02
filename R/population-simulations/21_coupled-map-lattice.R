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

## cellular lattice of 100 x 2700 habitat patches
lattice_N_it = matrix(ncol = 100, nrow = 2700)
lattice_r = matrix(ncol = 100, nrow = 2700)
lattice_E_it = matrix(ncol = 100, nrow = 2700)

r = 1 ## set maximum growth rate 
K = 100 ## set carrying capacity
p = 0.1 ## set prop of offspring dispersing 
## higher proportion dispersing = less pronounced effect of suitability gradient
lattice_r[1:2700,1:100] <- r ## start with growth rate = max growth rate 
lattice_N_it[1:2700,1:100] <- K ## start with population size = carrying capacity 

## work in latitudinal variation in climate
h = 10 ## half-saturation constant that defines distance at which Eit = 0.5
s = 1 ## shape parameter defining direction (negative = negative slope) and shape (s > 1 gives signmoid)
Emax = 1 ## set Emax to 1

## position optimum climatic conditions as a row on the lattice (Emax)
lattice_E_it[600,] = Emax

## assume that conditions decline sigmoidally away from this optimum in both directions
lattice_E_it[1:599,] = Emax*(1:599 / (1:599 + h^s))
lattice_E_it[601:1200,] = Emax*rev(1:600 / (1:600 + h^s))
lattice_E_it[is.na(lattice_E_it)] <- Emax*(1 / (1 + h^s))
# plot(x = 1:2700, y = lattice_E_it[1:2700,1])

## create autocorrelated time series for each row
lattice_ac_it = array(dim = c(2700,100,2000))

## create function to simulate autoregressive noise
#sim_ar_noise <- function(ac, length) {
#   wt <- rnorm(length, mean = 0, sd = 1)
#   
#   i = 2
#   ets <- c(wt[1])
#   while(i <= length) {
#     et0 <- ets[i-1]
#     et = ac*et0 + wt[i]
#     ets <- append(ets, et)
#     i = i + 1
#   }
#   
#   ## scale the variance
#   ets <- ets*1/sqrt(var(ets))*sqrt(109.086)
#   ets <- ets - mean(ets)
#   
#   return(ets)
# }

sim_ar_noise <- function(ac, length) {
  ## generate random component with mean 0 sd 1
  wt <- rnorm(length, mean = 0, sd = 0.05)
  ## note: changing sd of noise process changes the amplitude
  ## in the paper, they chose 0.05, 0.1, 0.15, 0.2
  
  i = 2
  ets <- c(wt[1])
  while(i <= 200000) {
    et0 <- ets[i-1] ## get previous element in series 
    et = ac*et0 + wt[i] ## calculate next element in series 
    ets <- append(ets, et)
    i = i + 1
  }
  
  ## get long-term mean
  lt_mean = mean(ets, na.rm=TRUE)
  
  ## crop series to desired length 
  ets = ets[1:length]
  
  ## scale the variance
  new_ets = (1/sd(wt))*ets - lt_mean
  
  return(ets)
}

ac = 0 # set autocorrelation coefficient

## environmental noise varies temporally, but not spatially 
## so, simulate one noise time series of length 2000 that all grid cells experience 
env_noise <- sim_ar_noise(length = 2000, ac = ac)

## calculate environmental variation in r:
## make into an array
noise_array <- array(dim = c(2700,100,2000))
env_noise_repped<- rep(env_noise, each = 2700*100)
noise_array[,,] = env_noise_repped
plot(noise_array[1,1,])

## saveRDS(noise_array, "data-processed/range-shift-sims/noise_array_100x2700_p0.1_K0_amp0.05.rds")
noise_array <- readRDS("data-processed/range-shift-sims/noise_array_100x2700_p0.1_K1_amp0.05.rds")

## replicate latitudinal gradient 2000 times 
lattice_E_it_array <- replicate(2000, lattice_E_it)
## add latitudinal gradient and temporal autocorrelation 
lattice_E_it_array <- lattice_E_it_array + noise_array
## multiply by r
lattice_r_array <- r*lattice_E_it_array

## each individual in patch i at time t produces offspring
## number drawn from poisson distribution with mean:
## mu = r / (1 + |1-r|*Nit/K)
lattice_off = matrix(ncol = 100, nrow = 2700)

## function to sample number of offspring produced in cell
sample_distribution <- function(N_it, r, K) {
  off <- sample(rpois(as.numeric(r/(1+abs(1-r)*N_it/K)), n = 1000), size = 1)
  ## note: if r is negative, returns NA
  if(is.na(off)){
    off = 0
  }
  return(off)
}

lattice_saved <- lattice_N_it
lattice_N_it <- lattice_saved

## run at stable climate conditions for 500 time steps
## then shift climate conditions for 1500 time steps 
x = 1 
while(x <= 500) {
  
  ## get grid of growth rates at time x
  lattice_r_curr <- lattice_r_array[,,x]
  
  ## calculate number of offspring produced per individual  
  all_offspring = foreach (i=1:length(lattice_N_it), .combine=rbind)  %dopar% {
    
    ## get number of individuals in cell
    num_ind = lattice_N_it[i]
    
    ## get current growth rate in cell
    curr_r = lattice_r_curr[i]
    
    # for each individual in the cell, add up their offspring
    num_off <- 0
    for(j in 1:lattice_N_it[i]) {
      num_off <- num_off + sample_distribution(N_it = num_ind, r = curr_r, K = K)
    }
    ## save number of offspring 
    num_off
  }
  
  ## each individual dies after reproducing, so population size now becomes number of offspring)
  lattice_N_it[,] <- all_offspring
  
  ## fixed proportion of offspring disperse into any unoccupied patch from nearest 8 neighbouring patches
  ## create grid to hold number of disperses, set to 0
  lattice_disp <- matrix(ncol = 100, nrow = 2700)
  lattice_disp[1:2700,1:100] <- 0
  lattice_non_disp <- lattice_N_it
  for(i in 1:nrow(lattice_N_it)) {
    for(j in 1:ncol(lattice_N_it)) {
      
      ## if cell has offspring 
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
        
        ## if there are unoccupied cells within neighbours 
        if(length(unocc) != 0) {
          ## get number of dispersers 
          disp = ceiling(lattice_N_it[i,j]*p)
          
          ## distribute evenly between unoccupied cells 
          remainder <- disp%%length(unocc)
          per_cell <- rep(floor(disp / length(unocc)), length(unocc))
          
          ## add remainder of dispersers to random empty neighbor cells 
          per_cell[sample(1:length(unocc), 1)] <- per_cell[sample(1:length(unocc), 1)] + remainder
          
          ## check that it adds up 
          disp == sum(per_cell)
          
          ## create new vector of # of dispersers per surrounding cell
          new <- rep(0, ncell)
          new[unocc] <- per_cell
          
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
  stats <- rbind(stats, data.frame(time_step = x,
                                   loc_Emax = loc_Emax,
                                   num_occupied = length(which(lattice_N_it != 0)),
                                   num_extinct = length(which(lattice_N_it == 0)),
                                   centroid = centroid,
                                   leading_edge = leading, 
                                   trailing_edge = trailing))
  
  print(paste("On loop number: ", x, sep = ""))
  x = x + 1
}

## save
## saveRDS(lattice_N_it, "data-processed/range-shift-sims/lattice_N_it_100x2700_p0.1_K0_amp0.05.rds")
lattice_N_it <- readRDS("data-processed/range-shift-sims/lattice_N_it_100x2700_p0.1_K0_amp0.05.rds")

## rasterize and plot the lattice 
rast <- raster(lattice_N_it)
plot(rast)

xyz = reshape2::melt(lattice_N_it)
colnames(xyz) <- c("y", "x", "z")

xyz %>%
  ggplot(aes(x = x, y = y, fill = z)) +
  geom_raster() + 
  coord_fixed()

## look for evidence of latitudinal gradient in suitability 
xyz %>%
  group_by(y) %>%
  mutate(sum = sum(z)) %>%
  select(-z) %>%
  distinct() %>%
  ggplot(aes(x = y, y = sum)) +
  geom_point()
## nice 


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
          
          ## distribute evenly between unoccupied cells 
          remainder <- disp%%length(unocc)
          per_cell <- rep(floor(disp / length(unocc)), length(unocc))
          
          ## add remainder of dispersers to random cell 
          per_cell[sample(1:length(unocc), 1)] <- per_cell[sample(1:length(unocc), 1)] + remainder
          
          ## check that is adds up 
          disp == sum(per_cell)
          
          ## create new vector of # of dispersers per surrounding cell
          new <- rep(0, ncell)
          new[unocc] <- per_cell
          
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
# write.csv(stats, "data-processed/range-shift-sims/stats_lattice_N_it_100x2700_p0.1_large_K0_amp0.05.csv", 
#         row.names = FALSE)

stats = read.csv("data-processed/range-shift-sims/stats_lattice_N_it_100x2700_p0.1_large_K0_amp0.05.csv")

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
## no

## with no noise added, no change in range size and near perfect climate tracking
## next: add white noise (K=0) and should see no extinction/change


## interpretation of first runs:
## - no extinction with white noise (matches expectation)
## - can't really interpret range size; need average of multiple runs 
## - under red noise, range limits shift less variably (consistent increase/decrease in position)
## - under white noise, range limits shift sporadically (highly variable increase/decrease)
## - found stable rear edges (expansion at trailing edge)

## next: 
## run 100 times
## - generate 100 time series and save in a file 
## - for K = 0, 0.5, 0.99



