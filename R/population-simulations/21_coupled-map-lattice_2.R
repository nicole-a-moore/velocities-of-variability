## building spatially explicit coupled map lattice 
library(tidyverse)
theme_set(theme_bew())

## cellular lattice of 100 x 600 habitat patches
lattice_N_it = matrix(ncol = 100, nrow = 1200)
lattice_r = matrix(ncol = 100, nrow = 1200)
lattice_E_it = matrix(ncol = 100, nrow = 1200)

r = 1 ## set maximum growth rate 
K = 100 ## set carrying capacity
p = 0.1 ## set prop of offspring dispersing 
## higher proportion dispersing = less pronounced effect of suitability gradient
lattice_r[1:1200,1:100] <- r ## start with growth rate = max growth rate 
lattice_N_it[1:1200,1:100] <- K ## start with population size = carrying capacity 

## work in latitudinal variation in climate
h = 10 ## half-saturation constant that defines distance at which Eit = 0.5
s = 1 ## shape parameter defining direction (negative = negative slope) and shape (s > 1 gives signmoid)
Emax = 1 ## set Emax to 1

## position optimum climatic conditions as a row on the lattice (Emax)
lattice_E_it[600,] = Emax

## assume that conditions decline sigmoidally away from this optimum in both directions
# plot(x = 1:299, y = rev(1:299 / (1:299 + h^s)))
lattice_E_it[1:599,] = Emax*(1:599 / (1:599 + h^s))
lattice_E_it[601:1200,] = Emax*rev(1:600 / (1:600 + h^s))
# plot(x = 1:600, y = lattice_E_it[1:600,1])

## create autcorrelated time series for each grid cell
lattice_ac_it = array(dim = c(1200,100,2000))

K = 0 # set autocorrelation coefficient 

## create function to simulate autoregressive noise
sim_ar_noise <- function(K, length) {
  wt <- rnorm(length, mean = 0, sd = 1)
  
  i = 2 
  ets <- c(wt[1])
  while(i <= length) {
    et0 <- ets[i-1]
    et = K*et0 + wt[i]
    ets <- append(ets, et)
    i = i + 1
  }
  
  ## scale the variance 
  ets <- ets*1/sqrt(var(ets))*sqrt(109.086)
  ets <- ets - mean(ets)
  
  return(ets)
}

for(i in 1:1200) {
  for(z in 1:100) {
    lattice_ac_it[i,z,] <- sim_ar_noise(K, 2000)
  }
}

## calculate environmental variation in r
lattice_r <- r*(lattice_E_it) 

## replicate and add temporal autocorrelation to time series 
lattice_r_array <- replicate(2000, lattice_r)
lattice_r_array <- lattice_r_array + lattice_ac_it

## each individual in patch i at time t produces offspring
## number drawn from poisson distribution with mean:
## mu = r / (1 + |1-r|*Nit/K)
lattice_off = matrix(ncol = 100, nrow = 1200)

sample_distribution <- function(N_it, r, K) {
  off <- sample(rpois(as.numeric(r/(1+abs(1-r)*N_it/K)), n = 1000), size = 1)
  return(off)
}

lattice_saved <- lattice_N_it


## run at stable climate conditions for 500 time steps
## then shift climate conditions for 1500 time steps 
x = 1
while(x <= 500) {
  
  lattice_r_curr <- lattice_r_array[,,x]
  
  ## calculate number of offspring produced per individual  
  for(i in 1:length(lattice_N_it)) {
    
    ## for each individual in the cell, add up offspring
    num_off <- 0
    for(j in 1:lattice_N_it[i]) {
      num_off <- num_off + sample_distribution(N_it = lattice_N_it[i], r = lattice_r_curr[i], K = K)
    }
    
    lattice_off[i] <- num_off
  }
  
  ## individual dies after reproducing
  ## (so population size beocmes number of offspring)
  lattice_N_it <- lattice_off
  
  ## fixed proportion of offspring disperse into any unoccupied patch from nearest 8 neighbouring patches
  lattice_disp <- matrix(ncol = 100, nrow = 1200)
  lattice_disp[1:1200,1:100] <- 0
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
  print(paste("On loop number: ", x, sep = ""))
  x = x + 1
}

## save
## saveRDS(lattice_N_it, "data-processed/range-shift-sims/lattice_N_it_100x1200_p0.1_ac0.rds")
lattice_N_it <- readRDS("data-processed/range-shift-sims/lattice_N_it_100x1200_p0.1_ac0.rds")
## saveRDS(lattice_N_it, "data-processed/range-shift-sims/lattice_N_it_100x1200_p0.7.rds")

## rasterize and plot the lattice 
r <- raster(lattice_N_it)
plot(r)

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
## defining "occupied cell" as cell in which suitability > 10

## turn into presence/absence using suitablity threshold of 10
prab <- xyz %>%
  mutate(PrAb = ifelse(z > 10, "1", "0")) 

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
## and measure what they did in the paper + the position of the leading edge / trailing edge / range centroid
## variable to keep track of climate optimum
loc_Emax = 601
shift_by = 1 
stats = data.frame()
x = 1
while(x <= 1500) {
  ## calculate climatic suitability  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~
  ## position optimum climatic conditions as a row on the lattice (Emax)
  lattice_E_it[loc_Emax,] = Emax
  
  ## assume that conditions decline sigmoidally away from this optimum in both directions
  lattice_E_it[1:(loc_Emax - 1),] = Emax*(1:(loc_Emax - 1) / (1:(loc_Emax - 1) + h^s))
  lattice_E_it[(loc_Emax+1):1200,] = Emax*rev(1:(1200-loc_Emax)/ (1:(1200-loc_Emax) + h^s))
  # plot(x = c(1:1200), y = c(lattice_E_it[1:loc_Emax,1], lattice_E_it[(loc_Emax+1):1200,1]))
  
  ## calculate environmental variation in r
  lattice_r <- r*(lattice_E_it) ## ADD AUTOCORRELATION HERE
  
  ## calculate number of offspring produced per individual  
  for(i in 1:length(lattice_N_it)) {
    
    ## for each individual in the cell, add up offspring
    num_off <- 0
    for(j in 1:lattice_N_it[i]) {
      num_off <- num_off + sample_distribution(N_it = lattice_N_it[i], r = lattice_r[i], K = K)
    }
    
    lattice_off[i] <- num_off
  }
  
  ## individual dies after reproducing
  ## (so population size beocmes number of offspring)
  lattice_N_it <- lattice_off
  
  ## fixed proportion of offspring disperse into any unoccupied patch from nearest 8 neighbouring patches
  lattice_disp <- matrix(ncol = 100, nrow = 1200)
  lattice_disp[1:1200,1:100] <- 0
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
 
  ## turn into presence/absence using suitablity threshold of 10
  prab <- xyz %>%
    mutate(PrAb = ifelse(z > 10, "1", "0")) 
  
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
  
  ## shift climatic optimum 
  loc_Emax = loc_Emax + 1
  
  print(paste("On loop number: ", x, sep = ""))
  x = x + 1
}

## are range limits shifting?
stats %>%
  gather(key = "range_edge", value = "latitude", c(centroid, trailing_edge, leading_edge)) %>%
  ggplot(aes(x = time_step, y = latitude, colour = range_edge)) +
  geom_point() +
  geom_point(aes(y = loc_Emax), colour = "black")

mod <- lm(centroid ~ time_step,
   data = stats)
summary(mod)
# no

## is range size changing?
stats %>%
  ggplot(aes(x = time_step, y = num_occupied)) +
  geom_point() +
  geom_smooth(method = "lm")

mod <- lm(num_occupied ~ time_step,
          data = stats)
summary(mod)
## yes it is increasing






plot(lattice_ac_it[1,1,])

## test
ac <- sim_ar_noise(1, 2000)
not_ac <- sim_ar_noise(0, 2000)

var(ac)
var(not_ac)

plot(ac)
plot(not_ac)


