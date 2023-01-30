## new C & Y simulations 
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

one_over_f <- function(beta){
  ## create 1/fB noise as described in Cuddington and Yodzis
  n = rnorm(524288, mean = 0, sd = 1) ## random numbers with 0 mean and unit variance 
  phases <- runif(524288, 0, 2*pi) ## random phases
  f = 1:524288 ## f
  
  a <- n*1/(f^(beta/2)) ### amplitudes = random normally distributed numbers * 1/f^beta/2
  
  complex <- a*cos(phases) + a*sin(phases) ## complex coefficients
  
  dft <- fft(complex, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
  noise = as.numeric(dft[1:51830]) ## crop the noise series to first 200,000 points
  
  ## remove mean and change variance to 109.086:
  noise <- noise*1/sqrt(var(noise))*sqrt(109.086)
  noise <- noise - mean(noise)
  
  ## estimate noise colour from a linear regression of pwoer spectrum:
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
  
  true_colour <- lm(data = spectral, log(power) ~ log(freq))
  
  return(list(noise, as.numeric(true_colour$coefficients[2])))
}

## write population model:
popmod <- function(t, # number of time steps
                   N0, # initial population size
                   K # carrying capacity 
){
  N <- numeric(t)
  N[1] <- N0
  for(i in 2:t){
    N[i] <- sample(rpois(as.numeric(N[i-1]*exp(1.5*(1 - (N[i-1]/K[i])^lambda[icp]))), n = 1000), 
                   size = 1)
  }
  return(N)
}

colours = seq(0, 3.2, by = 0.1)
all_results <- list()
lambda = c(0.1, 0.5, 0.9)

K0 = 100
N0 = 100
Tmax = 51830

col = 1
all_col <- list()
while (col <=  length(colours)) {
  
  ## run 100 simulations per model:
  icp = 1
  while (icp <= length(l)) {
    print(paste("On colour ", col, " and icp ", icp, sep = ""))
    
    ## run 100 simulations per colour:
    all <- foreach (z = 1:100, .combine=rbind)  %dopar% {
      
      ## generate coloured noise: 
      output <- one_over_f(beta = colours[col])
      noise = output[[1]]
      
      N = N0 ## set starting population size to 100
      
      K <- K0 + round(noise, digits = 0) ## generate new carrying capacity
      K[which(K <= 0)] <- 5 ## make sure none are negative
      
      ## run pop model
      Nts = popmod(t = 51830, N0 = N0, K = K) # C & Y model - change
      
      data.frame(t = 1:51830, sim = z, K = K, N = Nts, lambda = lambda[icp], 
                 noise = noise, colour = colours[col], true_colour = output[[2]])
    }
    ## write results from colour 
    write.csv(all, paste("data-processed/BerkeleyEarth/simulations/CYsims_icp-", lambda[icp], "_colour-", 
                         colours[col], ".csv", sep = ""), row.names = F)
    
    icp = icp + 1
  }
  ## move to next colour
  col = col + 1
}

## plot a few:
all %>%
  filter(sim %in% 1:10) %>%
  ggplot(., aes(x = t, y = N)) + geom_line() + facet_wrap(~sim)

########################################################
###            calculate extinction risk              ## 
########################################################
## read in 100 simulations for each location and calculate:
## 1) time to extinction (0, 10, 20, 30, 30, 50 individuals)
## 2) time to reach 50% population size
## 3) number of times 50% and 30% pop size is reached
## analyze change in noise colour over time windows 
icp = 1
lambda = c(0.1, 0.5, 0.9)
colours = seq(0, 3.2, by = 0.1)

MTE_all <- data.frame()

## loop through icps
icp = 1
while (icp <= length(l)) {
  
  ## loop through coloours
  col = 1
  while (col <= length(colours)) {
    
    filename = paste("data-processed/BerkeleyEarth/simulations/CYsims_icp-", lambda[icp], "_colour-", 
                     colours[col], ".csv", sep = "")
    
    ## read in data:
    sims <- read.csv(filename)
      
    MTE_50 <- sims %>%
      group_by(sim) %>%
      summarise(Nts_K = first(which(N <= 50)))
    
    MTE_50 <- gather(MTE_50, key = "pop_model", value = "extinction_time", c(Nts_K))
    MTE_50$extinction_metric <- "MTE_50"
    
    
    MTE_40 <- sims %>%
      group_by(sim) %>%
      summarise(Nts_K = first(which(N <= 40)))
    
    MTE_40 <- gather(MTE_40, key = "pop_model", value = "extinction_time", c(Nts_K))
    MTE_40$extinction_metric <- "MTE_40"
    
    
    MTE_30 <- sims %>%
      group_by(sim) %>%
      summarise(Nts_K = first(which(N <= 30)))
    
    MTE_30 <- gather(MTE_30, key = "pop_model", value = "extinction_time", c(Nts_K))
    MTE_30$extinction_metric <- "MTE_30"
    
    
    MTE_20 <- sims %>%
      group_by(sim) %>%
      summarise(Nts_K = first(which(N <= 20)))
    
    MTE_20 <- gather(MTE_20, key = "pop_model", value = "extinction_time", c(Nts_K))
    MTE_20$extinction_metric <- "MTE_20"
    
    MTE_10 <- sims %>%
      group_by(sim) %>%
      summarise(Nts_K = first(which(N <= 10)))
    
    MTE_10 <- gather(MTE_10, key = "pop_model", value = "extinction_time", c(Nts_K))
    MTE_10$extinction_metric <- "MTE_10"
    
    MTE_5 <- sims %>%
      group_by(sim) %>%
      summarise(Nts_K = first(which(N <= 5)))
    
    MTE_5 <- gather(MTE_5, key = "pop_model", value = "extinction_time", c(Nts_K))
    MTE_5$extinction_metric <- "MTE_5"
    
    ## number of times 50 or fewer individuals is reached 
    n50 <- sims %>%
      group_by(sim) %>%
      summarise(Nts_K = length(which(N <= 50)))
    
    n50 <- gather(n50, key = "pop_model", value = "extinction_time", c(Nts_K))
    n50$extinction_metric <- "n50"
    
    ## number of times 50 or fewer individuals is reached 
    n30 <- sims %>%
      group_by(sim) %>%
      summarise(Nts_K = length(which(N <= 30)))
    
    n30 <- gather(n30, key = "pop_model", value = "extinction_time", c(Nts_K))
    n30$extinction_metric <- "n30"
    
    ## bind all:
    MTE <- rbind(MTE_10, MTE_20, MTE_30, MTE_40, MTE_50, MTE_5, n50, n30)
    
    ## add other info
    MTE$lambda <- sims$l[1]
    MTE$colour <- colours[col]
    MTE$true_colour <- rep(unique(sims$true_colour), 8)
    
    ## bind to rest of the data:
    MTE_all <- rbind(MTE_all, MTE)
    
    print(paste("On icp: ", l[icp], " and colour: ", colours[col],
                sep = ""))
    
    col = col + 1
  }
  icp = icp + 1
}

#write.csv(MTE_all, "data-processed/BerkeleyEarth/simulations/CYsims_MTE.csv", row.names = F)
MTE_all <- read.csv("data-processed/BerkeleyEarth/simulations/CYsims_MTE.csv")

MTE_all$source = "simulated"

MTE_forced <- read.csv("data-processed/BerkeleyEarth/MTE_all_0.1.csv")%>%
  rename("lambda" = l) %>%
  select(-lat, -lon) %>%
  mutate(source = "forced")

MTE_all %>%
  ggplot(aes(x = extinction_time, fill = lambda)) + geom_histogram() + 
  facet_wrap(~extinction_metric)

MTE_all %>%
  filter(lambda == 0.1) %>%
  filter(extinction_metric == "n50") %>%
  ggplot(aes(x = extinction_time)) + geom_histogram()  +
  facet_wrap(~pop_model) 

MTE_forced %>%
  filter(extinction_metric == "n50", pop_model == "Nts_K_cy") %>%
  ggplot(aes(x = extinction_time)) + geom_histogram()  +
  facet_wrap(~pop_model) 

## plot together 
MTE_nocol <- select(MTE_all, -colour, - true_colour)
both <- filter(MTE_forced, pop_model == "Nts_K_cy") %>%
  rbind(MTE_nocol, .)

both %>%
  filter(lambda == 0.1) %>%
  filter(extinction_metric == "n50") %>%
  ggplot(aes(x = extinction_time, fill = source)) + geom_histogram()  
## interesting: distributions differ in shape

## plot across colour
spec <- readRDS("data-processed/BerkeleyEarth/BE_noise-colour.rds") 

MTE_gr <- read.csv("data-processed/BerkeleyEarth/MTE_all_0.1.csv") %>%
  rename("lambda" = l) %>%
  filter(pop_model == "Nts_K_cy", lambda == 0.1) %>%
  select(-sim, -pop_model) %>%
  unique() %>%
  mutate(extinction_time = ifelse(is.na(extinction_time), 51830, extinction_time)) %>%
  group_by(extinction_metric, lat, lon) %>%
  mutate(extinction_time = mean(extinction_time, na.rm = T)) %>%
  unique() %>%
  select(-lat, -lon) %>%
  mutate(source = "forced")

ggplot(MTE_gr, aes(x = extinction_time)) + geom_histogram() +
  facet_wrap(~extinction_metric)

join <- left_join(MTE_gr, spec) %>%
  gather(key = "slope_type", value = "spectral_slope", c("s_spec_exp_PSD_high_all","s_spec_exp_PSD_low_all", 
                                                         "s_spec_exp_PSD_all_all")) %>%
  gather(key = "change_slope_type", value = "change_spectral_slope", c("s_estimate_PSD_high","s_estimate_PSD_low", 
                                                         "s_estimate_PSD_all"))

join %>%
  filter(extinction_metric %in% c("n50")) %>%
 # filter(slope_type == "s_exp_PSD_low_all") %>%
  ggplot(aes(x = spectral_slope, y = extinction_time, shape = slope_type, 
             col = change_spectral_slope)) +
  geom_point() +
  facet_wrap(~change_slope_type)

join %>%
  filter(extinction_metric %in% c("n50")) %>%
  filter(slope_type == "s_spec_exp_PSD_low_all") %>%
  filter(change_slope_type == "s_estimate_PSD_low") %>%
  ggplot(aes(x = change_spectral_slope, y = extinction_time, shape = slope_type, 
             col = spectral_slope)) +
  geom_point() 
#+facet_wrap(~change_slope_type)

join %>%
  filter(extinction_metric %in% c("n50")) %>%
  filter(slope_type == "s_spec_exp_PSD_low_all") %>%
  filter(change_slope_type == "s_estimate_PSD_low") %>%
  ggplot(aes(x = change_spectral_slope, y = spectral_slope)) +
  geom_point() 

join <- select(ungroup(join), -lat, -lon)

## calculate avg
MTE_all <- MTE_all %>%
  filter(pop_model == "Nts_K", lambda == 0.1) %>%
  select(-sim, -pop_model) %>%
  unique() %>%
  mutate(extinction_time = ifelse(is.na(extinction_time), 51830, extinction_time)) %>%
  group_by(extinction_metric, colour) %>%
  mutate(extinction_time = mean(extinction_time, na.rm = T)) %>%
  unique() 
  
MTE_all <- select(MTE_all, -colour)

MTE_all <- rename(MTE_all, "spectral_slope" = true_colour)
MTE_all$slope_source = "sim"
MTE_all$spectral_slope = -MTE_all$spectral_slope

join$slope_source <- "real"

rbind(MTE_all, join) %>%
  filter(slope_type == "s_spec_exp_PSD_all_all" | is.na(slope_type)) %>%
  filter(extinction_metric %in% c("n50")) %>%
  ggplot(aes(x = spectral_slope, y = extinction_time, colour = slope_type)) +
  geom_point() +
  labs(x = "Spectral slope", y = "Number of times population size <50 ind.",
       colour = "Spectral slope type") 
  


rbind(MTE_all, join) %>%
  filter(extinction_metric %in% c("n50")) %>%
  ggplot(aes(x = extinction_time, fill = slope_type)) +
  geom_histogram(position = position_dodge()) 


## after accounting for baseline colour, is there still a relationship between extinction risk and spectral change?



all <- do.call("rbind", all_results) 

all <- all %>%
  mutate(unique_run = paste(true_colour, t)) 

unique_ten <- unique(all$unique_run)[1:10]
ten <- filter(all, unique_run %in% unique_ten)

ten %>%
  ggplot(., aes(x = time, y = K)) + geom_line(colour = "red", alpha = 0.5) +
  geom_line(data = ten, aes(x = time, y = N), inherit.aes = F) +
  facet_wrap(~unique_run)



series %>%
  ggplot(., aes(x = time, y = N)) + geom_line() +
  geom_line(data = series, aes(x = time, y = K), colour = "red", alpha = 0.5, inherit.aes = F) +
  theme_light()

## estimate noise colour from a linear regression of pwoer spectrum:
l <- nrow(series)
dft <- fft(series$N)/l
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

## low ICP (undercompensatory) = red 


red <- all %>%
  filter(colour == "2.5") %>%
  filter(l %in% c("0.1", "0.5", "0.9"))

one <- red %>% filter(l == 0.5) 

data %>%
  filter(time <= 365) %>%
  ggplot(., aes(x = time, y = K)) + geom_line() +
  geom_line(data = filter(data, time <= 365), aes(x = time, y = N), colour = "red", alpha = 0.5, inherit.aes = F) +
  theme_light()


