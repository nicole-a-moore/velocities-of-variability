### replicating Cuddington and Yodzis population model results 
## using multiple cores!!
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
Tmax = 83950 #length of time to run model
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

colours = seq(0, 3.2, by = 0.1)
all_results <- list()
l = seq(0.1, 1, by = 0.1)


icp = 1
while (icp <= length(l)) {
 
  col = 1
  all_col <- list()
  while (col <=  length(colours)) {
    print(paste("On colour ", col, " and icp ", icp, sep = ""))
    
    all <- foreach (1:500, .combine=rbind)  %dopar% {
      N = N0 ## set starting population size to 100
      
      ## generate coloured noise: 
      output <- one_over_f(beta = colours[col])
      noise = output[[1]]
      
      K <- K0 + round(noise, digits = 0) ## generate new carrying capacity
      K[which(K <= 0)] <- 5 ## make sure none are negative
      
      i=1
      while(i <= Tmax) {
        Nt = N*exp(1.5*(1 - (N/K[i])^l[icp]))
        Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
        
        ## stop if the pop goes extinct:
        if (is.na(Nt) | Nt == 0 | i == Tmax) {
          extinction_time = i
          i = Tmax+1
        }
        ## otherwise continue
        else {
          N = Nt   
          i = i + 1
        }
      }
      c(N = N0, K = K[extinction_time], l = l[icp], t = extinction_time, colour = colours[col],
        true_colour = output[[2]])
    }
    
    all_col[[col]] <- all
    results <- c()
    col = col + 1
  }
  
  ## bind results from icp together
  all_results[[icp]] <- do.call("rbind", all_col) 
  
  icp = icp + 1
}

saveRDS(all_results, "data-processed/ricker_parallel_all-results_sandbox.rds")


all_results <- readRDS("data-processed/ricker_parallel_all-results_sandbox.rds")

all <- as.data.frame(do.call("rbind", all_results))

all <- filter(all, t > 1)

all %>%
  ggplot(., aes(x = -true_colour, y = log(t), colour = l)) + geom_point()

all %>%
  ggplot(., aes(x = colour, y = t, colour = l, group = l)) + 
  geom_smooth(se = F) + scale_y_log10()

all %>%
  filter(l != 0) %>%
  group_by(colour, l) %>%
  mutate(mean_pers_time = mean(t))%>%
  select(-true_colour) %>%
  unique(.) %>%
  ggplot(., aes(y = mean_pers_time, x = colour, colour = l)) + geom_point() +
  scale_y_log10() +
  theme_light() +
  labs(x = "Noise colour", y = "Mean persistence time", colour = "Intraspecific\ncompetition\nparameter")


all %>%
  filter(l != 0) %>%
  group_by(colour, l) %>%
  mutate(cv = sd(t)/mean(t)) %>%
  ungroup() %>%
  select(-t) %>%
  unique(.) %>%
  ggplot(., aes(y = cv, x = colour, colour = l)) + geom_point()  +
  theme_light() +
  labs(x = "Noise colour", y = "CV persistence time", colour = "Intraspecific\ncompetition\nparameter")

all %>%
  filter(colour %in% c(0,1,2,3), l %in% c(0.2, 1)) %>%
  select(-true_colour) %>%
  unique(.) %>%
  ggplot(., aes(x = t, fill = l)) + geom_histogram(bins = 10)  + 
  facet_grid(rows = vars(colour), cols = vars(l)) +
  theme_light() +
  labs(x = "Population persistence", y = "No. of populations") +
  theme(legend.position = "none") +
  scale_x_log10(breaks = c(10^1, 10^3, 10^5))
  
  



library(plotly)

all %>%
  filter(l != 0) %>%
  group_by(colour, l) %>%
  mutate(mean_pers_time = mean(t))%>%
  select(-true_colour) %>%
  unique(.) %>% 
  plot_ly(x = .$colour, z = log(.$mean_pers_time),  y= .$l,
          type="scatter3d", mode="markers", color=.$l)






#### gams
library(mgcv)
smooth_interact <- all %>%
  filter(l != 0) %>%
   gam(log(t) ~ s(colour, l), data = .)

smooth_interact_summary <- summary(smooth_interact)
print(smooth_interact_summary$s.table)
plot(smooth_interact, page = 1, scheme = 3)
plot(smooth_interact, page=1, scheme=1, theta = 0, phi = 20) # will give a similar plot to the vis.gam()

vis.gam(smooth_interact, 
        view = c("colour", "l"), 
        theta = 40, 
        n.grid = 500,
        border = NA,
        type = "response")




sims %>%
  filter(l == "0.1") %>%
  ggplot(aes(x = s_spec_exp_PSD_low, y = log(t), colour = s_estimate_PSD_low, group = l)) +
  geom_point() +
  geom_smooth() +
  geom_point(data = as.data.frame(all), aes(x = colour, y = log(t)), inherit.aes = F)


