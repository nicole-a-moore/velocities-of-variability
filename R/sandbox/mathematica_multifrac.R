### trying to recreate mathematica script in R
library(gsignal)
library(broom)
library(ncdf4)
library(tidyverse)
library(foreach)
library(doParallel)
filter <- dplyr::filter

## C1 = intermittency exponent (H-q(1))
## when C1 = 0, non-intermittent 

## B = 1+2alpha-K(2)
## K(2) is the 'spectral intermittency correction' - for monofractal series, K(q) = 0, so K(2)=0 and B=1+2H

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
numCores

registerDoParallel(numCores)  # use multicore, set to the number of our cores


### functions:
## Levy
Levy = function (alpha) {
  phi = runif(1, -pi/2, pi/2)
  phi0 = -(pi/2)*((1-abs(1-alpha))/alpha)
  result = sign(alpha-1)*sin(alpha*(phi-phi0))*(cos(phi)*abs(alpha-1))^-(1/alpha)
  result = result*as.complex((cos(phi-alpha*(phi-phi0))/rexp(1,1)))^((1-alpha)/alpha)
  
  return(result)
}
  
## frac 
frac <- function (scal, alpha) {
  lamb = length(scal)
  rat = 2
  ap = 1/(1-1/alpha)
  lcut = pi/2
  lcut2 = lcut/rat
  
  exlambda = exp(1)^(-(scal/lcut)^4)
  ## if exlambda < -200, exlambda = -200
  exlambda <- ifelse(exlambda < -200, -200, exlambda)
  
  exlambda2 = exp(1)^(-(scal/lcut2)^4)
  ## if exlambda2 < -200, exlambda2 = -200
  exlambda2 <- ifelse(exlambda2 < -200, -200, exlambda2)
  
  ## chop
  mod <- ifelse(exlambda %% 1 < 10e-10, 0, exlambda %% 1)
  exlambda <- round(exlambda, 1) + mod  
  mod2 <- ifelse(exlambda2 %% 1 < 10e-10, 0, exlambda2 %% 1)
  exlambda2 <- round(exlambda2, 1) + mod2
  
  sing = scal^(-(1/ap))
  singsmooth = sing*exlambda
  t1 = sum(singsmooth)
  exlambda2 <- ifelse(exlambda2 < 1e-10, 0, exlambda2) ## chop
  singsmooth = sing*exlambda2
  t2 = sum(singsmooth)
  
  A = ((rat^(-(1/alpha)))*t1 - t2) / (rat^(-(1/alpha))-1)
  
  ff = exp(1)^(-(scal/3))
  ## if ff < -200, ff = -200
  ff <- ifelse(ff < -200, -200, ff)
  
  singsmooth = sing*ff
  G = sum(singsmooth)
  a = -A/G
  sing = (sing*(1+a*ff))^(1/(alpha-1))
  
  return(sing)
}

## eps1D
eps1D <- function(lambda, al, C1, switch) {
  
  NDf = ifelse(switch == 0, 1, 0.5)
  Heavi = 1
  s = abs(seq(-(lambda-1), lambda-2, 2))
  
  table <- sapply(rep(al, lambda), FUN = Levy)
  ggen1 = ((C1/NDf)^(1/al))*table ##100
  sing = frac(s, al)
  
  ggen1 = cconv(ggen1, sing) 
  
  ggen1 = ifelse(as.numeric(ggen1) < -200, -200, ggen1)
  ggen1 = exp(1)^ggen1
  ff = mean(ggen1)
  epsa = ggen1/ff
  
  return(epsa)
}

## function to calculate spectral exponent over a time series window
fft_calc <- function(ts_window) {
  
  l <- length(ts_window)
  
  # Fourier transform the time series window: 
  dft <- fft(ts_window)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral <- data.frame(freq = freq, power = amp^2)
  
  return(spectral)
}

### test:
tester <- eps1D(lambda = 83950, al = 1.3, C1 = 0.001, switch = 0)
#plot(ts(tester))
spec <- fft_calc(tester[41976:83950])
## plot spectrum:
spec %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  geom_vline(xintercept = 1/10)


## choose alpha = 1.8 and C1 = 0.1 since apparently typical for geophysical processes
mfs <- foreach (1:500, .combine=cbind)  %dopar% {
  sim <- eps1D(lambda = 83950, al = 1.8, C1 = 0.01, switch = 0)
  
  sim = sim[41976:(83950+41975)]
  sim <- sim[1:3650] #sim[41976:(83950+41975)]
  spec <- fft_calc(sim)
  
  # spec %>%
  #   ggplot(aes(x = freq, y = power)) + geom_line() +
  #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  #   theme_minimal()
  
  low <- spec %>%
    filter(freq < 1/16) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)") %>%
    .$estimate
  
  high <- spec %>%
    filter(freq > 1/16) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)") %>%
    .$estimate
  
  c(sim, low, high)
}

# exps <- data.frame(al = 1.9, c1 = 0.1, lows = mfs[83951, ],
#                    highs = mfs[83952,])
# 
# mfs <- mfs[-c(83951, 83952),]

## short
exps <- data.frame(al = 1.8, c1 = 0.01, lows = mfs[3651, ],
                   highs = mfs[3652,])

mfs <- mfs[-c(3651, 3652),]

## pick one out and perform fft:
multi <- mfs[,1]
plot(ts(multi))

spectral_multi <- fft_calc(multi)

spectral_multi %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

## generate 1/fB noise as described in Cuddington and Yodzis
n = rnorm(524288, mean = 0, sd = 1) ## random numbers with 0 mean and unit variance 
phases <- runif(524288, 0, 2*pi) ## random phases
f = 1:524288 ## f

a <- n*1/(f^(0.5/2)) ### amplitudes = random normally distributed numbers * 1/f^beta/2

complex <- a*cos(phases) + a*sin(phases) ## complex coefficients

dft <- fft(complex, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
noise = as.numeric(dft[1:41975]) ## crop the noise series to first 200,000 points
plot(ts(noise))

## plot spectra
spectral_mono <- fft_calc(noise)

spectral_mono %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

true_colour <- lm(data = spectral, log(power) ~ log(freq))

#### add them
## remove mean and change variance to 4140:
noise <- noise*1/sqrt(var(noise))*sqrt(var(multi))
noise <- noise - mean(noise)

plot(ts(noise + multi))

spectral_comb <- fft_calc(noise + multi)
spectral_comb %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

## plot them all together:
spectral_comb$type = "mono + multi"
spectral_mono$type = "mono"
spectral_multi$type = "multi"

spectral <- rbind(spectral_multi, spectral_mono) 
#%>%rbind(., spectral_multi)

spectral %>%
  ggplot(aes(x = freq, y = power, group = type, colour = type)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()


low <- spectral_multi %>%
  filter(freq < 1/16) %>%
  lm(., formula = log10(power) ~ log10(freq)) %>%
  tidy(.) %>%
  filter(term == "log10(freq)") %>%
  .$estimate

high <- spectral_multi %>%
  filter(freq > 1/16) %>%
  lm(., formula = log10(power) ~ log10(freq)) %>%
  tidy(.) %>%
  filter(term == "log10(freq)") %>%
  .$estimate




## next: 
## what features does real temp data have? how can i best simulate it?

## read in linearlly detrended temp data 
l_open = nc_open("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/l-detrended_lon-0-60_lat--30--90.nc")
l_detrended_tas = ncvar_get(l_open, "var1_1")
nc_close(l_open)

## look at a spectrum:
ts <- l_detrended_tas[1,30,] ##1,30
plot(ts(ts))
spec <- fft_calc(ts[1:3650])

## plot spectrum:
spec %>%
  #filter(freq > 1/16) %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")

spec$type = "real"

spectral %>%
  rbind(., spec) %>%
  filter(type %in% c("real", "multi")) %>%
  ggplot(aes(x = freq, y = power, group = type, colour = type)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()


## what range of slopes are there in real temp data?
spec_exp <- read.csv("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/spec-exp_long-0-60_lat-90-30.csv") %>%
  filter(time_window_width == "10 years")

hist(spec_exp$s_spec_exp_PSD_low)
hist(-exps$lows)

hist(spec_exp$s_spec_exp_PSD_high)
hist(-exps$high)

## next step:
## figure out a set of alpha and C1 value combos that get us the full variety of spectral exponents found in air temp data 




###### garbage  #####
## simulate some ts with the same parameters and see how much high and low frequency exponent varies
c1s <- c(0.000001, 0.00001, 0.0001, 0.001, 0.1)
als <- c(1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9)

all_mfs <- list()
all_exps <- list()
count = 1
x=1
while (x <= length(als)) {
  cur_al <- als[x]
  y = 1
  while(y <= length(c1s)) {
    cur_c1 <- c1s[y]
    mfs <- foreach (1:500, .combine=cbind)  %dopar% {
      sim <- eps1D(lambda = 83950, al = cur_al, C1 = cur_c1, switch = 0)
      
      spec <- fft_calc(sim[41976:83950])
      
      low <- spec %>%
        filter(freq < 1/16) %>%
        lm(., formula = log10(power) ~ log10(freq)) %>%
        tidy(.) %>%
        filter(term == "log10(freq)") %>%
        .$estimate
      
      high <- spec %>%
        filter(freq > 1/16) %>%
        lm(., formula = log10(power) ~ log10(freq)) %>%
        tidy(.) %>%
        filter(term == "log10(freq)") %>%
        .$estimate
      
      c(sim[41976:83950], low, high)
    }
    
    all_exps[[count]] = data.frame(al = cur_al, c1 = cur_c1, lows = mfs[41976,],
                                   highs = mfs[41976,])
    
    mfs <- mfs[-c(41976, 41977),]
    all_mfs[[count]] <- mfs
    count = count + 1
    print(paste0("On x ", x, " and y ", y))
    y = y + 1
  }
  x = x + 1
}

#saveRDS(all_mfs, "data-processed/mf-sims_all-mfs.rds")
#saveRDS(all_exps, "data-processed/mf-sims_all-exps.rds")



