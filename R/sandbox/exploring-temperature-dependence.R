library(tidyverse)
library(foreach)
library(doParallel)

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
numCores

registerDoParallel(numCores)  

#######################################
###          simple approach         ## 
#######################################
## plot y = exp(r + b*epsilont)
r = 1
b = c(0.05, 0.2, 0.3, 0.4)
eps = seq(-30, 30, by = 0.1)

i = 1
ys <- c()
while (i <= length(b)) {
  y = exp(r + b[i]*eps)
  ys <- append(ys, y)
  i = i+1
}

bs <- rep(b, each = length(eps))
df <- data.frame(ys = ys, bs = bs, eps = rep(eps, 4))

df %>%
  filter(bs =="0.05") %>%
  mutate(bs = as.factor(bs)) %>%
  ggplot(., aes(x = eps, y = ys, colour = bs)) + geom_point() + theme_bw()+
  labs(x = "epsilon t", y = "exp(b*epsilon t)", colour = "b")


##############################################
###          theory-driven approach         ## 
##############################################
########### varying carrying capacity ########
M = 100
Ea = 0.2 ## units: eV
k = 8.617e-5 ## unts: eV K-1
T = 273.15
Kt = M^(-3/4)*exp(Ea/(k*T))

range <- seq(from = -50, to = 40, by = 0.001)
temps <- c(273.15 + range)

Kt = M^(-3/4)*exp(Ea/(k*temps))

data = data.frame(Kt = Kt, temps = temps)

ggplot(data, aes(x = temps, y = Kt)) + geom_point() + theme_bw() +
  labs(x = "Temperature (K)", y = "Carrying capacity") + 
  geom_vline(xintercept = 223.15) + geom_vline(xintercept = 313.15)

## how to use an anomaly for temperature in this equation?
## one option: force populations with linearly detrended time series instead of seasonally detrended 
## also, how to shift curve down?

## vary activation energy of metabolism
M = 50
Ea = c(0.2, 0.4, 0.6) ## units: eV
k = 8.617e-5 ## unts: eV K-1
T = 273.15
Kt = M^(-3/4)*exp(Ea/(k*T))

range <- seq(from = -30, to = 40, by = 0.001)
temps <- c(273.15 + range)

i = 1
Kts <- c()
while (i <= length(Ea)) {
  Kt = M^(-3/4)*exp(Ea[i]/(k*temps))
  Kts <- append(Kts, Kt)
  i = i+1
}

Eas <- rep(Ea, each = length(temps))
df <- data.frame(Kts = Kts, Eas = Eas, temps = rep(temps, 3))

df %>%
  mutate(Eas = as.factor(Eas)) %>%
  ggplot(., aes(x = temps, y = Kts, colour = Eas)) + geom_point() + theme_bw() +
  labs(x = "Temperature (K)", y = "Carrying capacity", colour = "Ea (eV)")

## vary mass parameter 
M = c(10,20, 30, 50, 100)
Ea = 0.2 ## units: eV
k = 8.617e-5 ## unts: eV K-1
T = 273.15
Kt = M^(-3/4)*exp(Ea/(k*T))

range <- seq(from = -30, to = 40, by = 0.001)
temps <- c(273.15 + range)

i = 1
Kts <- c()
while (i <= length(M)) {
  Kt = M[i]^(-3/4)*exp(Ea/(k*temps))
  Kts <- append(Kts, Kt)
  i = i+1
}

Ms <- rep(M, each = length(temps))
df <- data.frame(Kts = Kts, Ms = Ms, temps = rep(temps, 5))

df %>%
  mutate(Ms = as.factor(Ms)) %>%
  ggplot(., aes(x = temps, y = Kts, colour = Ms)) + geom_point() + theme_bw() +
  labs(x = "Temperature (K)", y = "Carrying capacity", colour = "Mass (kg)")


########### varying growth rate ########
## write equation for temperature dependence of little r (equation 11 in paper): 
fitness <- function(temp) {
  TD = ((1/Tr) - (1/temp)) ## 
  bT = b_tr*exp(-((temp-Topt)^2/(2*s^2))) 
  inv_alphaT = alpha*exp(Aa*TD)
  term1 = d_a*exp(Ad_a*TD)
  term2 = (1/(alpha*exp(Aa*TD)))*(LambertW::W(b_tr*alpha*exp(Aa*TD - ((temp-Topt)^2/(2*s^2)) + alpha*exp(Aa*TD)*(d_a*exp(Ad_a*TD) - d_j*exp(Ad_j*TD)))))
  rT = -term1 + term2
  
  return(c(term1, term2, rT, bT, inv_alphaT))
}


## vary parameters:
s = 15
d_j = 0.03 # juvenile mortality at reference temperature 
d_a = 0.07 # adult mortality rate at reference temperature 
Ad_j = 7500 # Arrhenius constant for juvenile mortality  
Ad_a = 16000 # Arrhenius constant for adult mortality  
Aa = -4000 # Arrhenius constant for development
alpha = 1 # age at maturity
b_tr = 50 # average per capita fecundity at reference temperature, height and shape of term 2
Topt = 273 # temperature at which avg fecundity is maximal, shifts term 2 
Tr = 298 # reference temperature, changes height of term curves

## call function over range of temperatures
range <- seq(from = -50, to = 50, by = 0.001)
temps <- c(273.15 + range)
data <- sapply(FUN = fitness, temps)

## plot results
data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps,
                  bT = data[4,], inv_alphaT = data[5,]) %>%
  select(-bT, inv_alphaT) %>%
  gather(key = "term", value = "value", c(term1, term2, rT)) %>%
  mutate(facet = ifelse(term == "rT", "rT", "terms"))

data %>%
  ggplot(., aes(x = temps, y = value, colour = term)) + geom_point(size = 0.1) + theme_bw() +
  labs(x = "Temperature (C)", y = "rm(T) or components of rm(T)") +
  scale_y_continuous(limits = c(0, 2)) +
  facet_wrap(~facet, nrow = 2)  + 
  geom_vline(xintercept = 223.15) + geom_vline(xintercept = 313.15)




########### Figure out what parameters to use: ########
M = 100
Ea = 0.2 ## units: eV
k = 8.617e-5 ## unts: eV K-1
T = 273.15
Kt = M^(-3/4)*exp(Ea/(k*T))

range <- seq(from = -50, to = 40, by = 0.001)
temps <- c(273.15 + range)

Kt = M^(-3/4)*exp(Ea/(k*temps))

data = data.frame(Kt = Kt, temps = temps)

ggplot(data, aes(x = temps, y = Kt)) + geom_point() + theme_bw() +
  labs(x = "Temperature (K)", y = "Carrying capacity") + 
  geom_vline(xintercept = 223.15) + geom_vline(xintercept = 313.15)

## call function over range of temperatures
s = 15
d_j = 0.03 # juvenile mortality at reference temperature 
d_a = 0.07 # adult mortality rate at reference temperature 
Ad_j = 7500 # Arrhenius constant for juvenile mortality  
Ad_a = 16000 # Arrhenius constant for adult mortality  
Aa = -4000 # Arrhenius constant for development
alpha = 1 # age at maturity
b_tr = 50 # average per capita fecundity at reference temperature, height and shape of term 2
Topt = 273 # temperature at which avg fecundity is maximal, shifts term 2 
Tr = 298 # reference temperature, changes height of term curves

range <- seq(from = -50, to = 40, by = 0.001)
temps <- c(273.15 + range)
data2 <- sapply(FUN = fitness, temps)

## plot results
data2 = data.frame(term1 = data2[1,], term2 = data2[2,], rT = data2[3,], temps = temps,
                   bT = data2[4,], inv_alphaT = data2[5,])

data2 %>%
  select(-bT, inv_alphaT) %>%
  gather(key = "term", value = "value", c(term1, term2, rT)) %>%
  filter(term == "rT") %>%
  ggplot(., aes(x = temps, y = value, colour = term)) + geom_point(size = 0.1) + theme_bw() +
  labs(x = "Temperature (K)", y = "rm(T) or components of rm(T)") +
  geom_vline(xintercept = 223.15) + geom_vline(xintercept = 313.15) +
  geom_point(data = data, aes(x = temps, y = Kt/500), inherit.aes = F, size = 0.1) +
  scale_y_continuous(sec.axis = sec_axis(~ . * 500))


#####################################################
###        force models with coloured noise        ## 
#####################################################
#### write functions
one_over_f <- function(beta){
  ## create 1/fB noise as described in Cuddington and Yodzis
  n = rnorm(524288, mean = 0, sd = 1) ## random numbers with 0 mean and unit variance 
  phases <- runif(524288, 0, 2*pi) ## random phases
  f = 1:524288 ## f
  
  a <- n*1/(f^(beta/2)) ### amplitudes = random normally distributed numbers * 1/f^beta/2
  
  complex <- a*cos(phases) + a*sin(phases) ## complex coefficients
  
  dft <- fft(complex, inverse = T) ## inverse fast fourier transform the coefficients to get the temporal noise series
  noise = as.numeric(dft[1:89000]) ## crop the noise series to first 200,000 points
  
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
fitness <- function(temp) {
  TD = ((1/Tr) - (1/temp)) ## 
  bT = b_tr*exp(-((temp-Topt)^2/(2*s^2))) 
  inv_alphaT = alpha*exp(Aa*TD)
  term1 = d_a*exp(Ad_a*TD)
  term2 = (1/(alpha*exp(Aa*TD)))*(LambertW::W(b_tr*alpha*exp(Aa*TD - ((temp-Topt)^2/(2*s^2)) + alpha*exp(Aa*TD)*(d_a*exp(Ad_a*TD) - d_j*exp(Ad_j*TD)))))
  rT = -term1 + term2
  
  return(c(term1, term2, rT, bT, inv_alphaT))
}

## set model parameters:
##### carrying capacity #####
## theory-driven:
M = 100
Ea = 0.2 ## units: eV
k = 8.617e-5 ## unts: eV K-1
T = NA
Kt = M^(-3/4)*exp(Ea/(k*T))

## simple:
b = 0.05

##### growth rate #####
## theory-driven:
s = 15
d_j = 0.03 # juvenile mortality at reference temperature 
d_a = 0.07 # adult mortality rate at reference temperature 
Ad_j = 7500 # Arrhenius constant for juvenile mortality  
Ad_a = 16000 # Arrhenius constant for adult mortality  
Aa = -4000 # Arrhenius constant for development
alpha = 1 # age at maturity
b_tr = 50 # average per capita fecundity at reference temperature, height and shape of term 2
Topt = 273 # temperature at which avg fecundity is maximal, shifts term 2 
Tr = 298 # reference temperature, changes height of term curves

## simple:
b = 0.05

##### Ricker model #####
icp = 0.5
N0 = 100
Tmax = 89000
K0 = 100

##### noise colour #####
colours = seq(0, 3.2, by = 0.1)
col = 1
all_results <- list()
while (col <= length(colours)) {
  
  print(paste("On colour ", col, sep = ""))
  ## run 10 simulations per noise colour:
  all <- foreach (z = 1:10, .combine=rbind)  %dopar% {
    ## generate coloured noise:
    output <- one_over_f(beta = colours[col])
    noise = output[[1]]
    noise <- noise*1/sqrt(var(noise))*sqrt(150) + 273.15
    # plot(ts(noise))
    # max(noise)
    # min(noise)
    
    ## estimate noise colour from a linear regression of power spectrum:
    # L <- length(noise)
    # dft <- fft(noise)/L
    # amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
    # amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
    # freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
    # 
    # ## create periodogram data by squaring amplitude of FFT output
    # spectral <- data.frame(freq = freq, power = amp^2)
    # 
    # spectral %>%
    #   ggplot(aes(x = freq, y = power)) + geom_line() +
    #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
    #   theme_minimal()
    
    ## calculate new carrying capacity according to how it varies with temperature
    ## theory-driven:
    C = 273.15 - mean(noise) ## calculate constant needed 
    K <- M^(-3/4)*exp(Ea/(k*(noise + C))) 
    K[which(K <= 0)] <- 5 ## make sure none are negative
    K <- round(K, 0)
    
    ## simple:
    K_s <- mean(K)*exp(b*noise)
    K_s[which(K_s <= 0)] <- 5 ## make sure none are negative
    ## make variance equal
    K_s = round(K_s*1/sqrt(var(K_s))*sqrt(var(K)), 0) 
    
    ## calculate new growth rate according to how it varies with temperature
    ## theory-driven:
    C = 273.15 - mean(noise) ## calculate constant needed 
    r <- sapply(FUN = fitness, noise + C)[3,]
    #plot(ts(r))
    
    ## simple:
    r_s <- exp(b*(noise + C))
    #plot(ts(r_s))
    
    ## make variance equal
    r_s = r_s*1/sqrt(var(r_s))*sqrt(var(r))
    ## constrain r between 1 and 2
    fac = 0.5
    i=1
    while(max(r_s) > 2) {
      r_s = r_s*1/sqrt(var(r_s))*sqrt(fac - i*0.005)
      i = i+1
    }
    #plot(ts(r_s))
    
    N = N0 ## set starting population size to 100
    N_r = Nt = N_r_s = N_K_s = N_K = N_rK = N_rK_s = N
    Nts <- Nts_r <- Nts_r_s <- Nts_K_s <- Nts_K <- Nts_rK <- Nts_rK_s <- N
    i=1
    while(i < Tmax) {
      
      ## vary r, K, neither and both according to environment
      Nt = N*exp(1.5*(1 - (N/K0)^icp)) ## constant r, constant K
      Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
      
      ## simple r:
      Nt_r_s = N_r_s*exp(1.5*(1 - (N_r_s/K0)^icp) + r_s[i]) # variable r, constant K
      Nt_r_s = sample(rpois(as.numeric(Nt_r_s), n = 1000), size = 1)
      
      ## simple K:
      Nt_K_s = N_K_s*exp(1.5*(1 - (N_K_s/K_s[i])^icp)) # constant r, variable K
      Nt_K_s = sample(rpois(as.numeric(Nt_K_s), n = 1000), size = 1)
      
      ## theory-driven r:
      Nt_K = N_K*exp(1.5*(1 - (N_K/K[i])^icp)) # constant r, variable K
      Nt_K = sample(rpois(as.numeric(Nt_K), n = 1000), size = 1)
      
      ## theory-driven K:
      Nt_r = N_r*exp(r[i]*(1 - (N_r/K0)^icp)) # variable r, constant K
      Nt_r = sample(rpois(as.numeric(Nt_r), n = 1000), size = 1)
      
      ## simple both:
      Nt_rK_s = N_rK_s*exp(1.5*(1 - (N_rK_s/K_s[i])^icp) + r_s[i]) # variable r, variable K
      Nt_rK_s = sample(rpois(as.numeric(N_rK_s), n = 1000), size = 1)
      
      ## theory-driven both:
      Nt_rK = N_rK*exp(r[i]*(1 - (N_rK/K[i])^icp)) # variable r, variable K
      Nt_rK = sample(rpois(as.numeric(N_rK), n = 1000), size = 1)
      
      N = Nt
      N_r = Nt_r
      N_r_s = Nt_r_s
      N_K = Nt_K
      N_K_s = Nt_K_s
      N_rK = Nt_rK
      N_rK_s = Nt_rK_s
      Nts <- append(Nts, Nt)
      Nts_K <- append(Nts_K, Nt_K)
      Nts_K_s <- append(Nts_K_s, Nt_K_s)
      Nts_r <- append(Nts_r, Nt_r)
      Nts_r_s <- append(Nts_r_s, Nt_r_s)
      Nts_rK <- append(Nts_rK, Nt_rK)
      Nts_rK_s <- append(Nts_rK_s, Nt_rK_s)
      
      i = i + 1
    }
    
    ## return population time series, noise time series
    data.frame(sim = z, r = r, r_s = r_s, K = K, K_s = K_s,
      Nts = Nts,  Nts_K = Nts_K, Nts_K_s = Nts_K_s, Nts_r = Nts_r, Nts_r_s = Nts_r_s,
      Nts_rK = Nts_rK, Nts_rK_s = Nts_rK_s,
      l = icp, colour = colours[col], noise = noise)
  }
  
  ## bind results from icp together
  all_results[[col]] <- all
  
  col = col + 1
}

# saveRDS(all_results, "data-processed/temperature-dependent-model-out.rds")


all_results <- readRDS("data-processed/temperature-dependent-model-out.rds")


## analyze results 
## get one colour:
col1 <- all_results[[15]]

## plot 10 simulation results:
col1 %>%
  filter(sim ==1) %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = noise, group = sim)) + geom_line() +
  facet_wrap(~sim)

col1 %>%
  filter(sim ==1) %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = r, group = sim)) + geom_line() +
  facet_wrap(~sim)
col1 %>%
  filter(sim ==1) %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = K, group = sim)) + geom_line() +
  facet_wrap(~sim)

col1 %>%
  filter(sim ==1) %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = r_s, group = sim)) + geom_line() +
  facet_wrap(~sim)
col1 %>%
  filter(sim ==1) %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = K_s, group = sim)) + geom_line() +
  facet_wrap(~sim)


col1 %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = Nts_r, group = sim)) + geom_line() +
  facet_wrap(~sim)
col1 %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = Nts_K, group = sim)) + geom_line() +
  facet_wrap(~sim)
col1 %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = Nts_rk, group = sim)) + geom_line() +
  facet_wrap(~sim)

col1 %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = Nts_r_s, group = sim)) + geom_line() +
  facet_wrap(~sim)
col1 %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = Nts_K_s, group = sim)) + geom_line() +
  facet_wrap(~sim)
col1 %>%
  mutate(time = rep(1:89000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = Nts_rK_s, group = sim)) + geom_line() +
  facet_wrap(~sim)



## overcompensation experiment 
Nt1 = 1
Nt2 = 1
Nt1s <- Nt2s <- Nt1
i=1
while (i < 100) {
  Nt1 <- Nt1*exp(2*(1-(Nt1/100)^0.15))
  Nt2 <- Nt2*exp(2*(1-(Nt2/100)^0.16))  
  Nt1s <- append(Nt1s, Nt1)
  Nt2s <- append(Nt2s, Nt2)
  i = i +1
}
plot(ts(Nt1s))
plot(ts(Nt2s))
max(Nt1s)
max(Nt2s)
## overcompensation depends in icp and also r

## if we limit r from 0-2, then completely undercompensatory growth occurs when:
## Nt+1 = Nt*exp(r*(1 - (Nt/K)^icp))
##  (Nt/K)^icp <= 0.5
## maximum growth when Nt = 1, so:
## 1/K^icp <= 0.5
## icp = ln(0.5)/ln(1/K)


Nt = N*exp(1.5*(1 - (N/K0)^icp))


## look at populations!!!
plot(ts(Nts))
plot(ts(Nts_r))
plot(ts(Nts_r_s))
plot(ts(Nts_K))
plot(ts(Nts_K_s))
plot(ts(Nts_rK))
plot(ts(Nts_rK_s))
plot(ts(noise))

## look at spectral colour of dynamics!!
## estimate noise colour from a linear regression of pwoer spectrum:
L <- length(Nts)
dft <- fft(Nts)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()


## estimate noise colour from a linear regression of pwoer spectrum:
ext <- first(which(Nts_K == 0))
#L <- length(Nts_K[1:ext])
L <- length(Nts_K)
#dft <- fft(Nts_K[1:ext])/L
dft <- fft(Nts_K)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

## estimate noise colour from a linear regression of power spectrum:
ext <- first(which(Nts_K_s == 0))
#L <- length(Nts_K_s[1:ext])
L <- length(Nts_K_s)
#dft <- fft(Nts_K_s[1:ext])/L
dft <- fft(Nts_K_s)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()


## estimate noise colour from a linear regression of pwoer spectrum:
ext <- first(which(Nts_r == 0))
#L <- length(Nts_r[1:ext])
L <- length(Nts_r)
#dft <- fft(Nts_r[1:ext])/L
dft <- fft(Nts_r)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

## estimate noise colour from a linear regression of power spectrum:
ext <- first(which(Nts_r_s == 0))
#L <- length(Nts_r_s[1:ext])
L <- length(Nts_r_s)
#dft <- fft(Nts_r_s[1:ext])/L
dft <- fft(Nts_r_s)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()


## estimate noise colour from a linear regression of pwoer spectrum:
ext <- first(which(Nts_rK == 0))
#L <- length(Nts_rK[1:ext])
L <- length(Nts_rK)
#dft <- fft(Nts_rK[1:ext])/L
dft <- fft(Nts_rK)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()

## estimate noise colour from a linear regression of pwoer spectrum:
ext <- first(which(Nts_rK_s == 0))
#L <- length(Nts_rK_s[1:ext])
L <- length(Nts_rK_s)
#dft <- fft(Nts_rK_s[1:ext])/L
dft <- fft(Nts_rK_s)/L
amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
freq <- 1:(L/2)/L ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 

## create periodogram data by squaring amplitude of FFT output
spectral <- data.frame(freq = freq, power = amp^2)

spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()
