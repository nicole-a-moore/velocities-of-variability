library(tidyverse)
library(foreach)
library(doParallel)
library(broom)
theme_set(theme_bw())

## detect cores
starts <- rep(100, 40)
fx <- function(nstart) kmeans(Boston, 4, nstart=nstart)
numCores <- detectCores()
numCores

registerDoParallel(numCores)

## simple loop to test pop dynamics parameters
Nt <- 100
Nts <- 100
for(i in 1:100) {
  Nt <- Nt*exp(-1*(1 - (Nt/1000)^0.1)) ## constant r, constant K
  #Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
  Nts <- append(Nts,Nt)
}
plot(ts(Nts))


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
range <- seq(from = -150, to = 50, by = 0.001)
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
  scale_y_continuous(limits = c(-1, 2)) +
  facet_wrap(~facet, nrow = 2)  + 
  geom_vline(xintercept = 223.15) + geom_vline(xintercept = 313.15)




########### Figure out what parameters to use: ########
M = 50
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
M = 50 ## 20,1.5
Ea = 0.2  ## units: eV
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
Tmax = 200000 
K0 = 100
l = seq(0.1, 1, by = 0.1)

##### noise colour #####
colours = seq(0, 3.2, by = 0.1)

#### run simulations #####
##### keep track of population size for first 10 pops, but only persistence time for the rest 
icp = 1
all_results <- list()
while (icp <= length(l)) {
  
  col = 1
  all_col <- list()
  while (col <= length(colours)) {
    
    print(paste("On colour ", col, " and icp ", icp,  sep = ""))
    ## run 100 simulations per noise colour:
    all <- foreach (z = 1:100, .combine=rbind)  %dopar% {
      ## generate coloured noise:
      
      output <- one_over_f(beta = colours[col])
      noise = output[[1]]  
      #plot(ts(noise))
      
      noise <- noise*1/sqrt(var(noise))*sqrt(109.09) +  273.15
      
      ## calculate new carrying capacity according to how it varies with temperature
      ## theory-driven:
      C = 273.15 - mean(noise) ## calculate constant needed 
      K <- M^(-3/4)*exp(Ea/(k*(noise + C))) 
      K[which(K <= 0)] <- 5 ## make sure none are negative
      K <- round(K, 0)
      
      # ## simple:
      # K_s <- mean(K)*exp(b*noise)
      # K_s[which(K_s <= 0)] <- 5 ## make sure none are negative
      # ## make variance equal
      # K_s = round(K_s*1/sqrt(var(K_s))*sqrt(var(K)), 0) 
      
      ## calculate new growth rate according to how it varies with temperature
      ## theory-driven:
      C = 273.15 - mean(noise) ## calculate constant needed 
      r <- sapply(FUN = fitness, noise + C)[3,]
      r[which(r<(-1))] <- -1
      ## make sure none are below -2
      #plot(ts(r))
      
      # ## simple:
      # r_s <- exp(b*(noise + C))
      # #plot(ts(r_s))
      
      # ## make variance equal
      # r_s = r_s*1/sqrt(var(r_s))*sqrt(var(r))
      # ## constrain r between 1 and 2
      # fac = 0.5
      # i=1
      # while(max(r_s) > 2) {
      #   r_s = r_s*1/sqrt(var(r_s))*sqrt(fac - i*0.005)
      #   i = i+1
      # }
      # #plot(ts(r_s))
      
      N = N0 ## set starting population size to 100
      N_r = Nt = N_K = N_rK = N_r_sv = N  
      #N_r_s = N_K_s = N_rK_s = N
      Nts <- Nts_r <- Nts_K <- Nts_rK <- Nts_r_sv <- N
      #Nts_K_s <- Nts_r_s <- Nts_rK_s <- N
      i=1
      while(i < Tmax) {
        
        ## vary r, K, neither and both according to environment
        Nt = N*exp(1.5*(1 - (N/K0)^l[icp])) ## constant r, constant K
        Nt = sample(rpois(as.numeric(Nt), n = 1000), size = 1)
        
        # ## simple r:
        # Nt_r_s = N_r_s*exp(1.5*(1 - (N_r_s/K0)^l[icp]) + r_s[i]) # variable r, constant K
        # Nt_r_s = sample(rpois(as.numeric(Nt_r_s), n = 1000), size = 1)
        # 
        # ## simple K:
        # Nt_K_s = N_K_s*exp(1.5*(1 - (N_K_s/K_s[i])^l[icp])) # constant r, variable K
        # Nt_K_s = sample(rpois(as.numeric(Nt_K_s), n = 1000), size = 1)
        
        ## theory-driven K:
        Nt_K = N_K*exp(1.5*(1 - (N_K/K[i])^l[icp])) # constant r, variable K
        Nt_K = sample(rpois(as.numeric(Nt_K), n = 1000), size = 1)
        
        ## theory-driven r:
        Nt_r = N_r*exp(r[i]*(1 - (N_r/K0)^l[icp])) # variable r, constant K
        Nt_r = sample(rpois(as.numeric(Nt_r), n = 1000), size = 1)
        
        # ## simple both:
        # Nt_rK_s = N_rK_s*exp(1.5*(1 - (N_rK_s/K_s[i])^l[icp]) + r_s[i]) # variable r, variable K
        # Nt_rK_s = sample(rpois(as.numeric(N_rK_s), n = 1000), size = 1)
        
        ## theory-driven both:
        Nt_rK = N_rK*exp(r[i]*(1 - (N_rK/K[i])^l[icp])) # variable r, variable K
        Nt_rK = sample(rpois(as.numeric(N_rK), n = 1000), size = 1)
        
        N = Nt
        N_r = Nt_r
        # N_r_s = Nt_r_s
        N_K = Nt_K
        # N_K_s = Nt_K_s
        N_rK = Nt_rK
        # N_rK_s = Nt_rK_s
        Nts <- append(Nts, Nt)
        Nts_K <- append(Nts_K, Nt_K)
        # Nts_K_s <- append(Nts_K_s, Nt_K_s)
        Nts_r <- append(Nts_r, Nt_r)
        # Nts_r_s <- append(Nts_r_s, Nt_r_s)
        Nts_rK <- append(Nts_rK, Nt_rK)
        # Nts_rK_s <- append(Nts_rK_s, Nt_rK_s)
        
        i = i + 1
      }
      
      
      ## return population time series for the first ten series
      if (z %in% 1:10) {
        data.frame(sim = z, r = r, K = K,  Nts = Nts,   Nts_r = Nts_r, Nts_K = Nts_K, Nts_rK = Nts_rK,
                   #r_s = r_s, K_s = K_s, Nts_r_s = Nts_r_s, Nts_K_s = Nts_K_s, 
                   #Nts_rK_s = Nts_rK_s,
                   l = l[icp], colour = colours[col], noise = noise, 
                   max_temp = max(noise), min_temp = min(noise),
                   true_colour = output[[2]])
     }
     # otherwise return just the stats
     else {
       data.frame(sim = z, r = mean(r), K = mean(K),
                  Nts = first(which(Nts == 0)),
                  Nts_r = first(which(Nts_r == 0)),
                  Nts_K = first(which(Nts_K == 0)),
                  Nts_rK = first(which(Nts_rK == 0)),
                  #r_s = r_s, K_s = K_s, Nts_r_s = Nts_r_s, Nts_K_s = Nts_K_s,
                  #Nts_rK_s = Nts_rK_s,
                  l = l[icp], colour = colours[col], noise = NA,
                  max_temp = max(noise), min_temp = min(noise),
                  true_colour = output[[2]])
     }
    }
    
    ## save results
    write.csv(all,
              paste("data-processed/temperature-dependent-models_col-",col, "_icp-", l[icp], 
                    ".rds", 
                    sep = ""), row.names = F)
    
    col = col + 1
  }
  
  icp = icp + 1
}

all %>%
  filter(!sim %in% 1:10) %>%
  mutate(Nts = ifelse(is.na(Nts), 400000, Nts),
         Nts_r = ifelse(is.na(Nts_r), 400000, Nts_r),
         Nts_K = ifelse(is.na(Nts_K), 400000, Nts_K),
         Nts_rK = ifelse(is.na(Nts_rK), 400000, Nts_rK)) %>%
  gather(key = "pop_model", value = "ext_time", c(Nts, Nts_r, Nts_K, Nts_rK)) %>%
  group_by(colour, pop_model) %>%
  mutate(mean_pers_time = mean(ext_time)) %>%
  select(-true_colour) %>%
  unique(.)%>%
  ggplot(., aes(y = mean_pers_time, x = colour, colour = pop_model)) + geom_point() +
  #scale_y_log10() +
  theme_light() +
  labs(x = "Noise colour", y = "Mean persistence time", colour = "Population model") + 
  facet_wrap(~pop_model)

all %>%
  filter(!sim %in% 1:10) %>%
  mutate(Nts = ifelse(is.na(Nts), 400000, Nts),
         Nts_r = ifelse(is.na(Nts_r), 400000, Nts_r),
         Nts_K = ifelse(is.na(Nts_K), 400000, Nts_K),
         Nts_rK = ifelse(is.na(Nts_rK), 400000, Nts_rK)) %>%
  gather(key = "pop_model", value = "ext_time", c(Nts, Nts_r, Nts_K, Nts_rK)) %>%
  ggplot(., aes(y = ext_time, x = colour, colour = pop_model, group = colour)) + geom_boxplot() +
  #scale_y_log10() +
  theme_light() +
  labs(x = "Noise colour", y = "Mean persistence time", colour = "Population model") + 
  facet_wrap(~pop_model) + geom_hline(yintercept=400000)

all %>%
  filter(!sim %in% 1:10) %>%
  mutate(Nts = ifelse(is.na(Nts), 400000, Nts),
         Nts_r = ifelse(is.na(Nts_r), 400000, Nts_r),
         Nts_K = ifelse(is.na(Nts_K), 400000, Nts_K),
         Nts_rK = ifelse(is.na(Nts_rK), 400000, Nts_rK)) %>%
  gather(key = "pop_model", value = "ext_time", c(Nts, Nts_r, Nts_K, Nts_rK)) %>%
  ggplot(., aes(y = ext_time, x = colour, colour = pop_model, group = colour)) + geom_point() +
  #scale_y_log10() +
  theme_light() +
  labs(x = "Noise colour", y = "Mean persistence time", colour = "Population model") + 
  facet_wrap(~pop_model) + geom_hline(yintercept=400000)



## make sure data looks good 
col <- seq(0, 3.2, by = 0.1)
names <- paste("data-processed/temperature-dependent-models_col-", 1:length(col), 
               "_icp-0.1_adapted.rds", sep = "")

data <- c()
i=1
while(i <= length(names)) {
  data <- rbind(data, read.csv(names[i]))
  print(i)
  i=i+1
}
sims=data

## remake plot from c&y but for one ICP
# sims <- filter(data, !sim %in% c(1:10)) %>%
#   mutate(Nts = ifelse(is.na(Nts), 89000, Nts),
#          Nts_r = ifelse(is.na(Nts_r), 89000, Nts_r),
#          Nts_K = ifelse(is.na(Nts_K), 89000, Nts_K),
#          Nts_rK = ifelse(is.na(Nts_rK), 89000, Nts_rK))

# sims %>%
#   ggplot(., aes(x = -true_colour, y = Nts_rK)) + geom_point()

sims %>%
  filter(!sim %in% 1:10) %>%
  mutate(Nts = ifelse(is.na(Nts), 200000, Nts),
         Nts_r = ifelse(is.na(Nts_r), 200000, Nts_r),
         Nts_K = ifelse(is.na(Nts_K), 200000, Nts_K),
         Nts_rK = ifelse(is.na(Nts_rK), 200000, Nts_rK)) %>%
  gather(key = "pop_model", value = "ext_time", c(Nts, Nts_r, Nts_K, Nts_rK)) %>%
  group_by(colour, pop_model) %>%
  mutate(mean_pers_time = mean(ext_time)) %>%
  select(-true_colour) %>%
  unique(.) %>%
  ggplot(., aes(y = mean_pers_time, x = colour, colour = pop_model)) + geom_point() +
  #scale_y_log10() +
  theme_light() +
  labs(x = "Noise colour", y = "Mean persistence time", colour = "Population model") + 
  facet_wrap(~pop_model)

sims %>%
  filter(!sim %in% 1:10) %>%
  mutate(Nts = ifelse(is.na(Nts), 200000, Nts),
         Nts_r = ifelse(is.na(Nts_r), 200000, Nts_r),
         Nts_K = ifelse(is.na(Nts_K), 200000, Nts_K),
         Nts_rK = ifelse(is.na(Nts_rK), 200000, Nts_rK)) %>%
  gather(key = "pop_model", value = "ext_time", c(Nts, Nts_r, Nts_K, Nts_rK)) %>%
  group_by(colour, pop_model) %>%
  mutate(cv = sd(ext_time)/mean(ext_time)) %>%
  ungroup() %>%
  select(-ext_time) %>%
  unique(.) %>%
  ggplot(., aes(y = cv, x = colour, colour = pop_model)) + geom_point()  +
  theme_light() +
  labs(x = "Noise colour", y = "CV persistence time", colour = "Population model") + 
  facet_wrap(~pop_model)

sims %>%
  filter(!sim %in% 1:10) %>%
  mutate(Nts = ifelse(is.na(Nts), 200000, Nts),
         Nts_r = ifelse(is.na(Nts_r), 200000, Nts_r),
         Nts_K = ifelse(is.na(Nts_K), 200000, Nts_K),
         Nts_rK = ifelse(is.na(Nts_rK), 200000, Nts_rK)) %>%
  gather(key = "pop_model", value = "ext_time", c(Nts, Nts_r, Nts_K, Nts_rK)) %>%
  ggplot(., aes(y = ext_time, x = colour, colour = pop_model, group = colour)) + geom_boxplot() +
  #scale_y_log10() +
  theme_light() +
  labs(x = "Noise colour", y = "Mean persistence time", colour = "Population model") + 
  facet_wrap(~pop_model)

sims %>%
  filter(!sim %in% 1:10) %>%
  mutate(Nts = ifelse(is.na(Nts), 200000, Nts),
         Nts_r = ifelse(is.na(Nts_r), 200000, Nts_r),
         Nts_K = ifelse(is.na(Nts_K), 200000, Nts_K),
         Nts_rK = ifelse(is.na(Nts_rK), 200000, Nts_rK)) %>%
  gather(key = "pop_model", value = "ext_time", c(Nts, Nts_r, Nts_K, Nts_rK)) %>%
  group_by(colour, pop_model) %>%
  mutate(mean_pers_time = mean(ext_time),
         mean_max_temp = mean(max_temp)) %>%
  select(-true_colour) %>%
  unique(.) %>%
  ggplot(., aes(x = mean_max_temp, y = mean_pers_time, colour = pop_model)) + geom_point()

data %>%
  filter(sim == 8, colour == 3.1) %>%
  mutate(time = rep(1:200000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = Nts_r, group = sim)) + geom_line() + 
  geom_line(aes(y = r*100), colour = "red")

data %>%
  filter(sim == 2, colour == 0) %>%
  mutate(time = rep(1:200000, length(unique(.$sim)))) %>%
  ggplot(., aes(x = time, y = Nts_K, group = sim)) + geom_line() + 
  geom_line(aes(y = K), colour = "red")



## analyze change in noise colour over time windows 
sims <- filter(data, sim %in% c(1:10))

i = 1
col <- seq(0, 3.2, by = 0.1)
pop_exp_list <- list()
element = 1
## loop through colours
while (i <= length(col)) {
  
  ## loop through pop simulations
  x = 1
  while (x < 11) {
    ## pull out population simulation
    sim <- filter(sims, sim == x, colour == as.factor(col[i]))
    
    if(nrow(sim) == 0) {
      x = x + 1
    }
    else {
      ## for each time window, process time series and calculate colour
      n = 5
      while (n < 11) {
        year_start <- 0
        year_stop <- n*365 
        
        while (year_start <= (89000 - n*365)) {
          ## extract temps within time window
          ts_chunk_r <- sim$Nts_r[year_start:year_stop]
          ts_chunk_K <- sim$Nts_K[year_start:year_stop]
          noise <- sim$noise[year_start:year_stop]
          demsto <- sim$Nts[year_start:year_stop]
          
          ## if the window does not have NAs
          if (length(which(ts_chunk_r == 0)) == 0) {
            ## get length
            L = length(ts_chunk_r)
            
            ## preprocess the time series:
            ## a. subtracting mean
            ts_r <- ts_chunk_r - mean(ts_chunk_r)
            
            ## b. windowing - multiply by a parabolic window 
            window_r <- parabolic_window(series = ts_r, N = L)
            ts_r <- ts_r*window_r
            
            ## c. bridge detrending (endmatching)
            ## ie. subtracting from the data the line connecting the first and last points of the series
            ts_r <- bridge_detrender(windowed_series = ts_r, N = L)
            
            ## calculate spectral exponent in window using PSD and AWC methods
            r_exp_PSD <- spectral_exponent_calculator_PSD(ts_r, l = L)
            
            r_exp_PSD_low <- r_exp_PSD[[1]]
            r_exp_PSD_high <- r_exp_PSD[[2]]
            r_exp_PSD_all <- r_exp_PSD[[3]]
            # r_plot <- r_exp_PSD[[4]]
            # 
            # ggsave(r_plot, 
            #        path = "figures/pop-spectral-change", 
            #        filename = paste(" 0.9_r_", n, "_window-start-", year_start, 
            #                         ".png", sep = ""), 
            #        device = "png", width = 6, height = 4)
          }
          else {
            r_exp_PSD_low = r_exp_PSD_high = r_exp_PSD_all = NA
          }
          if (length(which(ts_chunk_K == 0)) == 0) {
            L = length(ts_chunk_K)
            window_K <- parabolic_window(series = ts_K, N = L)
            ts_K <- ts_chunk_K - mean(ts_chunk_K)
            ts_K <- ts_K*window_K
            ts_K <- bridge_detrender(windowed_series = ts_K, N = L)
            K_exp_PSD <- spectral_exponent_calculator_PSD(ts_K, l = L)
            
            K_exp_PSD_low <- K_exp_PSD[[1]]
            K_exp_PSD_high <- K_exp_PSD[[2]]
            K_exp_PSD_all <- K_exp_PSD[[3]]
            # K_plot <- K_exp_PSD[[4]]
            
            # ggsave(K_plot, 
            #        path = "figures/pop-spectral-change", 
            #        filename = paste(" 0.9_K_", n, "_window-start-", year_start, 
            #                         ".png", sep = ""), 
            #        device = "png", width = 6, height = 4)
          }
          else {
            K_exp_PSD_low = K_exp_PSD_high = K_exp_PSD_all = NA
          }
          if (length(which(noise == 0)) == 0) {
            ## get length
            L = length(noise)
            
            ## preprocess the time series:
            ## a. subtracting mean
            noise_ts <- noise - mean(noise)
            
            ## b. windowing - multiply by a parabolic window 
            window_noise <- parabolic_window(series = noise_ts, N = L)
            noise_ts <- noise_ts*window_noise
            
            ## c. bridge detrending (endmatching)
            ## ie. subtracting from the data the line connecting the first and last points of the series
            noise_ts <- bridge_detrender(windowed_series = noise_ts, N = L)
            
            ## calculate spectral exponent in window using PSD and AWC methods
            noise_exp_PSD <- spectral_exponent_calculator_PSD(noise_ts, l = L)
            
            noise_exp_PSD_low <- noise_exp_PSD[[1]]
            noise_exp_PSD_high <- noise_exp_PSD[[2]]
            noise_exp_PSD_all <- noise_exp_PSD[[3]]
          }
          else {
            noise_exp_PSD_low = noise_exp_PSD_high = noise_exp_PSD_all = NA
          }
          if (length(which(demsto == 0)) == 0) {
            ## get length
            L = length(demsto)
            
            ## preprocess the time series:
            ## a. subtracting mean
            demsto_ts <- demsto - mean(demsto)
            
            ## b. windowing - multiply by a parabolic window 
            window_demsto <- parabolic_window(series = demsto_ts, N = L)
            demsto_ts <- demsto_ts*window_demsto
            
            ## c. bridge detrending (endmatching)
            ## ie. subtracting from the data the line connecting the first and last points of the series
            demsto_ts <- bridge_detrender(windowed_series = demsto_ts, N = L)
            
            ## calculate spectral exponent in window using PSD and AWC methods
            demsto_exp_PSD <- spectral_exponent_calculator_PSD(demsto_ts, l = L)
            
            demsto_exp_PSD_low <- demsto_exp_PSD[[1]]
            demsto_exp_PSD_high <- demsto_exp_PSD[[2]]
            demsto_exp_PSD_all <- demsto_exp_PSD[[3]]
          }
          else {
            demsto_exp_PSD_low = demsto_exp_PSD_high = demsto_exp_PSD_all = NA
          }
          
          ## store:
          pop_exp_list[[element]] <- c(r_exp_PSD_low, r_exp_PSD_high, r_exp_PSD_all,
                                       K_exp_PSD_low, K_exp_PSD_high, K_exp_PSD_all,
                                       noise_exp_PSD_low, noise_exp_PSD_high, noise_exp_PSD_all,
                                       demsto_exp_PSD_low, demsto_exp_PSD_high, demsto_exp_PSD_all,
                                       year_start, year_stop, paste(n*365, "days"),
                                       0.1, col[i], x)
          
          
          
          
          ## move to next window
          year_start = year_stop + 1
          year_stop = year_stop + n*365 
          
          element = element + 1
        }
        
        print(paste("Calculating spectral exponent for population ", x, " window width ", n, 
                    " colour ", col[i], sep = ""))
        
        ## move to next window width
        n = n + 1
      }
      x = x + 1 
    }
  }
  
  i = i + 1
}

## bind rows in list into data frame
pop_exp_df <- data.frame(do.call(rbind, pop_exp_list), stringsAsFactors = FALSE)
colnames(pop_exp_df) <- c("r_spec_exp_PSD_low", "r_spec_exp_PSD_high", "r_spec_exp_PSD_all",
                           "K_spec_exp_PSD_low", "K_spec_exp_PSD_high", "K_spec_exp_PSD_all",
                          "noise_spec_exp_PSD_low", "noise_spec_exp_PSD_high", "noise_spec_exp_PSD_all",
                          "demsto_spec_exp_PSD_low", "demsto_spec_exp_PSD_high", 
                          "demsto_spec_exp_PSD_all",
                           "window_start_year",
                           "window_stop_year", "time_window_width", "icp", "colour", "sim")

write.csv(pop_exp_df, "data-processed/pop-dynamics-colour_icp-0.1.csv", row.names = F)

pop_exp_df <- read.csv("data-processed/pop-dynamics-colour_icp-0.1.csv")

## nice - low spec exp of population time series is correlated with noise colour
pop_exp_df %>%
  group_by(colour, time_window_width, sim) %>%
  mutate(mean = mean(as.numeric(r_spec_exp_PSD_low), na.rm=F)) %>%
  ggplot(., aes(x = colour, y = mean, colour = time_window_width)) + geom_point()

pop_exp_df %>%
  group_by(colour, time_window_width, sim) %>%
  mutate(mean = mean(as.numeric(r_spec_exp_PSD_all), na.rm=F)) %>%
  ggplot(., aes(x = colour, y = mean, colour = time_window_width)) + geom_point()

pop_exp_df %>%
  group_by(colour, time_window_width, sim) %>%
  mutate(mean = mean(as.numeric(r_spec_exp_PSD_high), na.rm=F)) %>%
  ggplot(., aes(x = colour, y = mean, colour = time_window_width)) + geom_point()

pop_exp_df %>%
  group_by(colour, time_window_width, sim) %>%
  mutate(mean = mean(as.numeric(K_spec_exp_PSD_low), na.rm=F)) %>%
  ggplot(., aes(x = colour, y = mean, colour = time_window_width)) + geom_point()

pop_exp_df %>%
  group_by(colour, time_window_width, sim) %>%
  mutate(mean = mean(as.numeric(K_spec_exp_PSD_all), na.rm=F)) %>%
  ggplot(., aes(x = colour, y = mean, colour = time_window_width)) + geom_point()

pop_exp_df %>%
  group_by(colour, time_window_width, sim) %>%
  mutate(mean = mean(as.numeric(K_spec_exp_PSD_high), na.rm=F)) %>%
  ggplot(., aes(x = colour, y = mean, colour = time_window_width)) + geom_point()

## next: does colour of noise in time window correlate with colour of pop dynam in time window?
pop_exp_df %>%
  mutate(colour = as.numeric(as.character(colour)),
         noise_spec_exp_PSD_low = as.numeric(as.character(noise_spec_exp_PSD_low)),
         K_spec_exp_PSD_low = as.numeric(as.character(K_spec_exp_PSD_low))) %>%
  ggplot(., aes(x = noise_spec_exp_PSD_low, y = K_spec_exp_PSD_low, 
                colour = colour)) + geom_point()

pop_exp_df %>%
  mutate(colour = as.numeric(as.character(colour)),
         noise_spec_exp_PSD_all = as.numeric(as.character(noise_spec_exp_PSD_all)),
         K_spec_exp_PSD_all = as.numeric(as.character(K_spec_exp_PSD_all))) %>%
  ggplot(., aes(x = noise_spec_exp_PSD_all, y = K_spec_exp_PSD_all, 
                colour = colour)) + geom_point()

pop_exp_df %>%
  mutate(colour = as.numeric(as.character(colour)),
         noise_spec_exp_PSD_high = as.numeric(as.character(noise_spec_exp_PSD_high)),
         K_spec_exp_PSD_high = as.numeric(as.character(K_spec_exp_PSD_high))) %>%
  ggplot(., aes(x = noise_spec_exp_PSD_high, y = K_spec_exp_PSD_high, 
                colour = colour)) + geom_point()

pop_exp_df %>%
  mutate(colour = as.numeric(as.character(colour)),
         noise_spec_exp_PSD_low = as.numeric(as.character(noise_spec_exp_PSD_low)),
         r_spec_exp_PSD_low = as.numeric(as.character(r_spec_exp_PSD_low))) %>%
  ggplot(., aes(x = noise_spec_exp_PSD_low, y = r_spec_exp_PSD_low, 
                colour = colour)) + geom_point()

pop_exp_df %>%
  mutate(colour = as.numeric(as.character(colour)),
         noise_spec_exp_PSD_high = as.numeric(as.character(noise_spec_exp_PSD_high)),
         r_spec_exp_PSD_high = as.numeric(as.character(r_spec_exp_PSD_high))) %>%
  ggplot(., aes(x = noise_spec_exp_PSD_high, y = r_spec_exp_PSD_high, 
                colour = colour)) + geom_point()

pop_exp_df %>%
  mutate(colour = as.numeric(as.character(colour)),
         noise_spec_exp_PSD_all = as.numeric(as.character(noise_spec_exp_PSD_all)),
         r_spec_exp_PSD_all = as.numeric(as.character(r_spec_exp_PSD_all))) %>%
  ggplot(., aes(x = noise_spec_exp_PSD_all, y = r_spec_exp_PSD_all, 
                colour = colour)) + geom_point()


## plot colour over time for a given sim x colour combo - should be straight line, how much variation?
pop_exp_df %>%
  filter(sim == 1, time_window_width == first(time_window_width)) %>%
  mutate(colour = as.numeric(as.character(colour)),
         noise_spec_exp_PSD_low = as.numeric(as.character(noise_spec_exp_PSD_low)),
         r_spec_exp_PSD_low = as.numeric(as.character(r_spec_exp_PSD_low)),
         window_start_year = as.numeric(as.character(window_start_year))) %>%
  ggplot(., aes(x = window_start_year, y = r_spec_exp_PSD_low, 
                colour = colour, group = colour)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_smooth(aes(y = noise_spec_exp_PSD_low), method = "lm") +
  geom_point(aes(y = noise_spec_exp_PSD_low)) +
  facet_wrap(~colour)

pop_exp_df %>%
  filter(sim == 1, time_window_width == first(time_window_width)) %>%
  mutate(colour = as.numeric(as.character(colour)),
         noise_spec_exp_PSD_low = as.numeric(as.character(noise_spec_exp_PSD_low)),
         K_spec_exp_PSD_low = as.numeric(as.character(K_spec_exp_PSD_low)),
         window_start_year = as.numeric(as.character(window_start_year))) %>%
  ggplot(., aes(x = window_start_year, y = K_spec_exp_PSD_low, 
                colour = colour, group = colour)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_smooth(aes(y = noise_spec_exp_PSD_low), method = "lm") +
  geom_point(aes(y = noise_spec_exp_PSD_low)) +
  facet_wrap(~colour)







##############################################
#####              FUNCTIONS:           ######
##############################################
## function that creates a parabolic window for a given series
## uses Eke 2000, Eq. 6: W(j) = 1 - (2j/(N+1) - 1)^2 for j = 1,...,N
parabolic_window <- function(series, N) {
  j = c(1:N)
  W = c()
  for (i in j) {
    W[i] = 1 - ((2*j[i])/(N+1) - 1)^2
  }
  #plot(W)
  return(W)
}

## function that bridge detrends a windowed series
## calculates line connecting the first and last points of the series
## then subtracts that line from data 
bridge_detrender <- function(windowed_series, N) {
  ## regress to get equation of line:
  data <- data.frame(x = c(1, N), y = windowed_series[c(1,N)])
  eq = lm(y ~ x, data = data) 
  coeffs = eq$coeff
  #plot(windowed_series)
  #abline(a = windowed_series[1], b = coeffs[2])
  
  ## subtract the line from the data
  df <- data.frame(x = 1:N)
  predictions <- predict(eq, df)
  windowed_series = windowed_series - predictions
  #plot(windowed_series) 
  
  return(windowed_series)
}

## function to calculate spectral exponent over a time series window uisng periodogram
spectral_exponent_calculator_PSD <- function(ts_window, l) {
  # Fourier transform the time series window: 
  dft <- fft(ts_window)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral <- data.frame(freq = freq, power = amp^2)
  
  high_freq <- filter(spectral, freq > 1/16)
  low_freq <- filter(spectral, freq < 1/16)
  
  ## plot spectrum:
  # plot <- spectral %>%
  #   ggplot(aes(x = freq, y = power)) + geom_line() +
  #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm", colour = "red") +
  #   theme_light() +
  #   geom_smooth(data = high_freq, method = "lm", colour = "blue") +
  #   geom_smooth(data = low_freq, method = "lm", colour = "green")

  ## get estimate of spectral exponent over time series window:
  ## fit slope to low freqeuncies
  model_output_low <- spectral %>%
    filter(freq < 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  ## fit slope to high frequencies
  model_output_high <- spectral %>%
    filter(freq >= 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  ## fit slope for species with generation times < 1year
  model_output_all <- spectral %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  return(list(-model_output_low$estimate, -model_output_high$estimate, -model_output_all$estimate 
              # ,plot
              ))
}

