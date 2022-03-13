library(tidyverse)


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

## parameters from figure 1c/d:
s = 2.5
d_j = 0.03 # juvenile mortality at reference temperature 
d_a = 0.05 # adult mortality rate at reference temperature 
Ad_j = 7500 # Arrhenius constant for juvenile mortality  
Ad_a = 10000 # Arrhenius constant for adult mortality  
Aa = -8000 # Arrhenius constant for development
alpha = 60 # age at maturity
b_tr = 50 # average per capita fecundity at reference temperature, height and shape of term 2
Topt = 298 # temperature at which avg fecundity is maximal, shifts term 2 
Tr = 294 # reference temperature, changes height of term curves

## vary parameters:
s = 4.8
d_j = 0.03 # juvenile mortality at reference temperature 
d_a = 0.05 # adult mortality rate at reference temperature 
Ad_j = 7500 # Arrhenius constant for juvenile mortality  
Ad_a = 10000 # Arrhenius constant for adult mortality  
Aa = -8000 # Arrhenius constant for development
alpha = 60 # age at maturity
b_tr = 50 # average per capita fecundity at reference temperature, height and shape of term 2
Topt = 298 # temperature at which avg fecundity is maximal, shifts term 2 
Tr = 294 # reference temperature, changes height of term curves

## call function over range of temperatures
range <- seq(from = -50, to = 40, by = 0.001)
temps <- c(273.15 + range)
data <- sapply(FUN = fitness, temps)

## plot results
data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps,
                  bT = data[4,], inv_alphaT = data[5,])

data %>%
  select(-bT, inv_alphaT) %>%
  gather(key = "term", value = "value", c(term1, term2, rT)) %>%
  mutate(facet = ifelse(term == "rT", "rT", "terms")) %>%
  mutate(facet = factor(.$facet, levels = c("terms", "rT"), ordered = T)) %>%
  ggplot(., aes(x = temps, y = value, colour = term)) + geom_point(size = 0.1) + theme_bw() +
  labs(x = "Temperature (C)", y = "rm(T) or components of rm(T)") +
  scale_y_continuous(limits = c(0, 0.4)) +
  facet_wrap(~facet, nrow = 2)  + 
  geom_vline(xintercept = 223.15) + geom_vline(xintercept = 313.15)





#### try adjusting time series of coloured noise for complex model:
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

## generate coloured noise: 
output <- one_over_f(beta = 1.5)
noise = output[[1]]
plot(ts(noise))

var(noise)
sd(noise)


noise_new <- noise*1/sqrt(var(noise))*sqrt(150)
plot(ts(noise_new + 273.15))
max(noise_new + 273.15)
min(noise_new + 273.15)


