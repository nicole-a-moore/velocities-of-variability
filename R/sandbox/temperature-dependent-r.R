### reproducing A&S results 
library(tidyverse)
theme_set(theme_minimal())

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

#################################################
###                 Figure 1c/d                ## 
#################################################
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

## call function over range of temperatures
range <- seq(from = -13, to = 40, by = 0.001)
temps <- c(273.15 + range)
data <- sapply(FUN = fitness, temps)

## plot results
data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps - 273.15,
                    bT = data[4,], inv_alphaT = data[5,])

basic_plot <- function(data)  {
  bp <- data %>%
    select(-bT, inv_alphaT) %>%
    gather(key = "term", value = "value", c(term1, term2, rT)) %>%
    ggplot(., aes(x = temps, y = value, colour = term)) + geom_point(size = 0.1) + 
    labs(x = "Temperature (C)", y = "rm(T) or components of rm(T)") +
    scale_y_continuous(limits = c(0, 0.7))
  
  return(bp)
}
  
basic_plot(data) +
  ## plot Tmax, Tmin, Topt from the real figure:
  geom_vline(xintercept = 35) + ## Tmax
  geom_vline(xintercept = 13) + ## Tmin
  geom_vline(xintercept = 28)   ## Topt

#################################################
###                 Figure 1e/f                ## 
#################################################
## parameters from figure 1e/f:
s = 4.8

## call function over range of temperatures
range <- seq(from = -13, to = 45, by = 0.001)
temps <- c(273.15 + range)
data <- sapply(FUN = fitness, temps)

## plot results 
data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps - 273.15,
                  bT = data[4,], inv_alphaT = data[5,])

basic_plot(data) +
  ## plot Tmax, Tmin, Topt from the real figure:
  geom_vline(xintercept = 42) + ## Tmax
  geom_vline(xintercept = -1) + ## Tmin
  geom_vline(xintercept = 33)   ## Topt

#################################################
###                 Figure 2a/c                ## 
#################################################
## parameters from figure 2a/c:
Ad_j = 9000 # Arrhenius constant for juvenile mortality  
Ad_a = 10000 # Arrhenius constant for adult mortality  
Aa = -2000 # Arrhenius constant for development

## call function over range of temperatures
range <- seq(from = -13, to = 45, by = 0.001)
temps <- c(273.15 + range)
data <- sapply(FUN = fitness, temps)

## plot results 
data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps - 273.15,
                  bT = data[4,], inv_alphaT = data[5,])

basic_plot(data) +
  ## plot Tmax, Tmin, Topt from the real figure:
  geom_vline(xintercept = 38) + ## Tmax
  geom_vline(xintercept = 12) + ## Tmin
  geom_vline(xintercept = 32.5)   ## Topt

#################################################
###                 Figure 2c/d                ## 
#################################################
## parameters from figure 2c/d:
Aa = -15000 # Arrhenius constant for development

## call function over range of temperatures
range <- seq(from = -13, to = 45, by = 0.001)
temps <- c(273.15 + range)
data <- sapply(FUN = fitness, temps)

## plot results 
data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps - 273.15,
                  bT = data[4,], inv_alphaT = data[5,])

basic_plot(data) +
  ## plot Tmax, Tmin, Topt from the real figure:
  geom_vline(xintercept = 38) + ## Tmax
  geom_vline(xintercept = 12) + ## Tmin
  geom_vline(xintercept = 32.5)   ## Topt

#################################################
###                 Figure 3a/b                ## 
#################################################
## parameters from figure 3a/b:
s = 2
d_j = 0.03 # juvenile mortality at reference temperature 
d_a = 0.05 # adult mortality rate at reference temperature 
Ad_j = 9000 # Arrhenius constant for juvenile mortality  
Ad_a = 10000 # Arrhenius constant for adult mortality  
Aa = -15000 # Arrhenius constant for development
alpha = 60
b_tr = 50 # average per capita fecundity at reference temperature 
Topt = 298 # temperature at which avg fecundity is maximal
Tr = 294 # reference temperature 

## call function over range of temperatures
range <- seq(from = -13, to = 40, by = 0.001)
temps <- c(273.15 + range)
data <- sapply(FUN = fitness, temps)

## plot results 
data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps - 273.15,
                  bT = data[4,], inv_alphaT = data[5,])

basic_plot(data) +
  ## plot Tmax, Tmin, Topt from the real figure:
  geom_vline(xintercept = 33) + ## Tmax
  geom_vline(xintercept = 17) + ## Tmin
  geom_vline(xintercept = 29)   ## Topt


#################################################
###                 Figure 3c/d                ## 
#################################################
## parameters from figure 3c/d:
s = 4

## call function over range of temperatures
range <- seq(from = -13, to = 40, by = 0.001)
temps <- c(273.15 + range)
data <- sapply(FUN = fitness, temps)

## plot results 
data = data.frame(term1 = data[1,], term2 = data[2,], rT = data[3,], temps = temps - 273.15,
                  bT = data[4,], inv_alphaT = data[5,])

basic_plot(data) +
  ## plot Tmax, Tmin, Topt from the real figure:
  geom_vline(xintercept = 43) + ## Tmax
  geom_vline(xintercept = 11) + ## Tmin
  geom_vline(xintercept = 36)   ## Topt


