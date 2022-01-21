## simulations to compare spectral exponent measuring methods 
## inspired by Eke 2000, 2002
library(tidyverse)
library(longmemo)
library(PNWColors)

#~~~~~~~~~~~~~~~~~~~#
#### functions: #####
#~~~~~~~~~~~~~~~~~~~#
## same procedure for calculating fft we are already using:
homemade_fft <- function(time_series) {
  # Fourier transform the simulated time series: 
  l <- length(time_series)
  dft <- fft(time_series)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral <- data.frame(freq = freq, power = amp^2)
  
  return(spectral)
}

## function that creates a parabolic window for a given series
## uses Eke 2000, Eq. 6: W(j) = 1 - (2j/(N+1) - 1)^2 for j = 1,...,N
parabolic_window <- function(series) {
  N = length(series)
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
bridge_detrender <- function(windowed_series) {
  N = length(windowed_series) # get length of series
  
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### simulations - playing around #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## simulate fractional Gaussian noise (fGn)
## of high autocorrelation:
sim <- simARMA0(n = 5000, H = 0.999) # H near 1 - highly self-similar
plot(sim)
beta = 2*(0.999) - 1  ## beta = 2H - 1
beta

# Fourier transform the simulated time series: 
spectral <- homemade_fft(time_series = sim)
spectral$H = 0.999

## of low autocorrelation:
sim <- simARMA0(n = 5000, H = 0.001) # H near 0 - not self-similar
plot(sim)
beta = 2*(0.001) - 1 
beta

# Fourier transform the simulated time series: 
spectral2 <- homemade_fft(time_series = sim)
spectral2$H = 0.001

spectral <- rbind(spectral, spectral2)

spectral %>%
  ggplot(aes(x = freq, y = power, group = H, colour = H)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  theme_minimal()


## try simulating fractional Brownian noise (fBn) from fGn
## from Eke 2000: "when the elements of an fGn or noise signal are cumulatively summed, the resultant series is composed of the elements of a discrete fBm or motion signal"
## also from Eke 2000: noises and motions of same H have B difference of 2

## simulate fGn:
fGn <- simARMA0(n = 50000, H = 0.01)
plot(fGn)

## cumulatively sum the elements:
fBn_1.02 <- ts(cumsum(fGn)) # beta = 2*(0.01) + 1 = 1.02
plot(fBn_1.02)

## repeat for 2 other betas
fGn <- simARMA0(n = 50000, H = 0.5)

fBn_2.0 <- ts(cumsum(fGn)) # beta = 2*(0.5) + 1 = 2 
plot(fBn_2.0)

fGn <- simARMA0(n = 50000, H = 0.99)

fBn_2.98 <- ts(cumsum(fGn)) # beta = 2*(0.99) + 1 = 2.98
plot(fBn_2.98)

fGn_0 <- simARMA0(n = 50000, H = 0.5) # beta = 2*(0.5) - 1 = 0 (random noise)
plot(fGn_0)

fGn_neg <- simARMA0(n = 50000, H = 0.01) # beta = 2*(0.01) - 1 = -0.98 (self dissimilarity)
plot(fGn_neg)


## Fourier transform
spectral <- homemade_fft(time_series = fBn_1.02)
spectral$beta = 1.02

spectral2 <- homemade_fft(time_series = fBn_2.0)
spectral2$beta = 2

spectral3 <- homemade_fft(time_series = fBn_2.98)
spectral3$beta = 2.98

spectral4 <- homemade_fft(time_series = fGn_0)
spectral4$beta = 0

spectral5 <- homemade_fft(time_series = fGn_neg)
spectral5$beta = -0.98

spectral <- rbind(spectral, spectral2, spectral3, spectral4, spectral5)

pal <- pnw_palette(name="Lake", 6,type="discrete")

spectral %>%
  mutate(beta = as.factor(beta)) %>%
  ggplot(aes(x = freq, y = power, group = beta, colour = beta)) + geom_line() +
  scale_y_log10() + scale_x_log10() + 
  #geom_smooth(method = "lm") +
  theme_light() +
  scale_color_manual(values = pal) + 
  labs(x = "Frequency", y = 'Power', colour = "Beta") 

## calculate regression slopes 
m_fGn <- filter(spectral, beta == 0) %>%
  lm(data = ., log10(power) ~ log10(freq))

m_fBn <- filter(spectral, beta == 2) %>%
  lm(data = ., log10(power) ~ log10(freq))
## observed B are not close to true B


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### wavelet analysis - playing around ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(wsyn)

## create sin wave
## a sine wave of amplitude 1 and period 15 that operates for t = 1, . . . , 100 but then disappears
time1<-1:100
time2<-101:200
times<-c(time1,time2)
ts1p1<-sin(2*pi*time1/15)
ts1p2<-0*time2
ts1<-c(ts1p1,ts1p2)
ts<-ts1

plot(ts(ts))

#Then add a sine wave of amplitude 1 and period 8 that operates for t = 101, . . . , 200 but before that is absent.
ts2p1<-0*time1
ts2p2<-sin(2*pi*time2/8)
ts2<-c(ts2p1,ts2p2)
ts<-ts+ts2

plot(ts(ts))

##  add normally distributed white noise of mean 0 and standard deviation 0.5.
ts3<-rnorm(200,mean=0,sd=0.5)
ts<-ts+ts3
plot(ts(ts))

ts<-cleandat(ts,times,clev=1)
wtres<-wt(ts$cdat,times)
class(wtres)
names(wtres)

## plot power 
wtres$values <- wtres$values^2

plotmag(wtres)
## we can see the oscillations at timescale 15 for the first hundred time steps, and the oscilaltions at timescale 8 for the last 100 time steps, as expected

ggdata <- data.frame(times = rep(1:200, n=77), timescale = rep(wtres$timescales, each = 200),
                     magnitude = as.vector(abs(wtres$values)))

ggplot(data = ggdata, aes(y = timescale, x = times, fill = magnitude)) +
  geom_tile() + scale_y_log10()

## try on our simulated time series:
ts <- simARMA0(n = 1000, H = 0.5)
plot(ts)

ts<-cleandat(ts, times = 1:1000, clev=1)
wtres<-wt(ts$cdat, 1:1000)
class(wtres)
names(wtres)

plotmag(wtres)

## now compare to one with higher autocorrelation 
ts <- simARMA0(n = 1000, H = 0.99)
plot(ts)

ts<-cleandat(ts, times = 1:1000, clev=1)
wtres<-wt(ts$cdat, 1:1000)
class(wtres)
names(wtres)

plotmag(wtres)

## :)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### wavelet analysis - getting serious ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## try the Average Wavelet Coefficient method:
## 1. simulate time series of known beta 
ts <- simARMA0(n = 3650, H = 0.5)
ts <- ts(cumsum(ts)) ## beta = 2
plot(ts)

## 2. compute wavelet transform
ts <- cleandat(ts, times = 1:3650, clev=1)
wavelets <- wt(ts$cdat, 1:3650)
plotmag(wavelets)

## 3. calculate arithmetic mean with respect to the translation coefficient (b)
df <- abs(wavelets$values) ## get magnitude of wavelet coeffs
rownames(df) <- wavelets$times
colnames(df) <- wavelets$timescales
avg <- colMeans(df, na.rm = T) ## get average
avg <- as.numeric(avg)

## recreate plot myself
ggdata <- data.frame(times = rep(1:2^17, n=158), timescale = rep(wavelets$timescales, 
                                                                 each = 2^17),
                     magnitude = as.vector(df))

ggplot(data = ggdata, aes(y = timescale, x = times, fill = magnitude)) +
  geom_tile() + scale_y_log10()

## :) 

## plot average coefficients versus period on a log–log plot
data <- data.frame(avg_wavelet_coeff = avg, period = wavelets$timescales)
ggplot(data = data, aes(x = period, y = avg_wavelet_coeff)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() + 
  theme_light() +
  labs(x = "Period", y = "Average wavelet coefficient")

## 4. calculate beta
## get slope
lm <- lm(log(avg_wavelet_coeff) ~ log(period), 
         data = data)

## slope equal to H + 1/2
H = as.numeric(lm$coefficients[2]) - 1/2
beta = 2*H + 1
beta
## should be 2!

## try Modified Average Wavelet Coefficient Method
## by calculating the median absolut deviation 
## and using it to weight regression points 
# MAD(W(a,b)) = Med( |W(a,b) - Med(W(a,b))| )
df <- wavelets$values
med <- median(df, na.rm = T) ## get median of detail coefficients 
MAD <- abs(df - med)
MAD <- robustbase::colMedians(MAD, na.rm = T) 

weight_data <- data.frame(weight = 1/as.vector(MAD), period = wavelets$timescales)

ggplot(data = weight_data, aes(x = period, y = weight)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10()

data <- left_join(data, weight_data)

## perform weighted regression 
lm = lm(log(avg_wavelet_coeff) ~ log(period), 
        weights = data$weight,
        data = data)

H = as.numeric(lm$coefficients[2]) - 1/2
beta = 2*H + 1
beta

## hmm ... not a better estimator 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### multitaper spectral estimation - playing around ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(multitaper)

## simulate time series of known beta 
ts <- simARMA0(n = 3650, H = 0.5)
ts <- ts(cumsum(ts)) ## beta = 2
plot(ts)

## preprocess
## a. subtracting mean
#plot(ts)
ts <- ts - mean(ts)
#plot(ts)

## b. windowing
## multiply by a parabolic window 
window <- parabolic_window(series = ts)
ts <- ts*window
#plot(ts)

## c. bridge detrending (endmatching)
## ie. subtracting from the data the line connecting the first and last points of the series
ts <- bridge_detrender(windowed_series = ts)
#plot(ts)

## normal Fourier transform:
fft <- homemade_fft(ts) %>%
  filter(freq <= 1/8*max(freq))

fft %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")

## multitaper 
mtm <- spec.mtm(ts, plot = F)

mtm <- data.frame(freq = mtm$freq, power = mtm$spec) %>%
  filter(freq <= max(fft$freq) & freq >= min(fft$freq))

mtm %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") + 
  geom_line(data = fft, aes(x = freq, y = power), colour = "blue") +
  scale_y_log10() + scale_x_log10() + 
  geom_smooth(data = fft, aes(x = freq, y = power), method = "lm")

lm(log(power) ~ log(freq),
   data = fft)

lm(log(power) ~ log(freq),
   data = mtm)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### simulations - getting serious #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
beta <- sample(seq(from = -0.999, to = 2.999, by = 0.001), size = 1000)
H <- ifelse(beta < 1, (beta + 1)/2,  
            ifelse(beta == 0, 0.001, 
                   (beta - 1)/2))
betas_df <- c()
lengths <- c(5*365, 10*365)

## for all values of H and beta 
z = 1
while (z <= length(H)) {
  ## for each time series length
  l = 1
  while (l <= length(lengths)) {
    N = lengths[l]
    
    ## simulate 10 time series and calculate beta using different methods and compare the outcome 
    x = 1
    while (x < 11) {
      ## 1. simulate time series of known beta (spectral exponent)
      ## simulate fGn of length N
      ts <- simARMA0(n = N, H = H[z])
      #plot(ts)
      
      ## if beta >= 1, brownian noise
      if (beta[z] >= 1) {
        ## cumulatively sum the elements to create fBn 
        ts <- ts(cumsum(ts))
        #plot(ts)
      }
      
      ## 2. calculate spectral exponent as per DiCecco & Gouier
      ## - pretend fBn is a seasonally- and linearly-detrended time series from our GCM
      ## perform fft:
      spec <- homemade_fft(time_series = ts)
      
      ## plot PSD:
      # spec %>%
      #   ggplot(aes(x = freq, y = power)) + geom_line() +
      #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
      #   theme_minimal()

      ## calculate beta as slope of PSD across all frequencies:
      regression <- lm(log10(power) ~ log10(freq), data = spec)
      ## save in df
      if (x == 1 & z == 1 & l == 1) {
        betas_df <- data.frame(obs_beta = -as.numeric(regression$coeff[2]), true_beta = beta[z],
                               sim = x, ts_length = N,
                               method = 'original')
      }
      else {
        betas_df <- rbind(betas_df, data.frame(obs_beta = -as.numeric(regression$coeff[2]), 
                                               true_beta = beta[z],
                                               sim = x, ts_length = N,
                                               method = 'original'))
      }
      
      ## 3. calculate spectral exponent as per Eke 2000, 2002
      ## preprocess the time series by:
      ## a. subtracting mean
      #plot(ts)
      pp_ts <- ts - mean(ts)
      #plot(pp_ts)
      
      ## b. windowing
      ## multiply by a parabolic window 
      window <- parabolic_window(series = pp_ts)
      pp_ts <- pp_ts*window
      #plot(pp_ts)
      
      ## c. bridge detrending (endmatching)
      ## ie. subtracting from the data the line connecting the first and last points of the series
      pp_ts <- bridge_detrender(windowed_series = pp_ts)
      #plot(pp_ts)
      
      ## Fourier transform:
      spec <- homemade_fft(time_series = pp_ts)
      
      ## d. removing high-frequency part of the spectrum 
      ## from Eke: remove frequencies < 1/8*fmax
      regression <- spec %>%
        filter(freq < 1/8*max(spec$freq)) %>%
        lm(log10(power) ~ log10(freq), data = .)
      
      betas_df <- rbind(betas_df, 
                        data.frame(obs_beta = -as.numeric(regression$coeff[2]), 
                                   true_beta = beta[z],
                                   sim = x, ts_length = N,
                                   method = 'Eke'))
      
      ## 4. calculate spectral exponent using the Average Wavelet Cofficient method
      ## a. compute wavelet transform
      wavelets <- biwavelet::wt(data.frame(time = 1:N, val = ts), do.sig = F)

      ## b. calculate arithmetic mean with respect to the translation coefficient (b)
      df <- sqrt(wavelets$power) ## get magnitude of wavelet coeffs
      avg <- rowMeans(df, na.rm = T) ## get average
      avg <- as.numeric(avg)
      
      data <- data.frame(avg_wavelet_coeff = avg, period = wavelets$period)
      ## plot average coefficients versus period on a log–log plot
      # ggplot(data = data, aes(x = period, y = avg_wavelet_coeff)) +
      #   geom_point() +
      #   scale_x_log10() + scale_y_log10() +
      #   theme_light() +
      #   labs(x = "Period", y = "Average wavelet coefficient")

      ## c. calculate beta
      ## get slope
      lm <- lm(log10(avg_wavelet_coeff) ~ log10(period),
               data = data)

      ## slope equal to H + 1/2
      H_calc = as.numeric(lm$coefficients[2]) - 1/2
      b = 2*H_calc + 1

      betas_df <- rbind(betas_df,
                        data.frame(obs_beta = b,
                                   true_beta = beta[z],
                                   sim = x, ts_length = N,
                                   method = 'AWC'))
      
      ## multitaper method
      mtm <- spec.mtm(pp_ts, plot = F)
      
      mtm <- data.frame(freq = mtm$freq, power = mtm$spec) %>%
        filter(freq < max(spec$freq) & freq > min(spec$freq))
      
      ## get slope
      lm <- lm(log10(power) ~ log10(freq),
               data = mtm)
      b <- as.numeric(-lm$coefficients[2])
      
      betas_df <- rbind(betas_df,
                        data.frame(obs_beta = b,
                                   true_beta = beta[z],
                                   sim = x, ts_length = N,
                                   method = 'multitaper'))
      
      
      ## move to next simulation
      print(paste("Simulating time series ", x, "of length ", N, " with true beta number ", z))
      x = x+1
    }
    ## move to next time series length
    l = l+1
  }
  ## move to next beta
  z = z + 1
}
#write.csv(betas_df, "data-processed/beta-simulations_with_awc.csv", row.names = F)

pal <- pnw_palette(name = "Sunset", n = 4, 'discrete')

## which method gives most accurate estimate of beta across all true beta?
betas_df %>%
  ggplot(., aes(x = true_beta, y = obs_beta, colour = method)) + 
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = 1) +
  theme_light() + theme(panel.grid = element_blank()) + 
  facet_wrap(~ts_length) +
  labs(x = "True beta", y = "Estimated beta", colour = "Method:") +
  scale_color_manual(values = pal, labels = c("Wavelet", "PSD - preprocessed", "Multitaper", 
                                              "PSD - not preprocessed"))

## mean diff:
betas_df %>%
  mutate(diff = (true_beta - obs_beta)) %>%
  group_by(true_beta, ts_length, method) %>%
  mutate(avg_diff = mean(diff)) %>%
  select(-diff, -sim, -obs_beta) %>%
  ungroup() %>%
  unique() %>%
  ggplot(., aes(y = avg_diff, x = true_beta, colour = method)) + geom_point() +
  theme_light() + theme(panel.grid = element_blank()) + 
  facet_wrap(~ts_length) +
  geom_abline(intercept = 0, slope = 0) +
  labs(x = 'True beta', y = "Average difference between true and estimated beta",
       colour = "Method:") +
  scale_color_manual(values = pal, labels = c("Wavelet", "PSD - preprocessed", "Multitaper", 
                                              "PSD - not preprocessed"))


betas_df %>%
  mutate(diff = (true_beta - obs_beta)) %>%
  ggplot(., aes(y = diff, x = method, fill = method)) + geom_boxplot() +
  stat_summary(fun = "mean", geom = "point") +
  theme_light() + theme(panel.grid = element_blank()) + 
  facet_wrap(~ts_length) +
  labs(y = "Abs(difference between true and estimated beta)", x = "")  +
  scale_fill_manual(values = pal, labels = c("Wavelet", "PSD - preprocessed", "Multitaper", 
                                             "PSD - not preprocessed"))+
  theme(legend.position = 'none') + 
  scale_x_discrete(labels = c("Wavelet", "Preprocessed", "Multitaper", "Not preprocessed"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### what frequency cut-off should we use? ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
beta <- c(0, 0.5, 1, 1.5, 2, 2.5, 2.98)
H <- ifelse(beta < 1, (beta + 1)/2,  (beta - 1)/2)
lengths <- c(5*365, 10*365)
freqs = seq(from = 1/64, to = 32/54, by = 1/64)
## for each unique beta x time series length, simulate 100 time series
z = 1
while (z <= length(H)) {
  ## for each time series length
  l = 1
  while (l <= length(lengths)) {
    N = lengths[l]
    
    ## simulate 100 time series 
    x = 1
    while (x < 100) {
      ## simulate fGn of length N
      ts <- simARMA0(n = N, H = H[z])
      #plot(ts)
      
      ## if beta is greater than or equal to 1, brownian noise 
      if (beta[z] >= 1) {
        ## cumulatively sum the elements to create fBn 
        ts <- ts(cumsum(ts))
        #plot(ts)
      }
     
      
      ## preprocess the time series by:
      ## a. subtracting mean
      #plot(ts)
      pp_ts <- ts - mean(ts)
      #plot(ts)
      
      ## b. windowing
      ## multiply by a parabolic window 
      window <- parabolic_window(series = pp_ts)
      pp_ts <- pp_ts*window
      #plot(pp_ts)
      
      ## c. bridge detrending (endmatching)
      ## ie. subtracting from the data the line connecting the first and last points of the series
      pp_ts <- bridge_detrender(windowed_series = pp_ts)
      #plot(pp_ts)

      ## perform fft:
      spec <- homemade_fft(time_series = pp_ts)
      
      ## calculate beta using different high frequency cutoffs
      y = 1
      while (y < length(freqs)) {
        regression <- spec %>%
          filter(freq < freqs[y]*max(spec$freq)) %>%
          lm(log10(power) ~ log10(freq), data = .)
        
        if (x == 1 & z == 1 & l == 1 & y == 1) {
          betas_df <- data.frame(obs_beta = as.numeric(regression$coeff[2]), 
                                 true_beta = beta[z],
                                 sim = x, ts_length = N,
                                 freq_cutoff = freqs[y])
        }
        else {
          betas_df <- rbind(betas_df, data.frame(obs_beta = as.numeric(regression$coeff[2]), 
                                                 true_beta = beta[z],
                                                 sim = x, ts_length = N,
                                                 freq_cutoff = freqs[y]))
        }
        
        y = y + 1
      }
      
      ## move to next simulation
      print(paste("Simulating time series ",x, "of length ", N, " with true beta number ", z))
      x = x+1
    }
    ## move to next time series length
    l = l+1
  }
  ## move to next beta
  z = z + 1
}

#write.csv(betas_df, "data-processed/cut-off-simulations.csv", row.names = F)

## which frequency cut-off  gives best estimate of beta across all ranges of beta?
unique_betas <- betas_df %>%
  group_by(true_beta, ts_length, freq_cutoff) %>%
  mutate(mean_obs_beta = mean(obs_beta)) %>%
  mutate(sd_obs_beta = sd(obs_beta)) %>%
  ungroup() %>%
  select(-obs_beta, -sim) %>%
  unique() 

pal <- pnw_palette(name="Lake", type="continuous")

unique_betas %>%
  mutate(diff = mean_obs_beta + true_beta) %>%
  mutate(true_beta = as.factor(true_beta)) %>%
  ggplot(., aes(x = freq_cutoff, y = diff, colour = true_beta)) + geom_point() +
  geom_linerange(aes(ymax = diff-sd_obs_beta, ymin = diff+sd_obs_beta)) +
  facet_wrap(~ts_length) +
  theme_light() + 
  labs(x = 'Frequency cut-off point', y = "Difference between true and estimated beta",
       colour = 'True beta') +
  scale_colour_manual(values = pal) +
  geom_vline(xintercept = 1/8) +
  scale_x_continuous(breaks = c(seq(from = 1/32, to = 17/32, by = 2/32)), 
                     labels = c("1/32", '3/32', "5/32",  "7/32", "9/32",
                                "11/32", '13/32', "15/32","17/32")) + 
  theme(axis.text.x = element_text(angle = 75, vjust = 0.5, hjust=1))

betas_df %>%
  filter(ts_length == 3650) %>%
  mutate(diff = obs_beta + true_beta) %>%
  mutate(true_beta = as.factor(true_beta)) %>%
  ggplot(., aes(x = freq_cutoff, y = diff, colour = true_beta)) + geom_point() + 
  labs(x = 'Frequency cut-off point', y = "Difference between true and estimated beta",
       colour = 'True beta') +
  scale_colour_manual(values = pal) +
  geom_vline(xintercept = 1/8)

## plot: for each beta, x = frequency cut-off, y = % of simulations within +/-0.1? With line showing 95%?

betas_df %>%
  filter(ts_length == 3650, true_beta == 1.5) %>%
  mutate(diff = obs_beta + true_beta) %>%
  ggplot(., aes(x = freq_cutoff, y = diff, colour = true_beta)) + geom_point() +
  geom_abline(intercept = 0.1, slope = 0) +
  geom_abline(intercept = -0.1, slope = 0) 

pal <- pnw_palette(name = "Sunset", n = 2, 'discrete')

betas_df %>%
  mutate(diff = obs_beta + true_beta) %>%
  group_by(freq_cutoff, true_beta, ts_length) %>%
  mutate(percent_within_0.1 = 100*(length(which(abs(diff) <= 0.1))/100)) %>%
  ungroup() %>%
  mutate(ts_length = as.factor(ts_length)) %>%
  ggplot(., aes(x = freq_cutoff, y = percent_within_0.1, colour = ts_length)) + geom_point() +
  geom_abline(intercept = 95, slope = 0) +
  facet_wrap(~true_beta) +
  geom_vline(xintercept = 1/8) +
  theme_light() +
  labs(x = 'Frequency cut-off point', y = "% of beta estimates +/-0.1 of true beta", 
       colour = "Time series length") +
  scale_color_manual(values = pal)  +
  scale_x_continuous(breaks = c(seq(from = 1/32, to = 17/32, by = 4/32)), 
                     labels = c("1/32",  "5/32",  "9/32",
                                '13/32',"17/32")) + 
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))

betas_df %>%
  mutate(diff = obs_beta + true_beta) %>%
  group_by(freq_cutoff, true_beta, ts_length) %>%
  mutate(percent_within_0.2 = 100*(length(which(abs(diff) <= 0.2))/100)) %>%
  ungroup() %>%
  mutate(ts_length = as.factor(ts_length)) %>%
  ggplot(., aes(x = freq_cutoff, y = percent_within_0.2, colour = ts_length)) + geom_point() +
  geom_abline(intercept = 95, slope = 0) +
  facet_wrap(~true_beta) +
  geom_vline(xintercept = 1/8) +
  theme_light() +
  labs(x = 'Frequency cut-off point', y = "% of beta estimates +/-0.2 of true beta",
       colour = "Time series length")  +
  scale_color_manual(values = pal)  +
  scale_x_continuous(breaks = c(seq(from = 1/32, to = 17/32, by = 4/32)), 
                     labels = c("1/32",  "5/32",  "9/32",
                                '13/32',"17/32")) + 
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1))



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### how does new methodology change previous results? ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## compare results from one 30deg x 30 deg latitude/longitude chunk

## make colour palette
pal = pnw_palette("Sunset",5, type = "discrete")
pal <- c(pal[1], pal[3], pal[4])

## read in old results and new results
old <- read.csv("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/spec-exp_long-0-60_lat-90-30.csv")
new <- read.csv("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/spec-exp_long-0-60_lat-90-30_new.csv")

## filter to only 10 year time windows 
new <- filter(new, time_window_width == "10 years")
old <- filter(old, time_window_width == "10 years")

## left join
data <- left_join(new, old) %>%
  select(lat, lon, s_spec_exp,  
         s_spec_exp_PSD_low, s_spec_exp_AWC,
         s_spec_exp_PSD_high, 
         window_start_year, window_stop_year) %>%
  gather(key = "Method",  value = "spectral_exponent", c(s_spec_exp, s_spec_exp_PSD_low, s_spec_exp_AWC)) %>%
  mutate(spectral_exponent = ifelse(Method %in% c("s_spec_exp"), 
                                    -spectral_exponent, 
                                    spectral_exponent)) %>% ## change sign of spectral exponent for PSD estimates
  mutate(Method = factor(.$Method, levels = c("s_spec_exp", "s_spec_exp_PSD_low", "s_spec_exp_AWC"),
                         ordered = T))

## compare estimates of spectral exponents from different method 
labs <- c("PSD (unprocessed)", "PSD (processed)",
               "AWC")
names(labs) <- c("s_spec_exp", "s_spec_exp_PSD_low", "s_spec_exp_AWC")

data %>%
  group_by(lat, lon, Method) %>%
  mutate(avg_spec_exp = mean(spectral_exponent)) %>%
  ungroup() %>%
  select(-spectral_exponent) %>%
  ggplot(., aes(x = lon, y = lat, fill = avg_spec_exp)) +
  geom_raster() +
  coord_fixed() +
  theme_minimal() +
  facet_wrap(~Method, labeller = labeller(Method = labs)) + 
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "#e7d8d3",
                                               midpoint = 0) +
  labs(fill = "Avg. spectral exponent", y = "Latitude", x = "Longitude") 

data %>%
  ggplot(., aes(x = Method, y = spectral_exponent, fill = Method)) +
  geom_boxplot() +
  theme_light() + 
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_fill_manual(values = pal) +
  scale_x_discrete(labels = c("PSD (unprocessed)", "PSD (processed)",
                                       "AWC")) +
  labs(y = "Spectral exponent")  

## previous method: higher spectral exponent = much steeper slopes
  

## how does change compare?
data %>%
  mutate(lat_lon_method = paste(lat, lon, Method),
         lat_lon = paste(lat, lon)) %>%
  filter(lat_lon %in% paste(unique(data$lat), unique(data$lon))[1:10]) %>%
  ggplot(., aes(x = window_start_year, y = spectral_exponent, 
                colour = lat, 
                group = lat_lon_method)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme_light() + 
  theme(panel.grid = element_blank()) +
  labs(y = "Spectral exponent", x = "Window start year")

data %>%
  group_by(Method, window_start_year) %>%
  mutate(global_avg = mean(spectral_exponent)) %>%
  ggplot(., aes(x = window_start_year, y = global_avg, 
                colour = Method, group = Method)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme_light() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = pal, labels = c("PSD (unprocessed)", "PSD (processed)",
                              "AWC")) +
  labs(y = "Average exponent across 30x30 degree area", x = "Window start year")

## calculate slopes
data %>%
  group_by(Method, window_start_year) %>%
  mutate(global_avg = mean(spectral_exponent)) %>%
  select(global_avg, window_start_year, Method) %>%
  group_by(Method) %>%
  unique() %>%
  do(tidy(lm(global_avg ~ window_start_year, data = .))) %>%
  filter(term == "window_start_year")



data <- left_join(new, old) %>%
  select(lat, lon, s_estimate, s_estimate_PSD_low, s_estimate_PSD_high, s_estimate_AWC, window_start_year, window_stop_year) %>%
  gather(key = "Method",  value = "slope_of_spec_exp", c(s_estimate, s_estimate_PSD_low, s_estimate_AWC)) %>%
  mutate(slope_of_spec_exp = ifelse(Method %in% c("s_estimate"), 
                                    -slope_of_spec_exp, 
                                    slope_of_spec_exp)) %>% ## change sign of slope for PSD estimates
  mutate(Method = factor(.$Method, levels = c("s_estimate", "s_estimate_PSD", "s_estimate_AWC"),
                         ordered = T)) %>%
  select(-window_start_year, - window_stop_year) %>%
  unique()

labs <- c("PSD (unprocessed)", "PSD (processed)",
          "AWC")
names(labs) <- c("s_estimate", "s_estimate_PSD_low", "s_estimate_AWC")

data %>%
  ggplot(., aes(x = lon, y = lat, fill = slope_of_spec_exp)) +
  geom_raster() +
  coord_fixed() +
  theme_minimal() +
  facet_wrap(~Method, labeller = labeller(Method = labs)) + 
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  labs(fill = "Slope of spectral exponent", y = "Latitude", x = "Longitude") 

data %>%
  ggplot(., aes(x = Method, y = slope_of_spec_exp, fill = Method)) +
  geom_boxplot() +
  theme_light() + 
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_fill_manual(values = pal) +
  scale_x_discrete(labels = c("PSD (unprocessed)", "PSD (processed)",
                              "AWC")) +
  labs(y = "Slope of spectral exponent")  


##  now for sea surface temperature:

## make colour palette
pal = pnw_palette("Sunset",5, type = "discrete")
pal <- c(pal[1], pal[3], pal[4])

## read in old results and new results
old <- read.csv("/Volumes/SundayLab/CMIP5-GCMs_tos/01_CMCC-CESM/spec-exp_long-0-60_lat-90-30.csv")
new <- read.csv("/Volumes/SundayLab/CMIP5-GCMs_tos/01_CMCC-CESM/spec-exp_long-0-60_lat-90-30_new.csv")

## filter to only 10 year time windows 
new <- filter(new, time_window_width == "10 years")
old <- filter(old, time_window_width == "10 years")

## left join
data <- left_join(new, old) %>%
  select(lat, lon, s_spec_exp,  
         s_spec_exp_PSD_low, s_spec_exp_AWC,
         s_spec_exp_PSD_high, 
         window_start_year, window_stop_year) %>%
  gather(key = "Method",  value = "spectral_exponent", c(s_spec_exp, s_spec_exp_PSD_low, s_spec_exp_AWC)) %>%
  mutate(spectral_exponent = ifelse(Method %in% c("s_spec_exp"), 
                                    -spectral_exponent, 
                                    spectral_exponent)) %>% ## change sign of spectral exponent for PSD estimates
  mutate(Method = factor(.$Method, levels = c("s_spec_exp", "s_spec_exp_PSD_low", "s_spec_exp_AWC"),
                         ordered = T))

## compare estimates of spectral exponents from different method 
labs <- c("PSD (unprocessed)", "PSD (processed)",
          "AWC")
names(labs) <- c("s_spec_exp", "s_spec_exp_PSD_low", "s_spec_exp_AWC")

data %>%
  group_by(lat, lon, Method) %>%
  mutate(avg_spec_exp = mean(spectral_exponent)) %>%
  ungroup() %>%
  select(-spectral_exponent) %>%
  ggplot(., aes(x = lon, y = lat, fill = avg_spec_exp)) +
  geom_raster() +
  coord_fixed() +
  theme_minimal() +
  facet_wrap(~Method, labeller = labeller(Method = labs)) + 
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  labs(fill = "Avg. spectral exponent", y = "Latitude", x = "Longitude") 

data %>%
  ggplot(., aes(x = Method, y = spectral_exponent, fill = Method)) +
  geom_boxplot() +
  theme_light() + 
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_fill_manual(values = pal) +
  scale_x_discrete(labels = c("PSD (unprocessed)", "PSD (processed)",
                              "AWC")) +
  labs(y = "Spectral exponent")  

## previous method: higher spectral exponent = much steeper slopes
data %>%
  #filter(Method != "s_spec_exp") %>%
  ggplot(., aes(x = spectral_exponent, fill = Method)) + geom_histogram()

## how does change compare?
data %>%
  mutate(lat_lon_method = paste(lat, lon, Method),
         lat_lon = paste(lat, lon)) %>%
  filter(lat_lon %in% paste(unique(data$lat), unique(data$lon))[1:10]) %>%
  ggplot(., aes(x = window_start_year, y = spectral_exponent, 
                colour = lat, 
                group = lat_lon_method)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme_light() + 
  theme(panel.grid = element_blank()) +
  labs(y = "Spectral exponent", x = "Window start year")

data %>%
  group_by(Method, window_start_year) %>%
  mutate(global_avg = mean(spectral_exponent)) %>%
  ggplot(., aes(x = window_start_year, y = global_avg, 
                colour = Method, group = Method)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  theme_light() + 
  theme(panel.grid = element_blank()) +
  scale_colour_manual(values = pal, labels = c("PSD (unprocessed)", "PSD (processed)",
                                               "AWC")) +
  labs(y = "Average exponent across 30x30 degree area", x = "Window start year")

## calculate slopes
data %>%
  group_by(Method, window_start_year) %>%
  mutate(global_avg = mean(spectral_exponent)) %>%
  select(global_avg, window_start_year, Method) %>%
  group_by(Method) %>%
  unique() %>%
  do(tidy(lm(global_avg ~ window_start_year, data = .))) %>%
  filter(term == "window_start_year")



data <- left_join(new, old) %>%
  select(lat, lon, s_estimate, s_estimate_PSD_low, s_estimate_PSD_high, s_estimate_AWC, window_start_year, window_stop_year) %>%
  gather(key = "Method",  value = "slope_of_spec_exp", c(s_estimate, s_estimate_PSD_low, s_estimate_AWC)) %>%
  mutate(slope_of_spec_exp = ifelse(Method %in% c("s_estimate"), 
                                    -slope_of_spec_exp, 
                                    slope_of_spec_exp)) %>% ## change sign of slope for PSD estimates
  mutate(Method = factor(.$Method, levels = c("s_estimate", "s_estimate_PSD_low", "s_estimate_AWC"),
                         ordered = T)) %>%
  select(-window_start_year, - window_stop_year) %>%
  unique()

labs <- c("PSD (unprocessed)", "PSD (processed)",
          "AWC")
names(labs) <- c("s_estimate", "s_estimate_PSD", "s_estimate_AWC")

data %>%
  ggplot(., aes(x = lon, y = lat, fill = slope_of_spec_exp)) +
  geom_raster() +
  coord_fixed() +
  theme_minimal() +
  facet_wrap(~Method, labeller = labeller(Method = labs)) + 
  scale_fill_gradient2(low = "darkblue", high = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  labs(fill = "Slope of spectral exponent", y = "Latitude", x = "Longitude") 

data %>%
  ggplot(., aes(x = Method, y = slope_of_spec_exp, fill = Method)) +
  geom_boxplot() +
  theme_light() + 
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_fill_manual(values = pal) +
  scale_x_discrete(labels = c("PSD (unprocessed)", "PSD (processed)",
                              "AWC")) +
  labs(y = "Slope of spectral exponent")  


