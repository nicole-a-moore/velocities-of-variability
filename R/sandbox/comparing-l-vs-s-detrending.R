### comparing s and l detrended data 

se_filenames = c("/Volumes/Nikki 6TB/new_spec_exp/spec-exp_long-0-60_lat-30--30.csv",
                 "/Volumes/Nikki 6TB/new_spec_exp/spec-exp_long-0-60_lat-90-30.csv",   
                 "/Volumes/Nikki 6TB/new_spec_exp/spec-exp_long-60-120_lat-30--30.csv",
                 "/Volumes/Nikki 6TB/new_spec_exp/spec-exp_long-60-120_lat-90-30.csv", 
                 "/Volumes/Nikki 6TB/new_spec_exp/spec-exp_long-120-180_lat-90-30.csv",
                 "/Volumes/Nikki 6TB/new_spec_exp/spec-exp_long-180-240_lat-90-30.csv",
                 "/Volumes/Nikki 6TB/new_spec_exp/spec-exp_long-240-300_lat-90-30.csv",
                 "/Volumes/Nikki 6TB/new_spec_exp/spec-exp_long-300-360_lat-90-30.csv")

file = 1
while (file < length(se_filenames) + 1) {
  if (file == 1) {
    spec_exp <- read.csv(se_filenames[file])
  }
  else {
    spec_exp <- rbind(spec_exp, read.csv(se_filenames[file]))
  }
  print(paste("Reading file #", file, "/", length(se_filenames), sep = ""))
  file = file + 1
}

r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
            res = 1)

## reorder time_window_width so list elements are in order of increasing time window widths:
spec_exp$time_window_width <- factor(spec_exp$time_window_width, levels = 
                                       c("5 years", "6 years", "7 years", "8 years",
                                         "9 years", "10 years"))

tas <- spec_exp %>%
  filter(time_window_width == "10 years") %>%
  select(lon, lat, l_spec_exp_PSD_low, s_spec_exp_PSD_low, l_spec_exp_PSD_high, s_spec_exp_PSD_high) %>%
  unique()

tas %>%
  ggplot(., aes(x = l_spec_exp_PSD_low, y = s_spec_exp_PSD_low, col = lat)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)

## low PSD estimate is smaller when seasonal trend is removed - makes sense, since seasonal cycle adds low frequency variation

tas %>%
  ggplot(., aes(x = l_spec_exp_PSD_high, y = s_spec_exp_PSD_high, col = lat)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)

## high PSD estimate does not change, as expected


## how do estimates of spectral change change?
tas <- spec_exp %>%
  filter(time_window_width == "10 years") %>%
  select(lon, lat, l_estimate_PSD_low, s_estimate_PSD_low, l_estimate_PSD_high, s_estimate_PSD_high) %>%
  unique()

tas %>%
  ggplot(., aes(x = l_estimate_PSD_low, y = s_estimate_PSD_low, col = lat)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)

## seems to matter for high latitudes, where seasonality is greatest 
# (ie. where annual peak in linearly-detrended data has the most power)

tas %>%
  ggplot(., aes(x = l_estimate_PSD_high, y = s_estimate_PSD_high, col = lat)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)



tas <- spec_exp %>%
  filter(time_window_width == "10 years")%>%
  select(lon, lat, l_spec_exp_PSD_low, s_spec_exp_PSD_low, l_spec_exp_PSD_high, s_spec_exp_PSD_high,
         window_start_year) %>%
  unique()

tas %>%
  group_by(window_start_year) %>%
  mutate(mean_l = mean(l_spec_exp_PSD_low), 
         mean_s = mean(s_spec_exp_PSD_low)) %>%
  select(mean_s, mean_l, window_start_year) %>%
  unique(.)%>%
  ggplot(., aes(x = window_start_year, y = mean_l), col = "blue") + geom_point() +
  geom_smooth(aes(x = window_start_year, y = mean_l), method = "lm") +
  geom_point(aes(x = window_start_year, y = mean_s),  col = "red") +
  geom_smooth(aes(x = window_start_year, y = mean_s), method = "lm")

tas %>%
  group_by(window_start_year) %>%
  mutate(mean_l = mean(l_spec_exp_PSD_high), 
         mean_s = mean(s_spec_exp_PSD_high)) %>%
  select(mean_s, mean_l, window_start_year) %>%
  unique(.)%>%
  ggplot(., aes(x = window_start_year, y = mean_l), col = "blue") + geom_point()  +
  geom_smooth(aes(x = window_start_year, y = mean_l), method = "lm") +
  geom_point(aes(x = window_start_year, y = mean_s),  col = "red") +
  geom_smooth(aes(x = window_start_year, y = mean_s), method = "lm") 

tas <- spec_exp %>%
  filter(time_window_width == "10 years")

## how does fit of linear regression differ?
ggplot(tas, aes(x = l_std.error_PSD_low)) + geom_histogram(colour = "blue") +
  geom_histogram(aes(x = s_std.error_PSD_low), colour = "red") 
  
## linear regression of linearly detrended spectral exponents has lower standard error... interesting  



## look at AWC:
tas <- spec_exp %>%
  filter(time_window_width == "10 years") %>%
  select(lon, lat, l_spec_exp_AWC, s_spec_exp_AWC) %>%
  unique()

tas %>%
  ggplot(., aes(x = l_spec_exp_AWC, y = s_spec_exp_AWC, col = lat)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)

## how do estimates of spectral change change?
tas <- spec_exp %>%
  filter(time_window_width == "10 years") %>%
  select(lon, lat, l_estimate_AWC, s_estimate_AWC) %>%
  unique()

tas %>%
  ggplot(., aes(x = l_estimate_AWC, y = s_estimate_AWC, col = lat)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)


tas <- spec_exp %>%
  filter(time_window_width == "10 years")%>%
  select(lon, lat, l_spec_exp_AWC, s_spec_exp_AWC,
         window_start_year) %>%
  unique()

tas %>%
  group_by(window_start_year) %>%
  mutate(mean_l = mean(l_spec_exp_AWC), 
         mean_s = mean(s_spec_exp_AWC)) %>%
  select(mean_s, mean_l, window_start_year) %>%
  unique(.)%>%
  ggplot(., aes(x = window_start_year, y = mean_l), col = "blue") + geom_point() +
  geom_smooth(aes(x = window_start_year, y = mean_l), method = "lm") +
  geom_point(aes(x = window_start_year, y = mean_s),  col = "red") +
  geom_smooth(aes(x = window_start_year, y = mean_s), method = "lm")

tas <- spec_exp %>%
  filter(time_window_width == "10 years")%>%
  select(lon, lat, s_spec_exp_AWC, s_spec_exp_PSD_low,
         window_start_year) %>%
  unique()


## plot a power spectrum comparison 
ts <- readRDS("data-processed/local-time-series_lat-60.5_lon-32.5.rds") 

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

window_l <- parabolic_window(ts$l_temp)
window_s <- parabolic_window(ts$s_temp)
ts_l <- ts$l_temp*window_l
ts_s <- ts$s_temp*window_s

ts_l <- bridge_detrender(ts_l)
ts_s <- bridge_detrender(ts_s)

## calculate spectral exponent over a time series window uisng periodogram
ts_window = ts_l
  l <- length(ts_window)
  
  # Fourier transform the time series window: 
  dft <- fft(ts_window)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral_l <- data.frame(freq = freq, power = amp^2)
  
  ## get estimate of spectral exponent over time series window:
  ## fit slope to low frequencies
  model_output_low <- spectral_l %>%
    filter(freq < 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  ## fit slope to high frequencies
  model_output_high <- spectral_l %>%
    filter(freq >= 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  
  ts_window = ts_s
  l <- length(ts_window)
  
  # Fourier transform the time series window: 
  dft <- fft(ts_window)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## sampling frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral_s <- data.frame(freq = freq, power = amp^2)
  
  ## get estimate of spectral exponent over time series window:
  ## fit slope to low frequencies
  model_output_low <- spectral_s %>%
    filter(freq < 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  ## fit slope to high frequencies
  model_output_high <- spectral_s %>%
    filter(freq >= 1/8*max(spectral$freq)) %>%
    lm(., formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  

  # ## plot spectrum:
  
  spectral_l$type = "linearly detrended"
  spectral_s$type = "seasonally detrended"
  
  spectral_l %>%
    rbind(., spectral_s) %>%
    ggplot(aes(x = freq, y = power, group = type, colour = type)) + 
    geom_vline(xintercept = 1/365, colour = "black") +
    geom_line() +
    scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")





