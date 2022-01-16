library(longmemo)
library(wsyn)
library(biwavelet)

i = 1 
while (i < 100) {
  ts <- simARMA0(n = 3650, H = 0.5)
  
  clean_ts <- cleandat(ts, 1:3650, clev = 2)
  wave_wsyn <- wsyn::wt(t.series = clean_ts$cdat, times = 1:3650)
  
  df <- abs(wave_wsyn$values) ## get magnitude of wavelet coeffs
  avg <- colMeans(df, na.rm = T) ## get average
  
  avg <- as.numeric(avg)
  
  data <- data.frame(avg_wavelet_coeff = avg, period = wave_wsyn$timescales)
  
  # plot average coefficients versus period on a log–log plot
  # plot <- ggplot(data = data, aes(x = period, y = avg_wavelet_coeff)) +
  #   geom_point() +
  #   scale_x_log10() + scale_y_log10() +
  #   theme_light() +
  #   labs(x = "Period", y = "Average wavelet coefficient")

  ## c. calculate beta
  ## get slope
  lm <- lm(log(avg_wavelet_coeff) ~ log(period), 
           data = data)
  
  ## slope equal to H + 1/2
  H_calc = as.numeric(lm$coefficients[2]) - 1/2
  b_wsyn = 2*H_calc + 1
  
  
  wave_bi <- biwavelet::wt(data.frame(time = 1:3650, val = ts), do.sig = F)
  
  df <- sqrt(wave_bi$power) ## get magnitude of wavelet coeffs
  avg <- rowMeans(df, na.rm = T) ## get average
  avg <- as.numeric(avg)
  
  data_2 <- data.frame(avg_wavelet_coeff = avg, period = wave_bi$period)
  
  ## plot average coefficients versus period on a log–log plot
  # plot + geom_point(data = data_2, aes(x = period, y = avg_wavelet_coeff), colour= "red") +
  #   theme_light() +
  #   labs(x = "Period", y = "Average wavelet coefficient")
  
  ## c. calculate beta
  ## get slope
  lm <-  lm(log(avg_wavelet_coeff) ~ log(period), data = data_2) 
  
  ## slope equal to H + 1/2
  H_calc = as.numeric(lm$coefficients[2]) - 1/2
  b_bi = 2*H_calc + 1
  
  if (i ==1) {
    bs<- data.frame(method = c('b_wsyn', "b_bi"),
                    b = c(b_wsyn, b_bi))
  }
  else {
    bs <-rbind(bs, data.frame(method = c('b_wsyn', "b_bi"),
                              b = c(b_wsyn, b_bi)))
               
  }
  
  i = i + 1
}

ggplot(bs, aes(x = method, y = b)) + geom_boxplot() +
  stat_summary(geom = 'point', fun = mean)


