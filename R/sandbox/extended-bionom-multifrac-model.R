
xi = rep(1, 100000)
xi = runif(100000, -1000, 1000)
a = 0.5
b = 0.6
l = 100000

ab <- c(a,b)

nsub = 2
while (nsub <= l) {
  
  ## divide into subseries 
  xi <- split(xi, cut(seq_along(xi), nsub, labels = FALSE)) 
  ab_mult <- rep(ab, length(xi)/2)
  a_or_b <- 1:length(xi) %% 2
  
  num <- 1
  while (num <= length(xi)) {
    if(a_or_b[num] == 1) {
      xi[[num]] = xi[[num]]*a
    }
    else {
      xi[[num]] = xi[[num]]*b
    }
    num =  num + 1
  }
  
  xi <- as.vector(unlist(xi))
  
  nsub = nsub*2
}

plot(ts(xi))

spectral <- fft_calc(xi)
spectral %>%
  ggplot(aes(x = freq, y = power)) + geom_line() +
  scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm") +
  geom_vline(xintercept = 1/10)

(2*log(a+b)-log(a^2+b^2))/log(2)


q = 2
equ <- 1/q - (log(a^q + b^q))/(q*log(2)) + log(a+b)/log(2) 

a = 0.1
b = 0.1
log(a+b)/log(2)
