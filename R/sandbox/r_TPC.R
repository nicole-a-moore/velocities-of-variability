

temp_cmin <- 5
temp_cmax <- 26

# Temperature at which performance is at its maximum value.
temp_op <- (temp_cmax+temp_cmin)/3+sqrt(((temp_cmax+temp_cmin)/3)^2-
                                          (temp_cmax*temp_cmin)/3)

ro <- 2

# Temperature that occurs in the minimum time of the simulation.

temp <- seq(from = 4, to = 27, by = 0.1)

rate <- rate_TPC(temp,ro,temp_cmin,temp_cmax,temp_op)

plot(1:length(temp),rate, type="l")
