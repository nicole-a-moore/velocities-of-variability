## calculating local climate velocities 

## first you will need to install the package VoCC from github:
if (!"remotes" %in% rownames(installed.packages())) {
  install.packages("remotes")
  Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
  remotes::install_github("JorGarMol/VoCC")
}

library(VoCC)
library(tidyverse)

## read in spectral exponent rasterStack data
## spectral exponent will be used as our climate variable 


## assess temporal gradient in spectral exponent using tempTrend():
## modify VoCC's tempTrend function so that slope has correct time value increments:
tempTrend <- function(r, th, t_window) {
  y <- getValues(r)
  goodtogo <- which(rowSums(is.na(y))!= ncol(y)) ## remove cells that contain all NA values
  y <- t(y[goodtogo, ])
  
  N <- apply(y, 2, function(x) sum(!is.na(x))) ## count number of non-NA observations in each cell
  ind <- which(N >= th) ## which cells have enough observations?
  y <- y[,ind] ## drop cells with less than th observations
  N <- apply(y, 2, function(x) sum(!is.na(x))) ## update count of non-NA observations
  
  x <- matrix(nrow = nlayers(r), ncol = ncol(y))
  x[] <- 1:nlayers(r)*t_window+1865
  
  ## create x values and put NA values into the x values so they correspond with y
  x1 <- y
  x1[!is.na(x1)] <- 1
  x <- x*x1
  
  # calculate the sum terms
  sx <- apply(x, 2, sum, na.rm = T)
  sy <- apply(y, 2, sum, na.rm = T)
  sxx <- apply(x, 2, function(x) sum(x^2, na.rm = T))
  syy <- apply(y, 2, function(x) sum(x^2, na.rm = T))
  xy <- x*y
  sxy <- apply(xy, 2, sum, na.rm = T)
  
  # Estimate slope coefficients and associated standard errors and p-values
  slope <- (N*sxy-(sx*sy))/(N*sxx-sx^2)
  sres <- (N*syy-sy^2-slope^2*(N*sxx-sx^2))/(N*(N-2))
  SE <- suppressWarnings(sqrt((N*sres)/(N*sxx-sx^2)))
  Test <- slope/SE
  
  p <- mapply(function(x,y) (2*pt(abs(x), df = y-2, lower.tail = FALSE)), x = Test, y = N)
  
  slpTrends <- sigTrends <- seTrends <- raster(r[[1]])
  slpTrends[goodtogo[ind]] <- slope
  
  seTrends[goodtogo[ind]] <- SE
  
  sigTrends[goodtogo[ind]] <- p
  
  output <- stack(slpTrends,seTrends,sigTrends)
  names(output) <- c("slpTrends", "seTrends", "sigTrends")
  
  return(output)
}

calculate_tempTrend <- function (stacks) {
  ## import spectral exponent rasterStacks
  l_stack_list <- stacks[[1]]
  s_stack_list <- stacks[[2]]
  
  ## run temporal trend function on each stack of varying window widths
  ltt5 <- tempTrend(l_stack_list[[1]], th = 20, 5)
  ltt6 <- tempTrend(l_stack_list[[2]], th = 20, 6)
  ltt7 <- tempTrend(l_stack_list[[3]], th = 20, 7)
  ltt8 <- tempTrend(l_stack_list[[4]], th = 20, 8)
  ltt9 <- tempTrend(l_stack_list[[5]], th = 20, 9)
  ltt10 <- tempTrend(l_stack_list[[6]], th = 20, 10)
  stt5 <- tempTrend(s_stack_list[[1]], th = 20, 5)
  stt6 <- tempTrend(s_stack_list[[2]], th = 20, 6)
  stt7 <- tempTrend(s_stack_list[[3]], th = 20, 7)
  stt8 <- tempTrend(s_stack_list[[4]], th = 20, 8)
  stt9 <- tempTrend(s_stack_list[[5]], th = 20, 9)
  stt10 <- tempTrend(s_stack_list[[6]], th = 20, 10)
  
  ## combine the output into a list
  ltt <- list(ltt5, ltt6, ltt7, ltt8, ltt9, ltt10)
  names(ltt) <- names(l_stack_list)
  stt <- list(stt5, stt6, stt7, stt8, stt9, stt10)
  names(stt) <- names(s_stack_list)
  
  tt <- list(ltt, stt)
  names(tt) <- c("linearly_detrended", "seasonally_detrended")
  
  return(tt)
}

calculate_spatGrad <- function (stacks) {
  l_stack_list <- stacks[[1]]
  s_stack_list <- stacks[[2]]
  
  ## assess spatial gradient in spectral exponent using spatGrad():
  lsg5 <- spatGrad(l_stack_list[[1]], projected = TRUE)
  lsg6 <- spatGrad(l_stack_list[[2]], projected = TRUE)
  lsg7 <- spatGrad(l_stack_list[[3]], projected = TRUE)
  lsg8 <- spatGrad(l_stack_list[[4]], projected = TRUE)
  lsg9 <- spatGrad(l_stack_list[[5]], projected = TRUE)
  lsg10 <- spatGrad(l_stack_list[[6]], projected = TRUE)
  ssg5 <- spatGrad(s_stack_list[[1]], projected = TRUE)
  ssg6 <- spatGrad(s_stack_list[[2]], projected = TRUE)
  ssg7 <- spatGrad(s_stack_list[[3]], projected = TRUE)
  ssg8 <- spatGrad(s_stack_list[[4]], projected = TRUE)
  ssg9 <- spatGrad(s_stack_list[[5]], projected = TRUE)
  ssg10 <- spatGrad(s_stack_list[[6]], projected = TRUE)
  
  ## combine the output into a list
  lsg <- list(lsg5, lsg6, lsg7, lsg8, lsg9, lsg10)
  names(lsg) <- names(s_stack_list)
  ssg <- list(ssg5, ssg6, ssg7, ssg8, ssg9, ssg10)
  names(ssg) <- names(s_stack_list)
  
  sg <- list(lsg, ssg)
  names(sg) <- c("linearly_detrended", "seasonally_detrended")
  
  return(sg)
}

calculate_VoCC <- function (tt, sg) {
  ltt <- tt[[1]]
  stt <- tt[[2]]
  lsg <- sg[[1]]
  ssg <- sg[[2]]
  
  ## calculate gradient-based climate velocities using gVoCC(): 
  lvocc5 <- gVoCC(tempTrend = ltt[[1]], spatGrad = lsg[[1]])
  lvocc6 <- gVoCC(tempTrend = ltt[[2]], spatGrad = lsg[[2]])
  lvocc7 <- gVoCC(tempTrend = ltt[[3]], spatGrad = lsg[[3]])
  lvocc8 <- gVoCC(tempTrend = ltt[[4]], spatGrad = lsg[[4]])
  lvocc9 <- gVoCC(tempTrend = ltt[[5]], spatGrad = lsg[[5]])
  lvocc10 <- gVoCC(tempTrend = stt[[6]], spatGrad = ssg[[6]])
  svocc5 <- gVoCC(tempTrend = stt[[1]], spatGrad = ssg[[1]])
  svocc6 <- gVoCC(tempTrend = stt[[2]], spatGrad = ssg[[2]])
  svocc7 <- gVoCC(tempTrend = stt[[3]], spatGrad = ssg[[3]])
  svocc8 <- gVoCC(tempTrend = stt[[4]], spatGrad = ssg[[4]])
  svocc9 <- gVoCC(tempTrend = stt[[5]], spatGrad = ssg[[5]])
  svocc10 <- gVoCC(tempTrend = stt[[6]], spatGrad = ssg[[6]])
  
  lvocc <- list(lvocc5, lvocc6, lvocc7, lvocc8, lvocc9, lvocc10)
  names(lvocc) <- names(ltt)
  svocc <- list(svocc5, svocc6, svocc7, svocc8, svocc9, svocc10)
  names(svocc) <- names(ltt)
  
  vocc <- list(lvocc, svocc)
  names(vocc) <- c("linearly_detrended", "seasonally_detrended")
  
  return(vocc)
}


## the slpTrends should be identical to the l_estimate and s_estimate column that I calculated for each cell in spec_exp
## check using 5 year window width: 
# spec_exp <- read.csv("data-processed/spec-exp.csv")

# l_estimate_raster <- filter(spec_exp, time_window_width == "5 years") %>%
#   select(long, lat, l_estimate) %>%
#   unique(.) %>%
#   rasterFromXYZ(.)

#plot(ltt5$slpTrends)
#plot(l_estimate_raster)

## YAAAAAYYYYYYY!!!! i did it!!
## this means I can take this step out later (potentially) to reduce computation time
