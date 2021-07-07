## trying my hand at animation!!!
library(tidyverse)
library(animation)

## read in my rasterStacks:
l_stack_list <- readRDS("data-processed/l_stack_list.rds")
s_stack_list <- readRDS("data-processed/s_stack_list.rds")

width = 5
while (width < 11) {
  l_data <- l_stack_list[[width - 4]]
  names(l_data) <- c(paste(rep(paste("L_detrended_", 
                               width, "y_window", sep = ""), 
                               nlayers(l_data)), 1870+0:nlayers(l_data)*width, sep = '_'))
  
  animation::saveGIF(raster::animate(l_data, pause = 1, n = 1),
                     movie.name = paste("animation_L", width, "-y-window.gif", sep = ""))
  
  s_data <- s_stack_list[[width - 4]]
  names(s_data) <- c(paste(rep(paste("S_detrended_", 
                                     width, "y window", sep = ""), 
                               nlayers(s_data)), 1870+0:nlayers(s_data)*width, sep = "_"))
  
  animation::saveGIF(raster::animate(s_data, pause = 1, n = 1), 
                     movie.name = paste("animation_S", width, "-year-window.gif", sep = ""))
  
  width = width + 1
}




