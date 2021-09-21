## trying my hand at animation!!!
library(tidyverse)
library(animation)
library(raster)
library(gganimate)
library(av)

## read in my rasterStacks:
l_stack_list <- readRDS(paste(path, "l_stack_list.rds", sep = ""))
s_stack_list <- readRDS(paste(path, "s_stack_list.rds", sep = ""))

width = 5
while (width < 11) {
  ## linearly detrended 
  xyz <- purrr::map_dfr(
    as.list(l_stack_list[[width - 4]]), 
    ~setNames(as.data.frame(as(., "SpatialPixelsDataFrame")), c('z', 'x', 'y')), 
    .id = 'year'
  )
  xyz$x <- ifelse(xyz$x >= 180, xyz$x - 358, xyz$x)
  xyz$year = as.numeric(xyz$year)*width + 1871 - width
  
  countries <- map_data("world")
  
  # draw the map
  gg = xyz %>%
    ggplot(., aes(x = x, y = y, fill = z)) +
    geom_raster() + coord_fixed() + 
    theme_minimal() + 
    theme(panel.grid = element_blank()) +
    labs(x = "Longitude", y = "Latitude", fill = "Spectral exponent") +
    scale_fill_gradient2(high = "darkblue", mid = "white", low = "darkred", 
                         midpoint = 0) +
    geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
                 aes(x=long, y=lat, group = group)) 
  
  gganim <- gg + transition_states(year) + labs(title = "Start year: {closest_state}") 
  
  # animate(gganim, fps = 10, height = 800, width = 800, 
  #         start_pause = 10, end_pause = 10,
  #         renderer = av_renderer(paste("animation_L", width, "-y-window.mp4", sep = "")))
  anim_save(animation = animate(gganim, fps = 10, height = 800, width = 800, 
                                start_pause = 10, end_pause = 10), 
            paste("animation_L", width, "-y-window.gif", sep = ""))
  
  ## seasonally detrended
  xyz <- purrr::map_dfr(
    as.list(s_stack_list[[width - 4]]), 
    ~setNames(as.data.frame(as(., "SpatialPixelsDataFrame")), c('z', 'x', 'y')), 
    .id = 'year'
  )
  xyz$x <- ifelse(xyz$x >= 180, xyz$x - 358, xyz$x)
  xyz$year = as.numeric(xyz$year)*width + 1871 - width
  
  countries <- map_data("world")
  
  # draw the map
  gg = xyz %>%
    ggplot(., aes(x = x, y = y, fill = z)) +
    geom_raster() + coord_fixed() + 
    theme_minimal() + 
    theme(panel.grid = element_blank()) +
    labs(x = "Longitude", y = "Latitude", fill = "Spectral exponent") +
    scale_fill_gradient2(high = "darkblue", mid = "white", low = "darkred", 
                         midpoint = 0) +
    geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
                 aes(x=long, y=lat, group = group)) 
  
  gganim <- gg + transition_states(year) + labs(title = "Start year: {closest_state}") 
  
  # animate(gganim, fps = 10, height = 800, width = 800, 
  #         start_pause = 10, end_pause = 10,
  #         renderer = av_renderer(paste("animation_S", width, "-y-window.mp4", sep = "")))
  anim_save(animation = animate(gganim, fps = 10, height = 800, width = 800, 
                    start_pause = 10, end_pause = 10), 
            paste("animation_S", width, "-y-window.gif", sep = ""))
  
  width = width + 1
}







