## spatial autocorrelation of spectral exponent 
library(gstat)
library(geosphere)
library(tidyverse)
library(raster)
select <- dplyr::select

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Calculate spatial autocorrelation of spectral exponent  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## extract 10 year time window 
l_stack_tas <- readRDS("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/l_stack_list.rds")[[6]] 
s_stack_tas <- readRDS("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/s_stack_list.rds")[[6]] 
extent(l_stack_tas) <- c(-179, 179, -89, 89)
extent(s_stack_tas) <- c(-179, 179, -89, 89)

l_stack_tos <- readRDS("data-raw/01_CMCC-CESM/l_stack_list_tos.rds")[[6]] 
s_stack_tos <- readRDS("data-raw/01_CMCC-CESM/s_stack_list_tos.rds")[[6]] 

## make sure that extents match by cropping
l_stack_tas <- crop(l_stack_tas, l_stack_tos)
s_stack_tas <- crop(s_stack_tas, s_stack_tos)
l_stack_tos <- crop(l_stack_tos, l_stack_tas)
s_stack_tos <- crop(s_stack_tos, s_stack_tas)

## combine tos and tas data into one raster stack
isna <- is.na(l_stack_tos[[1]])
l_stack_tos[isna == 1] <- l_stack_tas[isna == 1]
plot(l_stack_tos[[1]])

isna <- is.na(s_stack_tos[[1]])
s_stack_tos[isna == 1] <- s_stack_tas[isna == 1]
plot(s_stack_tos[[1]])

## rename
s_mosaic <- s_stack_tos
l_mosaic <- l_stack_tos

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Plot a semi-variogram        #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## convert rasterstacks to points 
s_points <- data.frame(rasterToPoints(s_mosaic))
l_points <- data.frame(rasterToPoints(l_mosaic))

## calculate variogram 
var <- variogram(window_1 ~ 1, locations = ~ x + y, cutoff = 300,
                 data = s_points)
var
plot(var, col = "black")

## fit relationship
vfit <- fit.variogram(var, vgm("Sph"))
vfit
plot(var, vfit, col = "black")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  How does the spatial range of spectral exponent change over time?  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## for each time window:
col = 3
while (col < ncol(l_points)+1) {
  
  ## calculate variogram for window 
  var <- variogram(l_points[,col] ~ 1, locations = ~ x + y, cutoff = 300,
                   data = l_points)
  var
  plot(var, col = "black")
  
  ## fit relationship
  vfit <- fit.variogram(var, vgm("Sph"))
  vfit
  plot(var, vfit, col = "black")
  
  ## store model parameters:
  if (col == 3) {
    models <- data.frame(window = col-2, 
                         nug = vfit$psill[1], 
                         psill = vfit$psill[2],
                         range = vfit$range[2])
  }
  else {
    models <- rbind(models,
                    data.frame(window = col-2, 
                               nug = vfit$psill[1],
                               psill = vfit$psill[2],
                               range = vfit$range[2]))
  }
  
  col = col + 1
}

## decrease in autocorrelation at window ~9? 
ggplot(models, aes(x = window, y = range)) + geom_point()



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  How does the spatial range of spectral exponent compare to spatial range of temperature?  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## calculate mean detrended temperature in each 10y time series window
## this will take time:
#   - read in spatial chunks, group by time windows and calculate mean, save
#   - combine all spatial chunks 
#   - do the same for tos / tos 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Land vs. ocean: how does the spatial range of spectral exponent change over time?         #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~##
#####  Land  #####
##~~~~~~~~~~~~~~##


##~~~~~~~~~~~~~~~##
#####  Ocean  #####
##~~~~~~~~~~~~~~~##








##### garbage #####
s_points$id_let <- paste(1:nrow(s_points), "a", sep = "")
s_points$id_num <- 1:nrow(s_points)

## get all combinations of grid cell coordinates
grid <- expand.grid(s_points$id_let, s_points$id_num)
colnames(grid) <- c("id_let", "id_num")
points_let <- select(s_points, x, y, id_let)
points_num <- select(s_points, x, y, id_num)
all_combs <- left_join(grid, points_let)
names(all_combs)[3:4] <- c("lon_1", "lat_1")
all_combs <- left_join(all_combs, points_num)
names(all_combs)[5:6] <- c("lon_2", "lat_2")
all_combs <- select(all_combs, -id_let, -id_num)

## compute pairwise distances between points 
distm(s_points[1,1:2], s_points[1,1:2], fun = distHaversine)