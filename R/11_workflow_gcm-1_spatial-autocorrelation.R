## spatial autocorrelation of spectral exponent 
library(gstat)
library(geosphere)
library(tidyverse)
library(raster)
select <- dplyr::select
library(PNWColors)


## create colour palette:=
pal = pnw_palette("Sunset",5, type = "discrete")
pal_realm <- c(pal[1], pal[3])
pal_lat <- c(pal[2], pal[4])

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####   Simulate varying degrees of spatial autocorrelation   #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
library(gstat)
list_of_grids <- list()
ranges <- c(1, 10)
i=1
while (i < length(ranges)+1) {
  gridDim <- 60
  xy <- expand.grid(y=1:gridDim, x=1:gridDim)
  
  # Variogram model, with defined sill and range
  varioMod <- vgm(psill=0.05, range=ranges[i], model='Exp')
  
  # Set up an additional variable from simple kriging
  zDummy <- gstat(formula=z~1, locations = ~x+y, dummy=TRUE, 
                  beta=1, model=varioMod, nmax=20)
  
  # Generate 2 randomly autocorrelated predictor data fields
  set.seed(3)
  xyz <- predict(zDummy, newdata=xy, nsim=2)
  
  # Generate an autocorrelated response variable:
  # Deterministic on predictors, with normal error
  e <- rnorm(nrow(xyz), sd=.5)
  xyz$resp <- .7*xyz$sim1 + 1.4*xyz$sim2 + e
  
  raster <- rasterFromXYZ(xyz)
  plot(raster[[1]])
  
  list_of_grids <- append(list_of_grids, raster[[1]])
  
  i=i+1
}
plot(list_of_grids[[2]])

## convert to data frame
stack <- stack(list_of_grids[1:2])
names(stack) <- c("Range_1", "Range_10")
stack <- as.data.frame(rasterToPoints(stack))
stack <- gather(stack, key = "range", value = "fake_variable", c(3:4))
stack$range <- factor(stack$range, 
                         levels =c("Range_1", "Range_10"),
                         ordered = T)

ggplot(stack, aes(x = x, y = y, fill = fake_variable)) + coord_fixed() +
  geom_raster() +
  facet_grid(rows = ~range) +
  theme_void() +
  labs(fill = "") +
  scale_fill_gradient(high = pal[4], low = pal[1])
  

## now calculate and plot empirical variogram for each
stack <- stack(list_of_grids[1:2])
names(stack) <- c("Range_1", "Range_10")
stack <- as.data.frame(rasterToPoints(stack))
i=1
while (i < length(ranges) + 1) {
  pts <- stack[,c(1:2,i+2)]
  colnames(pts)[3] <- "fake_variable"
  
  var <- variogram(fake_variable ~ 1, locations = ~ x + y, 
                   data = pts)
  var
  plot(var, col = "black")

  vfit <- fit.variogram(var, vgm("Exp"))
  vfit
  plot(var, vfit, col = "black")
  
  i = i+1
}


  

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Calculate spatial autocorrelation of spectral exponent  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
mosaic_specexp <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        Plot a semi-variogram        #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## convert rasterstacks to points 
points <- data.frame(rasterToPoints(mosaic_specexp))

## calculate variogram 
var <- variogram(window_1 ~ 1, locations = ~ x + y, cutoff = 300,
                 data = points)
var
plot(var, col = "black")

## fit relationship
vfit <- fit.variogram(var, vgm("Sph"))
vfit
plot(var, vfit, col = "black")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  3a. How does the spatial range of spectral exponent change over time?  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## for each time window:
col = 3
while (col < ncol(points)+1) {
  
  ## calculate variogram for window 
  var <- variogram(points[,col] ~ 1, locations = ~ x + y, cutoff = 300,
                   data = points)
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

models$year <- seq(1871, 2081, by = 10)

ggplot(models, aes(x = year, y = range)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spatial range of spectral exponent (km)") +
  geom_smooth(colour = "black")



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  3b. Does change in spatial range of spectral exponent differ between ocean/land?  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
mosaic_specexp <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent.rds")
land_mask <- readRDS("data-processed/01_tos-tas-land-mask.rds")

## split data into land and ocean
land <- mask(mosaic_specexp, land_mask)
plot(land[[1]])
ocean <- mask(mosaic_specexp, land_mask, inverse = T)
plot(ocean[[1]])

#### LAND
## convert rasterstacks to points 
points <- data.frame(rasterToPoints(land))
## for each time window:
col = 3
while (col < ncol(points)+1) {
  
  ## calculate variogram for window 
  var <- variogram(points[,col] ~ 1, locations = ~ x + y, cutoff = 300,
                   data = points)
  var
  plot(var, col = "black")
  
  ## fit relationship
  vfit <- fit.variogram(var, vgm("Sph"))
  vfit
  plot(var, vfit, col = "black")
  
  ## store model parameters:
  if (col == 3) {
    models_land <- data.frame(window = col-2, 
                         nug = vfit$psill[1], 
                         psill = vfit$psill[2],
                         range = vfit$range[2])
  }
  else {
    models_land <- rbind(models_land,
                    data.frame(window = col-2, 
                               nug = vfit$psill[1],
                               psill = vfit$psill[2],
                               range = vfit$range[2]))
  }
  
  col = col + 1
}

models_land$year <- seq(1871, 2081, by = 10)

ggplot(models_land, aes(x = year, y = range)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spatial range of spectral exponent (km)") +
  geom_smooth(colour = pal_realm[1], fill = pal_realm[1])

#### OCEAN
## convert rasterstacks to points
points <- data.frame(rasterToPoints(ocean))
## for each time window:
col = 3
while (col < ncol(points)+1) {
  
  ## calculate variogram for window 
  var <- variogram(points[,col] ~ 1, locations = ~ x + y, cutoff = 300,
                   data = points)
  var
  plot(var, col = "black")
  
  ## fit relationship
  vfit <- fit.variogram(var, vgm("Sph"))
  vfit
  plot(var, vfit, col = "black")
  
  ## store model parameters:
  if (col == 3) {
    models_ocean <- data.frame(window = col-2, 
                              nug = vfit$psill[1], 
                              psill = vfit$psill[2],
                              range = vfit$range[2])
  }
  else {
    models_ocean <- rbind(models_ocean,
                         data.frame(window = col-2, 
                                    nug = vfit$psill[1],
                                    psill = vfit$psill[2],
                                    range = vfit$range[2]))
  }
  
  col = col + 1
}

models_ocean$year <- seq(1871, 2081, by = 10)

ggplot(models_ocean, aes(x = year, y = range)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spatial range of spectral exponent (km)") +
  geom_smooth(colour = pal_realm[2], fill = pal_realm[2])

models_ocean$realm = "ocean"
models_land$realm = "land"

## plot them together
models_ocean %>%
  rbind(., models_land) %>%
  ggplot(., aes(x = year, y = range, colour = realm, fill = realm)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spatial range of spectral exponent (km)") +
  geom_smooth() +
  scale_color_manual(values = pal_realm) +
  scale_fill_manual(values = pal_realm) +
  theme(legend.position = "none")

## boxplot comparing spatial range
models_ocean %>%
  rbind(., models_land) %>%
  ggplot(., aes(x = realm, y = range, fill = realm)) + 
  geom_boxplot() +
  scale_fill_manual(values = pal_realm) + 
  guides(fill = "none") +
  labs(x = "", y = "Spatial range of spectral exponent across all years (km)") +
  theme_light() +
  scale_x_discrete(labels = c("Land", "Ocean"))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  4a. How autocorrelated is the change in spectral exponent? #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specchange
mosaic_specchange <- readRDS("data-processed/01_tos-tas-mosaic_spectral-change.rds")

## convert rasterstacks to points 
points <- data.frame(rasterToPoints(land[[2]]))

## calculate variogram for window 
var <- variogram(points[,3] ~ 1, locations = ~ x + y, cutoff = 300,
                 data = points)
var
plot(var, col = "black")

## fit relationship
vfit <- fit.variogram(var, vgm("Sph"))
vfit
plot(var, vfit, col = "black")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  4b. How autocorrelated is the change in spectral exponent on land vs. ocean? #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specchange
mosaic_specchange <- readRDS("data-processed/01_tos-tas-mosaic_spectral-change.rds")

## split data into land and ocean
land <- mask(mosaic_specchange, land_mask)
plot(land[[1]])
ocean <- mask(mosaic_specchange, land_mask, inverse = T)
plot(ocean[[1]])

#### LAND
## convert rasterstacks to points 
points <- data.frame(rasterToPoints(land[[2]]))

## calculate variogram for window 
var <- variogram(points[,3] ~ 1, locations = ~ x + y,
                   data = points)
var
plot(var, col = "black")
  
## fit relationship
vfit <- fit.variogram(var, vgm("Sph"))
vfit
plot(var, vfit, col = "black")
  
 
#### OCEAN
## convert rasterstacks to points
points <- data.frame(rasterToPoints(ocean[[2]]))

## calculate variogram for window 
var <- variogram(points[,3] ~ 1, locations = ~ x + y, cutoff = 300,
                 data = points)
var
plot(var, col = "black")

## fit relationship
vfit <- fit.variogram(var, vgm("Sph"))
vfit
plot(var, vfit, col = "black")







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



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 3d.Does change in spatial range of spectral exponent differ between tropical/temperate ocean/land?  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specchange
mosaic_specexp <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent.rds")

## split data into tropics and temperate latitudes
tropical <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)

temperate <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y > 23.5 | y < -23.5) %>%
  rasterFromXYZ(.)

## split tropical and temperate into land and ocean
tropical_land <- crop(land_mask, tropical)
tropical_land <- mask(tropical, tropical_land)
tropical_ocean <- mask(tropical, tropical_land, inverse = T)
plot(tropical_land[[1]])
plot(tropical_ocean[[1]])

temperate_land <- crop(land_mask, temperate)
temperate_land <- mask(temperate, temperate_land)
temperate_ocean <- mask(temperate, temperate_land, inverse = T)
plot(temperate_land[[1]])
plot(temperate_ocean[[1]])


#### TROPICAL LAND
## convert rasterstacks to points 
points <- data.frame(rasterToPoints(temperate_land))
## for each time window:
col = 3
while (col < ncol(points)+1) {
  
  ## calculate variogram for window 
  var <- variogram(points[,col] ~ 1, locations = ~ x + y, cutoff = 300,
                   data = points)
  var
  plot(var, col = "black")
  
  ## fit relationship
  vfit <- fit.variogram(var, vgm("Sph"))
  vfit
  plot(var, vfit, col = "black")
  
  ## store model parameters:
  if (col == 3) {
    models_trop_land <- data.frame(window = col-2, 
                                   nug = vfit$psill[1], 
                                   psill = vfit$psill[2],
                                   range = vfit$range[2])
  }
  else {
    models_trop_land <- rbind(models_trop_land,
                              data.frame(window = col-2, 
                                         nug = vfit$psill[1],
                                         psill = vfit$psill[2],
                                         range = vfit$range[2]))
  }
  
  col = col + 1
}

models_trop_land$year <- seq(1871, 2081, by = 10)

ggplot(models_trop_land, aes(x = year, y = range)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spatial range of spectral exponent (km)") +
  geom_smooth(colour = pal_lat[2], fill = pal_lat[2])



#### temperate
## convert rasterstacks to points
points <- data.frame(rasterToPoints(temperate))
## for each time window:
col = 3
while (col < ncol(points)+1) {
  
  ## calculate variogram for window 
  var <- variogram(points[,col] ~ 1, locations = ~ x + y, cutoff = 300,
                   data = points)
  var
  plot(var, col = "black")
  
  ## fit relationship
  vfit <- fit.variogram(var, vgm("Sph"))
  vfit
  plot(var, vfit, col = "black")
  
  ## store model parameters:
  if (col == 3) {
    models_temp <- data.frame(window = col-2, 
                              nug = vfit$psill[1], 
                              psill = vfit$psill[2],
                              range = vfit$range[2])
  }
  else {
    models_temp <- rbind(models_temp,
                         data.frame(window = col-2, 
                                    nug = vfit$psill[1],
                                    psill = vfit$psill[2],
                                    range = vfit$range[2]))
  }
  
  col = col + 1
}

models_temp$year <- seq(1871, 2081, by = 10)

ggplot(models_temp, aes(x = year, y = range)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spatial range of spectral exponent (km)") +
  geom_smooth(colour = pal_lat[1], fill = pal_lat[1])


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 3c.Does change in spatial range of spectral exponent differ between tropical/temperate?  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
mosaic_specexp <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent.rds")

## split data into tropics and temperate latitudes
tropical <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)

temperate <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y > 23.5 | y < -23.5) %>%
  rasterFromXYZ(.)

#### TROPICAL
## convert rasterstacks to points 
points <- data.frame(rasterToPoints(tropical))
## for each time window:
col = 3
while (col < ncol(points)+1) {
  
  ## calculate variogram for window 
  var <- variogram(points[,col] ~ 1, locations = ~ x + y, cutoff = 300,
                   data = points)
  var
  plot(var, col = "black")
  
  ## fit relationship
  vfit <- fit.variogram(var, vgm("Sph"))
  vfit
  plot(var, vfit, col = "black")
  
  ## store model parameters:
  if (col == 3) {
    models_trop <- data.frame(window = col-2, 
                              nug = vfit$psill[1], 
                              psill = vfit$psill[2],
                              range = vfit$range[2])
  }
  else {
    models_trop <- rbind(models_trop,
                         data.frame(window = col-2, 
                                    nug = vfit$psill[1],
                                    psill = vfit$psill[2],
                                    range = vfit$range[2]))
  }
  
  col = col + 1
}

models_trop$year <- seq(1871, 2081, by = 10)

ggplot(models_trop, aes(x = year, y = range)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spatial range of spectral exponent (km)") +
  geom_smooth(colour = pal_lat[2], fill = pal_lat[2])

#### temperate
## convert rasterstacks to points
points <- data.frame(rasterToPoints(temperate))
## for each time window:
col = 3
while (col < ncol(points)+1) {
  
  ## calculate variogram for window 
  var <- variogram(points[,col] ~ 1, locations = ~ x + y, cutoff = 300,
                   data = points)
  var
  plot(var, col = "black")
  
  ## fit relationship
  vfit <- fit.variogram(var, vgm("Sph"))
  vfit
  plot(var, vfit, col = "black")
  
  ## store model parameters:
  if (col == 3) {
    models_temp <- data.frame(window = col-2, 
                              nug = vfit$psill[1], 
                              psill = vfit$psill[2],
                              range = vfit$range[2])
  }
  else {
    models_temp <- rbind(models_temp,
                         data.frame(window = col-2, 
                                    nug = vfit$psill[1],
                                    psill = vfit$psill[2],
                                    range = vfit$range[2]))
  }
  
  col = col + 1
}

models_temp$year <- seq(1871, 2081, by = 10)

ggplot(models_temp, aes(x = year, y = range)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Spatial range of spectral exponent (km)") +
  geom_smooth(colour = pal_lat[1], fill = pal_lat[1])





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






