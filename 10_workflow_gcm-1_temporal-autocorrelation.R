## whole workflow for first GCM  !!!
library(tidyverse)
library(raster)
library(PNWColors)
select <- dplyr::select

## create colour palette:
pal = pnw_palette("Sunset",5, type = "discrete")
pal_realm <- c(pal[1], pal[3])
pal_lat <- c(pal[2], pal[4])

## for now: just do for seasonally detrended data, not linearly detrended
## create data that has tas on land and tos in the ocean
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       Create mosaics of tos and tas       #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       spectral exponent      #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## extract 10 year time window 
s_stack_tas <- readRDS("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/s_stack_list.rds")[[6]] 
extent(s_stack_tas) <- c(-179, 179, -89, 89)
s_stack_tos <- readRDS("data-raw/01_CMCC-CESM/s_stack_list_tos.rds")[[6]] 

s_stack_tas <- crop(s_stack_tas, s_stack_tos)
s_stack_tos <- crop(s_stack_tos, s_stack_tas)

## mosaic rasters together 
isna <- is.na(s_stack_tos[[1]])
s_stack_tos[isna == 1] <- s_stack_tas[isna == 1]

mosaic_specexp <- s_stack_tos
plot(mosaic_specexp[[1]])

## keep track of land and ocean by making a mask:
land_mask <- s_stack_tas
land_mask[land_mask != mosaic_specexp] <- NA
land_mask[land_mask == mosaic_specexp] <- 1
land_mask <- land_mask[[1]]

## save:
saveRDS(mosaic_specexp, "data-processed/01_tos-tas-mosaic_spectral-exponent.rds")
saveRDS(land_mask, "data-processed/01_tos-tas-land-mask.rds")


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        spectral change       #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
### get tas data
path = "/Volumes/SundayLab/CMIP5-GCMs/" 

## create vector of file folders to put data into:
gcm_models <- c("01_CMCC-CESM", "02_CMCC-CM", '03_CMCC-CMS', '04_MPI-ESM-LR', '05_MPI-ESM-MR',
                "06_GFDL-ESM2G", '07_GFDL-CM3', '08_GFDL-ESM2M', '09_HadGEM2-CC', '10_HadGEM2-ES',
                "11_HadGEM2-AO", '12_IPSL-CM5A-LR', '13_IPSL-CM5B-LR', '14_MIROC5', '15_MIROC5-ESM-CHEM',
                '16_MIROC5-ESM', "17_inmcm4", '18_CNRM-CM5', "19_MRI-CGCM3", '20_MRI-ESM1',
                '21_IPSL-CM5A-MR')

folders <- paste(path, gcm_models, "/", sep = "")

path = folders[1]

## calculate average change in spectral exponent for each time series, across all sliding window widths 
se_filenames <- readRDS(paste(path, "se_filenames.rds",  sep = ""))

## combine all spectral exponent csvs into one big dataframe
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
  select(lon, lat, l_estimate, s_estimate, l_p.value, s_p.value) %>%
  unique()

### get tos data
path = "data-raw/" 

## create vector of file folders to put data into:
gcm_models <- c("01_CMCC-CESM", "02_CMCC-CM", '03_CMCC-CMS', '04_MPI-ESM-LR', '05_MPI-ESM-MR',
                "06_GFDL-ESM2G", '07_GFDL-CM3', '08_GFDL-ESM2M', '09_HadGEM2-CC', '10_HadGEM2-ES',
                "11_HadGEM2-AO", '12_IPSL-CM5A-LR', '13_IPSL-CM5B-LR', '14_MIROC5', '15_MIROC5-ESM-CHEM',
                '16_MIROC5-ESM', "17_inmcm4", '18_CNRM-CM5', "19_MRI-CGCM3", '20_MRI-ESM1',
                '21_IPSL-CM5A-MR')

folders <- paste(path, gcm_models, "/", sep = "")

path = folders[1]

## calculate average change in spectral exponent for each time series, across all sliding window widths 
se_filenames <- readRDS(paste(path, "se_filenames_tos.rds",  sep = ""))

## combine all spectral exponent csvs into one big dataframe
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

tos <- spec_exp %>%
  filter(time_window_width == "10 years") %>%
  select(lon, lat, l_estimate, s_estimate, l_p.value, s_p.value) %>%
  unique()


## mosaic together tas and tos 
raster_tas <- rasterFromXYZ(tas)
raster_tos <- rasterFromXYZ(tos)
extent(raster_tas) <- c(-179, 179, -89, 89)

## crop
raster_tas <- crop(raster_tas, raster_tos)
raster_tos <- crop(raster_tos, raster_tas)

isna <- is.na(raster_tos[[1]])
raster_tos[isna == 1] <- raster_tas[isna == 1]

mosaic_specchange <- raster_tos
plot(mosaic_specchange$s_estimate)
plot(mosaic_specchange$s_estimate)

## save:
saveRDS(mosaic_specchange, "data-processed/01_tos-tas-mosaic_spectral-change.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 1a. Analyze spectral exponent on land vs. ocean #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specexp
mosaic_specexp <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent.rds")

## split data into land and ocean
land <- mask(mosaic_specexp, land_mask)
plot(land[[1]])
ocean <- mask(mosaic_specexp, land_mask, inverse = T)
plot(ocean[[1]])

## summary stats and figures:
## mean
mean(values(land), na.rm = T)
mean(values(ocean), na.rm = T)

## standard deviation
sd(values(land), na.rm = T)
sd(values(ocean), na.rm = T)

## histogram:
mean_land <- calc(land, mean)
mean_ocean <- calc(ocean, mean)
df <- data.frame(land = c(values(mean_land)), ocean = c(values(mean_ocean)))
df <- gather(df, key = "realm", value = "mean_spec_exp", land, ocean)
df <- filter(df, !is.na(mean_spec_exp))

df %>%
  ggplot(., aes(x = mean_spec_exp, fill = realm)) + 
  geom_histogram(position = position_dodge(), binwidth = 0.05) +
  theme_light() +
  labs(x = "Mean local spectral exponent",
       y = "Frequency", fill = "Realm") +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none")

legend_rast <- land_mask
legend_rast[is.na(land_mask)] <- 2

map_legend <- data.frame(rasterToPoints(legend_rast)) %>%
  mutate(window_1 = as.factor(window_1)) %>%
  ggplot(., aes(x = x, y = y, fill = window_1)) +
  geom_raster() +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none")

## boxplot:
df %>%
  ggplot(., aes(y = mean_spec_exp, x = realm, fill = realm)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none") +
  scale_x_discrete(labels = c("Land", "Ocean"))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 1b. Analyze spectral exponent across tropical vs. temperate latitudes #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specexp
mosaic_specexp <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent.rds")

## split data into tropics and temperate latitudes
tropical <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)

temperate <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y > 23.5 | y < -23.5) %>%
  rasterFromXYZ(.)
  
## summary stats and figures:
## mean
mean(values(tropical), na.rm = T)
mean(values(temperate), na.rm = T)

## standard deviation
sd(values(tropical), na.rm = T)
sd(values(temperate), na.rm = T)

## histogram:
mean_trop <- calc(tropical, mean)
mean_temp <- calc(temperate, mean)
df <- data.frame(latitudinal_region = rep("tropical", length(values(mean_trop))),
                 mean_spec_exp = c(values(mean_trop)))
df <- rbind(df, data.frame(latitudinal_region = rep("temperate", length(values(mean_temp))),
                 mean_spec_exp = c(values(mean_temp))))
df <- filter(df, !is.na(mean_spec_exp))

df %>%
  ggplot(., aes(x = mean_spec_exp, fill = latitudinal_region)) + 
  geom_histogram(position = position_dodge(), binwidth = 0.05) +
  theme_light() +
  labs(x = "Mean local spectral exponent",
       y = "Frequency", fill = "Latitudinal region") +
  scale_fill_manual(values = pal_lat) +
  guides(fill = "none")

legend_rast <- land_mask
legend_rast[is.na(temperate[[1]])] <- 2
legend_rast[!is.na(temperate[[1]])] <- 1

map_legend <- data.frame(rasterToPoints(legend_rast)) %>%
  mutate(window_1 = as.factor(window_1)) %>%
  ggplot(., aes(x = x, y = y, fill = window_1)) +
  geom_raster() +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = pal_lat) +
  guides(fill = "none") +
  geom_polygon(data = countries, col="transparent", size = 0.01, fill = "transparent", alpha = 0.2,
               aes(x=long, y=lat, group = group)) 

## boxplot:
df %>%
  ggplot(., aes(y = mean_spec_exp, x = latitudinal_region, fill = latitudinal_region)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  scale_fill_manual(values = pal_lat) +
  guides(fill = "none") +
  scale_x_discrete(labels = c("Temperate", "Tropical"))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####           1c. Tropical vs. temperate latitudes x land vs. ocean                 #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specexp
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

## calculate means


## summary stats and figures:
## mean
mean(values(tropical_land), na.rm = T)
mean(values(temperate_land), na.rm = T)
mean(values(tropical_ocean), na.rm = T)
mean(values(temperate_ocean), na.rm = T)

## standard deviation
sd(values(tropical_land), na.rm = T)
sd(values(temperate_land), na.rm = T)
sd(values(tropical_ocean), na.rm = T)
sd(values(temperate_ocean), na.rm = T)

## histogram - tropical vs. temperate 
## calculate mean
tropical_land <- calc(tropical_land, mean) 
tropical_ocean <- calc(tropical_ocean, mean) 
temperate_land <- calc(temperate_land, mean) 
temperate_ocean <- calc(temperate_ocean, mean) 

df1 <- data.frame(Land = c(values(tropical_land)), 
                  Ocean = c(values(tropical_ocean)))
df1 <- gather(df1, key = "realm", value = "mean_spec_exp", Land, Ocean)
df1$trop_or_temp = "Tropical"
df2 <- data.frame(Land = c(values(temperate_land)),
                  Ocean = c(values(temperate_ocean)))
df2 <- gather(df2, key = "realm", value = "mean_spec_exp", Land, Ocean)
df2$trop_or_temp = "Temperate"
df <- rbind(df1, df2) %>% 
  filter(!is.na(mean_spec_exp))
df$lat_and_realm = paste(df$trop_or_temp, df$realm, sep = "_")

df %>%
  ggplot(., aes(x = mean_spec_exp, fill = lat_and_realm)) + 
  geom_histogram(position = position_dodge()) +
  theme_light() +
  labs(x = "Mean local spectral exponent",
       y = "Frequency",
       fill = "Latitude and realm") +
  facet_wrap(~realm) +
  scale_fill_discrete(labels = c("Temperate land", "Temperate ocean", "Tropical land", "Tropical ocean"))

## map it!
trop_land_df <- data.frame(rasterToPoints(tropical_land))
trop_land_df$lat_and_realm = "Tropical_land"
temp_land_df <- data.frame(rasterToPoints(temperate_land))
temp_land_df$lat_and_realm = "Temperate_land"
trop_ocean_df <- data.frame(rasterToPoints(tropical_ocean))
trop_ocean_df$lat_and_realm = "Tropical_ocean"
temp_ocean_df <- data.frame(rasterToPoints(temperate_ocean))
temp_ocean_df$lat_and_realm = "Temperate_ocean"
map_df <- rbind(temp_land_df,trop_ocean_df, trop_land_df, temp_ocean_df)

map_df %>%
  ggplot(., aes(x = x, y = y, fill = lat_and_realm)) + 
  geom_raster() +
  coord_fixed() +
  theme_void() +
  labs(fill = "Mean local\nspectral\nexponent") +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))  

## boxplot:
df %>%
  ggplot(., aes(y = mean_spec_exp, x = lat_and_realm, fill = lat_and_realm)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  guides(fill = "none") +
  scale_x_discrete(labels = 
                     c("Temperate land", "Temperate ocean", "Tropical land", "Tropical ocean"))



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 2a. Analyze change in spectral exponent on land vs. ocean #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specchange
mosaic_specchange <- readRDS("data-processed/01_tos-tas-mosaic_spectral-change.rds")

## split data into land and ocean
land <- mask(mosaic_specchange, land_mask)
plot(land[[1]])
ocean <- mask(mosaic_specchange, land_mask, inverse = T)
plot(ocean[[1]])

## summary stats and figures:
## mean
mean(values(land[[2]]), na.rm = T)
mean(values(ocean[[2]]), na.rm = T)

## standard deviation
sd(values(land[[2]]), na.rm = T)
sd(values(ocean[[2]]), na.rm = T)

## histogram:
df <- data.frame(Land = c(values(land[[2]])), Ocean = c(values(ocean[[2]])))
df <- gather(df, key = "realm", value = "slope_spec_exp", Land, Ocean)
df <- filter(df, !is.na(slope_spec_exp))

df %>%
  ggplot(., aes(x = slope_spec_exp, group = realm, fill = ..x..)) + 
  geom_histogram(position = position_dodge()) +
  theme_light() +
  labs(x = "Slope of spectral exponent",
       y = "Frequency",
       fill = "Slope of\nspectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  facet_wrap(~realm)

## map it!
land_df <- data.frame(rasterToPoints(land))
land_df$realm = "Land"
ocean_df <- data.frame(rasterToPoints(ocean))
ocean_df$realm = "Ocean"
map_df <- rbind(ocean_df, land_df)

map_df %>%
  ggplot(., aes(x = x, y = y, fill = s_estimate)) + 
  geom_raster() +
  coord_fixed() +
  theme_void() +
  labs(fill = "Slope of\nspectral\nexponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  facet_wrap(~realm) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))

## boxplot:
df %>%
  ggplot(., aes(y = slope_spec_exp, x = realm, fill = realm)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 2b. Analyze change in spectral exponent across tropical vs. temperate latitudes #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specchange
mosaic_specchange <- readRDS("data-processed/01_tos-tas-mosaic_spectral-change.rds")

## split data into tropics and temperate latitudes
tropical <- data.frame(rasterToPoints(mosaic_specchange)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)

temperate <- data.frame(rasterToPoints(mosaic_specchange)) %>%
  filter(y > 23.5 | y < -23.5) %>%
  rasterFromXYZ(.)

## summary stats and figures:
## mean
mean(values(tropical[[2]]), na.rm = T)
mean(values(temperate[[2]]), na.rm = T)

## standard deviation
sd(values(tropical[[2]]), na.rm = T)
sd(values(temperate[[2]]), na.rm = T)

## histogram:
df <- data.frame(latitudinal_region = rep("Tropical", length(values(tropical[[2]]))),
                 slope_spec_exp = c(values(tropical[[2]])))
df <- rbind(df, data.frame(latitudinal_region = rep("Temperate", length(values(temperate[[2]]))),
                           slope_spec_exp = c(values(temperate[[2]]))))
df <- filter(df, !is.na(slope_spec_exp))

df %>%
  ggplot(., aes(x = slope_spec_exp, group = latitudinal_region, fill = ..x..)) + 
  geom_histogram(position = position_dodge()) +
  theme_light() +
  labs(x = "Slope of spectral exponent",
       y = "Frequency",
       fill = "Slope of\nspectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  facet_wrap(~latitudinal_region)

## map it!
trop_df <- data.frame(rasterToPoints(tropical))
trop_df$realm = "Tropical"
temp_df <- data.frame(rasterToPoints(temperate))
temp_df$realm = "Temperate"
map_df <- rbind(trop_df, temp_df)

map_df %>%
  ggplot(., aes(x = x, y = y, fill = s_estimate)) + 
  geom_raster() +
  coord_fixed() +
  theme_void() +
  labs(fill = "Slope of\nspectral\nexponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  facet_wrap(~realm) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))  +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.1,
               aes(x=long, y=lat, group = group)) 

## boxplot:
df %>%
  ggplot(., aes(y = slope_spec_exp, x = latitudinal_region, fill = latitudinal_region)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  scale_fill_manual(values = pal_lat) +
  guides(fill = "none")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####           2c. Tropical vs. temperate latitudes x land vs. ocean                 #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specchange
mosaic_specchange <- readRDS("data-processed/01_tos-tas-mosaic_spectral-change.rds")

## split data into tropics and temperate latitudes
tropical <- data.frame(rasterToPoints(mosaic_specchange)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)

temperate <- data.frame(rasterToPoints(mosaic_specchange)) %>%
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

## summary stats and figures:
## mean
mean(values(tropical_land[[2]]), na.rm = T)
mean(values(temperate_land[[2]]), na.rm = T)
mean(values(tropical_ocean[[2]]), na.rm = T)
mean(values(temperate_ocean[[2]]), na.rm = T)

## standard deviation
sd(values(tropical_land[[2]]), na.rm = T)
sd(values(temperate_land[[2]]), na.rm = T)
sd(values(tropical_ocean[[2]]), na.rm = T)
sd(values(temperate_ocean[[2]]), na.rm = T)

## histogram - tropical vs. temperate 
df1 <- data.frame(Land = c(values(tropical_land[[2]])), 
                  Ocean = c(values(tropical_ocean[[2]])))
df1 <- gather(df1, key = "realm", value = "slope_spec_exp", Land, Ocean)
df1$trop_or_temp = "Tropical"
df2 <- data.frame(Land = c(values(temperate_land[[2]])),
                 Ocean = c(values(temperate_ocean[[2]])))
df2 <- gather(df2, key = "realm", value = "slope_spec_exp", Land, Ocean)
df2$trop_or_temp = "Temperate"
df <- rbind(df1, df2) %>% 
  filter(!is.na(slope_spec_exp))
df$lat_and_realm = paste(df$trop_or_temp, df$realm, sep = "_")

df %>%
  ggplot(., aes(x = slope_spec_exp, fill = lat_and_realm)) + 
  geom_histogram(position = position_dodge()) +
  theme_light() +
  labs(x = "Slope of spectral exponent",
       y = "Frequency",
       fill = "Latitude and realm") +
  facet_wrap(~realm) +
  scale_fill_discrete(labels = c("Temperate land", "Temperate ocean", "Tropical land", "Tropical ocean"))

## map it!
trop_land_df <- data.frame(rasterToPoints(tropical_land))
trop_land_df$lat_and_realm = "Tropical_land"
temp_land_df <- data.frame(rasterToPoints(temperate_land))
temp_land_df$lat_and_realm = "Temperate_land"
trop_ocean_df <- data.frame(rasterToPoints(tropical_ocean))
trop_ocean_df$lat_and_realm = "Tropical_ocean"
temp_ocean_df <- data.frame(rasterToPoints(temperate_ocean))
temp_ocean_df$lat_and_realm = "Temperate_ocean"
map_df <- rbind(temp_land_df,trop_ocean_df, trop_land_df, temp_ocean_df)

map_df %>%
  ggplot(., aes(x = x, y = y, fill = lat_and_realm)) + 
  geom_raster() +
  coord_fixed() +
  theme_void() +
  labs(fill = "Slope of\nspectral\nexponent") +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6))  

## boxplot:
df %>%
  ggplot(., aes(y = slope_spec_exp, x = lat_and_realm, fill = lat_and_realm)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  guides(fill = "none") +
  scale_x_discrete(labels = 
                        c("Temperate land", "Temperate ocean", "Tropical land", "Tropical ocean"))




