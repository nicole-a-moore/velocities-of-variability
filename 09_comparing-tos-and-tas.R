## compare sea surface vs air surface temperature analysis for one gcm 
library(tidyverse)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####  Compare air and sea surface temperature  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in mosaics:
tas_l <- readRDS("vov-shiny/l_mosaic-01_CMCC-CESM.rds")
tas_s <- readRDS("vov-shiny/s_mosaic-01_CMCC-CESM.rds")
tos_l <- readRDS("vov-shiny/l_mosaic-01_CMCC-CESM_tos.rds")
tos_s <- readRDS("vov-shiny/s_mosaic-01_CMCC-CESM_tos.rds")
extent(tas_l) <- c(-180, 180, -90, 90)
extent(tas_s) <- c(-180, 180, -90, 90)

## crop tas data to ocean only
terr <- raster("data-processed/raster_terr_mask.nc") 
terr <- crop(terr, tas_l)
tas_l <- mask(tas_l, terr, inverse = T)
tas_s <- mask(tas_s, terr, inverse = T)
plot(tas_l[[1]])
plot(tas_s[[1]])

## plot different between tas and tos
mean_tas_l <- calc(tas_l, mean)
plot(mean_tas_l)
mean_tas_s <- calc(tas_s, mean)
plot(mean_tas_s)

mean_tos_l <- calc(tos_l, mean)
plot(mean_tos_l)
mean_tos_s <- calc(tos_s, mean)
plot(mean_tos_s)

## mask so missing areas in tos are missing from tas
mean_tas_l <- mask(mean_tas_l, mean_tos_l)
mean_tas_s <- mask(mean_tas_s, mean_tos_s)

## sea surface temperature less variable - as predicted



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#####  Compare spectral exponent of air and sea surface temperature  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
terr <- raster("data-processed/raster_terr_mask.nc") 
## extract 10 year time window 
l_stack_tas <- readRDS("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/l_stack_list.rds")[[6]] 
s_stack_tas <- readRDS("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/s_stack_list.rds")[[6]] 
extent(l_stack_tas) <- c(-179, 179, -89, 89)
extent(s_stack_tas) <- c(-179, 179, -89, 89)

l_stack_tos <- readRDS("data-raw/01_CMCC-CESM/l_stack_list_tos.rds")[[6]] 
s_stack_tos <- readRDS("data-raw/01_CMCC-CESM/s_stack_list_tos.rds")[[6]] 

l_stack_tas <- crop(l_stack_tas, l_stack_tos)
s_stack_tas <- crop(s_stack_tas, s_stack_tos)
l_stack_tos <- crop(l_stack_tos, l_stack_tas)
s_stack_tos <- crop(s_stack_tos, s_stack_tas)

terr <- crop(terr, l_stack_tas[[1]])

## mask to only include ocean, only areas in tos
l_stack_tas <- mask(l_stack_tas, terr, inverse = T)
s_stack_tas <- mask(s_stack_tas, terr, inverse = T)
plot(l_stack_tas[[1]])
plot(s_stack_tas[[1]])
l_stack_tas <- mask(l_stack_tas, l_stack_tos)
s_stack_tas <- mask(s_stack_tas, s_stack_tos)
plot(l_stack_tas[[1]])
plot(s_stack_tas[[1]])

## take a looksie:
plot(l_stack_tas[[1]]) ## higher spectral exponent (shallower, less autocorrelated, white noise)
plot(l_stack_tos[[1]]) ## lower spectral exponent (steeper, more autocorrelated, red noise)

## subtract spectral exponent of tos from spectral exponent of tas
diff_l <- l_stack_tas - l_stack_tos
diff_s <- s_stack_tas - s_stack_tos

plot(diff_l[[1]])
plot(diff_s[[1]])
## spec exp of tas > spec exp of tos
## smaller spectral exponent = more autocorrelated

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#####  Compare change in spectral exponent of air and sea surface temperate  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###

##### Organize tas data in ocean #####
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

## crop data to ocean only
terr <- raster("data-processed/raster_terr_mask.nc") 
tas <- select(tas, lon, lat, everything())
raster_tas <- rasterFromXYZ(tas)
extent(raster_tas) <- c(-179, 179, -89, 89)
terr <- crop(terr, raster_tas)
raster_tas <- mask(raster_tas, terr, inverse = T)
plot(raster_tas[[1]])

## turn back into data frame for ggplot
tas <- data.frame(rasterToPoints(raster_tas))
colnames(tas)[1:2] <- c("lon", "lat")

l <- tas %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of\nspectral exponent (tas)") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

## remove non-significant slopes
l_sig <- tas %>%
  mutate(l_estimate = ifelse(l_p.value < 0.05, l_estimate, NA)) %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of\nspectral exponent (tas)") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

s <- tas %>%
  ggplot(., aes(x = lon, y = lat, fill = s_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of\nspectral exponent (tas)") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

s_sig <- tas %>%
  mutate(s_estimate = ifelse(l_p.value < 0.05, s_estimate, NA)) %>%
  ggplot(., aes(x = lon, y = lat, fill = s_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of\nspectral exponent (tas)") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

##### Organize tos data  #####
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

l <- tos %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of\nspectral exponent (tos)") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

## remove non-significant slopes
l_sig <- tos %>%
  mutate(l_estimate = ifelse(l_p.value < 0.05, l_estimate, NA)) %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of\nspectral exponent (tos)") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

s <- tos %>%
  ggplot(., aes(x = lon, y = lat, fill = s_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of\nspectral exponent (tos)") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

s_sig <- tos %>%
  mutate(s_estimate = ifelse(l_p.value < 0.05, s_estimate, NA)) %>%
  ggplot(., aes(x = lon, y = lat, fill = s_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of\nspectral exponent (tos)") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 




##### Compare magnitude and direction of spectral change #####
colnames(tos)[3:6] <- paste(colnames(tos)[3:6], "_tos", sep = "")
colnames(tas)[3:6] <- paste(colnames(tas)[3:6], "_tas", sep = "")
colnames(tas)[1:2] <- c("lon", "lat")
combined <- left_join(tos, tas)

## where does direction agree between air and sea surface temperature?
combined$s_same_direction <- ifelse(sign(combined$s_estimate_tos) == sign(combined$s_estimate_tas),
                                    "Same", "Different")
combined$l_same_direction <- ifelse(sign(combined$l_estimate_tos) == sign(combined$l_estimate_tas),
                                    "Same", "Different")

combined %>%
  ggplot(., aes(x = lon, y = lat, fill = l_same_direction)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Does direction\nof magnitude agree?")

combined %>%
  ggplot(., aes(x = lon, y = lat, fill = s_same_direction)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Does direction of magnitude agree?")

## how different is the magnitude?
combined$s_tas_minus_tos <- combined$s_estimate_tas - combined$s_estimate_tos
combined$l_tas_minus_tos <- combined$l_estimate_tas - combined$l_estimate_tos

combined %>%
  ggplot(., aes(x = lon, y = lat, fill = s_tas_minus_tos)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", 
       fill = "Difference\nbetween air\nand sea surface\ntemperature\nspectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") 

combined %>%
  ggplot(., aes(x = lon, y = lat, fill = l_tas_minus_tos)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", 
       fill = "Difference\nbetween air\nand sea surface\ntemperature\nspectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") 


###### Create mosaic of tos and tas ######
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

## mosaic together tas and tos 
raster_tas <- rasterFromXYZ(tas)
raster_tos <- rasterFromXYZ(tos)
extent(raster_tas) <- c(-179, 179, -89, 89)

isna <- is.na(raster_tos[[1]])
raster_tos[isna == 1] <- raster_tas[isna == 1]

mosaic <- raster_tos
plot(mosaic$l_estimate_tos)
plot(mosaic$s_estimate_tos)

## save:
saveRDS(mosaic, "data-processed/01_tos-tas-mosaic.rds")
