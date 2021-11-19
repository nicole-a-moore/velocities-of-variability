## compare sea surface vs air surface temperature analysis for one gcm 
library(tidyverse)
library(raster)
select <- dplyr::select

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

## plot mean difference between tas and tos
mean_tas_l <- calc(tas_l, mean)
mean_tas_s <- calc(tas_s, mean)

mean_tos_l <- calc(tos_l, mean)
mean_tos_s <- calc(tos_s, mean)

## plot mean daily tas and tos in ggplot
pts_tas <- data.frame(rasterToPoints(mean_tas_s))
pts_tos <- data.frame(rasterToPoints(mean_tos_s))
colnames(pts_tas)[3] <- "Air surface temperature"
colnames(pts_tos)[3] <- "Sea surface temperature"
pts <- left_join(pts_tas, pts_tos)

pts <- gather(pts, key = 'temp_type', value = 'mean_temp', "Air surface temperature", 
              "Sea surface temperature") %>%
  filter(!is.na(mean_temp))

pts %>%
  ggplot(., aes(x = x, y = y, fill = mean_temp)) + geom_raster() +
  coord_fixed() +
  facet_wrap(~temp_type) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Mean seasonally detrended\ndaily temperature\n(1871-2100)") +
  theme(panel.grid = element_blank())

## plot standard dev in daily temp for tas and tos
sd_tas_l <- calc(tas_l, sd)
sd_tas_s <- calc(tas_s, sd)

sd_tos_l <- calc(tos_l, sd)
sd_tos_s <- calc(tos_s, sd)

## plot difference in standard dev between tas and tos in ggplot
pts_tas <- data.frame(rasterToPoints(sd_tas_s))
pts_tos <- data.frame(rasterToPoints(sd_tos_s))
colnames(pts_tas)[3] <- "Air surface temperature"
colnames(pts_tos)[3] <- "Sea surface temperature"
pts <- left_join(pts_tas, pts_tos)

pts <- gather(pts, key = 'temp_type', value = 'sd', "Air surface temperature", 
              "Sea surface temperature") %>%
  filter(!is.na(sd))
  
pts %>%
  ggplot(., aes(x = x, y = y, fill = sd)) + geom_raster() +
  coord_fixed() +
  facet_wrap(~temp_type) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Standard deviation\nof seasonally detrended\ndaily temperature\n(1871-2100)") +
  theme(panel.grid = element_blank())


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

## calculate mean spectral exponent for tos and tas
mean_spec_tas_l <- calc(l_stack_tas, mean)
mean_spec_tas_s <- calc(s_stack_tas, mean)

mean_spec_tos_l <- calc(l_stack_tos, mean)
mean_spec_tos_s <- calc(s_stack_tos, mean)

## plot mean spectral exponent for tos and tas
pts_tas <- data.frame(rasterToPoints(mean_spec_tas_s))
pts_tos <- data.frame(rasterToPoints(mean_spec_tos_s))
colnames(pts_tas)[3] <- "Air surface temperature"
colnames(pts_tos)[3] <- "Sea surface temperature"
pts <- left_join(pts_tas, pts_tos)

pts <- gather(pts, key = 'temp_type', value = 'mean_spec_exp', "Air surface temperature", 
              "Sea surface temperature") %>%
  filter(!is.na(mean_spec_exp))

pts %>%
  ggplot(., aes(x = x, y = y, fill = mean_spec_exp)) + geom_raster() +
  coord_fixed() +
  facet_wrap(~temp_type) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Mean spectral exponent\n of seasonally detrended\ndaily temperature\n(1871-2100)") +
  theme(panel.grid = element_blank())
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
all_tas <- tas[,c(1,2,4)]
colnames(all_tas)[3] <- "Air surface temperature"
all_tos <- tos[,c(1,2,4)]
colnames(all_tos)[3] <- "Sea surface temperature"

all <- left_join(all_tas, all_tos)

all <- gather(all, key = 'temp_type', value = 'spec_exp_slope', "Air surface temperature", 
              "Sea surface temperature") %>%
  filter(!is.na(spec_exp_slope))

all %>%
  ggplot(., aes(x = lon, y = lat, fill = spec_exp_slope)) + geom_raster() +
  coord_fixed() +
  facet_wrap(~temp_type) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Slope of spectral\nexponent") +
  theme(panel.grid = element_blank()) +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") 

tas_sig <- mutate(tas, s_estimate = ifelse(s_p.value < 0.05, s_estimate, NA))
tos_sig <- mutate(tos, s_estimate = ifelse(s_p.value < 0.05, s_estimate, NA))
sig_tas <- tas_sig[,c(1,2,4)]
colnames(sig_tas)[3] <- "Air surface temperature"
sig_tos <- tos_sig[,c(1,2,4)]
colnames(sig_tos)[3] <- "Sea surface temperature"

sig <- left_join(sig_tas, sig_tos)

sig <- gather(sig, key = 'temp_type', value = 'sig_spec_exp_slope', "Air surface temperature", 
              "Sea surface temperature")

sig %>%
  ggplot(., aes(x = lon, y = lat, fill = sig_spec_exp_slope)) + geom_raster() +
  coord_fixed() +
  facet_wrap(~temp_type) +
  theme_minimal() +
  labs(x = "", y = "", fill = "Slope of spectral\nexponent") +
  theme(panel.grid = element_blank()) +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white")  +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group))


## where does direction agree between air and sea surface temperature?
combined_tos <- tos %>%
  rename("s_estimate_tos" = s_estimate) %>%
  select(lon, lat, s_estimate_tos) 
combined_tas <- tas %>%
  rename("s_estimate_tas" = s_estimate) %>%
  select(lon, lat, s_estimate_tas) 

combined <- left_join(combined_tas, combined_tos)

combined$same_direction <- ifelse(sign(combined$s_estimate_tos) == sign(combined$s_estimate_tas),
                                    "Same", "Different")

combined %>%
  filter(!is.na(same_direction)) %>%
  ggplot(., aes(x = lon, y = lat, fill = same_direction)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Direction\nof change")


## how different is the global average spectral change across the ocean when you use air vs. sea temp?
## convert raster layer to data frame
df_tas <- data.frame(rasterToPoints(s_stack_tas))
df_tas <- gather(df_tas, key = "window_number", value = "spec_exp", c(3:ncol(df_tas)))

tas_average <- df_tas %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Air surface temperature" = mean_spec_exp) 
tas_average$year = seq(1871, 2081, by = 10)

df_tos <- data.frame(rasterToPoints(s_stack_tos))
df_tos <- gather(df_tos, key = "window_number", value = "spec_exp", c(3:ncol(df_tos)))

tos_average <- df_tos %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Sea surface temperature" = mean_spec_exp)
tos_average$year = seq(1871, 2081, by = 10)

average <- left_join(tas_average, tos_average)
average <- gather(average, key = 'temp_type', value = 'mean_spec_exp', "Air surface temperature", 
                  "Sea surface temperature") %>%
  filter(!is.na(mean_spec_exp))

average %>%
  ggplot(., aes(x = year, y = mean_spec_exp)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent across ocean") + 
  facet_wrap(~temp_type) +
  geom_smooth(method = "lm", colour = "black")

average %>%
  filter(temp_type == "Sea surface temperature") %>% 
  lm(mean_spec_exp ~ year,
   data = .)

average %>%
  filter(temp_type == "Air surface temperature") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

