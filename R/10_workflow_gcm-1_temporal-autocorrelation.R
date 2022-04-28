## whole workflow for first GCM  !!!
library(tidyverse)
library(raster)
library(PNWColors)
library(gridExtra)
library(cowplot)
select <- dplyr::select
countries <- map_data("world")

## create colour palette:
pal = pnw_palette("Sunset",5, type = "discrete")
pal_realm <- c(pal[1], pal[3])
pal_lat <- c(pal[2], pal[4])
pal = pnw_palette("Shuksan2",5, type = "discrete")
pal_lat_and_realm <-  c(pal[5], pal[4], pal[1], pal[2])

## for now: just do for seasonally detrended data, not linearly detrended
## create data that has tas on land and tos in the ocean
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       Create mosaics of tos and tas       #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####       spectral exponent      #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## extract 10 year time window 
s_stack_tas_PSD_low <- readRDS("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/s_stack_list_PSD_low.rds")[[6]] 
s_stack_tas_AWC <- readRDS("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/s_stack_list_AWC.rds")[[6]] 
s_stack_tas_PSD_high <- readRDS("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/s_stack_list_PSD_high.rds")[[6]] 

s_stack_tos_PSD_low <- readRDS("/Volumes/SundayLab/CMIP5-GCMs_tos/01_CMCC-CESM/s_stack_list_PSD_low_tos.rds")[[6]] 
s_stack_tos_AWC <- readRDS("/Volumes/SundayLab/CMIP5-GCMs_tos/01_CMCC-CESM/s_stack_list_AWC_tos.rds")[[6]] 
s_stack_tos_PSD_high <- readRDS("/Volumes/SundayLab/CMIP5-GCMs_tos/01_CMCC-CESM/s_stack_list_PSD_high_tos.rds")[[6]] 

## crop to same extent
list_tas <- c(s_stack_tas_PSD_low, s_stack_tas_AWC, s_stack_tas_PSD_high)
list_tas <- sapply(list_tas, FUN = crop, s_stack_tos_PSD_low)
list_tas <- sapply(list_tas, FUN = stack)

list_tos <- c(s_stack_tos_PSD_low, s_stack_tos_AWC, s_stack_tos_PSD_high)
list_tos <- sapply(list_tos, FUN = crop, list_tas[[1]])
list_tos <- sapply(list_tos, FUN = stack)

land <- raster("data-processed/masks/cmip5-land.grd") 
extent(land) <- c(0, 360, -90, 90)
land <- crop(land, list_tas[[1]])
land[land != 1] <- NA

list_tas <- sapply(list_tas, FUN = mask, land)
list_tas <- sapply(list_tas, FUN = stack)

## mosaic rasters together 
s_stack_tas_PSD_low <- list_tas[[1]]
s_stack_tas_AWC <- list_tas[[2]]
s_stack_tas_PSD_high <- list_tas[[3]]

s_stack_tos_PSD_low <- list_tos[[1]]
s_stack_tos_AWC <- list_tos[[2]]
s_stack_tos_PSD_high <- list_tos[[3]]

isna <- is.na(s_stack_tos_PSD_low[[1]])
s_stack_tos_PSD_low[isna == 1] <- s_stack_tas_PSD_low[isna == 1]
s_stack_tos_AWC[isna == 1] <- s_stack_tas_AWC[isna == 1]
s_stack_tos_PSD_high[isna == 1] <- s_stack_tas_PSD_high[isna == 1]

mosaic_specexp_low <- s_stack_tos_PSD_low
mosaic_specexp_AWC <- s_stack_tos_AWC
mosaic_specexp_high <- s_stack_tos_PSD_high

## keep track of land and ocean by making a mask:
land_mask <- s_stack_tos_PSD_low
land_mask[land_mask != mosaic_specexp_low] <- NA
land_mask[land_mask == mosaic_specexp_low] <- 1
land_mask <- land_mask[[1]]

## save:
#saveRDS(mosaic_specexp, "data-processed/01_tos-tas-mosaic_spectral-exponent.rds")
#saveRDS(land_mask, "data-processed/01_tos-tas-land-mask.rds")

saveRDS(mosaic_specexp_low, "data-processed/01_tos-tas-mosaic_spectral-exponent_low.rds")
saveRDS(mosaic_specexp_high, "data-processed/01_tos-tas-mosaic_spectral-exponent_high.rds")
saveRDS(mosaic_specexp_AWC, "data-processed/01_tos-tas-mosaic_spectral-exponent_AWC.rds")

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
  select(lon, lat, l_estimate_PSD_low, s_estimate_PSD_low, l_estimate_PSD_high, s_estimate_PSD_high) %>%
  unique() %>%
  mutate(lon = ifelse(lon >= 180, lon - 180, lon + 178))

### get tos data
path = "/Volumes/SundayLab/CMIP5-GCMs_tos/" 

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
  select(lon, lat, l_estimate_PSD_low, s_estimate_PSD_low, l_estimate_PSD_high, s_estimate_PSD_high) %>%
  unique()

## mosaic together tas and tos 
raster_tas <- rasterFromXYZ(tas)
raster_tos <- rasterFromXYZ(tos)
raster_tas <- extend(raster_tas, c(0, 360, -90, 90))
raster_tos <- extend(raster_tos, c(0, 360, -90, 90))

## crop
raster_tas <- crop(raster_tas, raster_tos)
raster_tos <- crop(raster_tos, raster_tas)

isna <- is.na(raster_tos[[1]])
raster_tos[isna == 1] <- raster_tas[isna == 1]

mosaic_specchange <- raster_tos
plot(mosaic_specchange$s_estimate_PSD_low)

## save:
saveRDS(mosaic_specchange, "data-processed/01_tos-tas-mosaic_spectral-change_multifrac.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 1a. Analyze spectral exponent on land vs. ocean #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specexp
mosaic_specexp_low <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent_low.rds")
mosaic_specexp_high <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent_high.rds")
mosaic_specexp <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent.rds")
land_mask <-readRDS("data-processed/01_tos-tas-land-mask.rds")

legend_rast <- land
legend_rast[is.na(land)] <- 2

land <- extend(land, c(0, 358, -78, 89))
land <- crop(land, c(0, 358, -78, 89))


## split data into land and ocean
land_low <- mask(mosaic_specexp_low, land)
plot(land_low[[1]])
ocean_low <- mask(mosaic_specexp_low, land, inverse = TRUE)
plot(ocean_low[[1]])

land_high <- mask(mosaic_specexp_high, land)
plot(land_high[[1]])
ocean_high <- mask(mosaic_specexp_high, land, inverse = TRUE)
plot(ocean_high[[1]])

land_both <- mask(mosaic_specexp, land_mask)
plot(land_both[[1]])
ocean_both <- mask(mosaic_specexp, land_mask, inverse = TRUE)
plot(ocean_both[[1]])

## baseline:
max(values(land_low[[1]]), na.rm=TRUE)
min(values(land_low[[1]]),  na.rm=TRUE)
mean(values(land_low[[1]]), na.rm=TRUE)
max(values(ocean_low[[1]]), na.rm=TRUE)
min(values(ocean_low[[1]]),  na.rm=TRUE)
mean(values(ocean_low[[1]]), na.rm=TRUE)


## summary stats and figures:
## mean
mean(values(land_low), na.rm=TRUE)
mean(values(ocean_low), na.rm=TRUE)
mean(values(land_high), na.rm=TRUE)
mean(values(ocean_high), na.rm=TRUE)
mean(values(land_both), na.rm=TRUE)
mean(values(ocean_both), na.rm=TRUE)

## standard deviation
sd(values(land_low), na.rm=TRUE)
sd(values(ocean_low), na.rm=TRUE)
sd(values(land_high), na.rm=TRUE)
sd(values(ocean_high), na.rm=TRUE)
sd(values(land_both), na.rm=TRUE)
sd(values(ocean_both), na.rm=TRUE)

## histogram:
mean_land_low <- calc(land_low, mean)
mean_ocean_low <- calc(ocean_low, mean)
mean_land_high <- calc(land_high, mean)
mean_ocean_high <- calc(ocean_high, mean)
mean_land_both <- calc(land_both, mean)
mean_ocean_both <- calc(ocean_both, mean)


low_or_high = append(append(rep("low", length(values(mean_land_low))), 
                            rep("high", length(values(mean_land_high)))),
                     append(rep("low", length(values(mean_ocean_low))), 
                            rep("high", length(values(mean_ocean_high)))))
low_or_high <- append(low_or_high, append(rep("both", length(values(mean_land_both))), 
                                          rep("both", length(values(mean_ocean_both)))))
realm =  append(rep("Land", length(values(mean_land_low)) + length(values(mean_land_high))), 
                rep("Ocean", length(values(mean_ocean_low)) + length(values(mean_ocean_high))))
realm <- append(realm, append(rep("Land", length(values(mean_land_both))), 
                       rep("Ocean", length(values(mean_ocean_both)))))
values = append(append(values(mean_land_low), values(mean_land_high)),
                append(values(mean_ocean_low), values(mean_ocean_high)))
values = append(values,  append(-values(mean_land_both), -values(mean_ocean_both)))

df <- data.frame(mean_spec_exp = values,
                 low_or_high = low_or_high,
                 realm = realm)
df <- filter(df, !is.na(mean_spec_exp)) %>%
  mutate(group = paste(realm, low_or_high, sep = ""))

hist <- df %>%
  filter(low_or_high %in% c("both", "low")) %>%
  ggplot(., aes(x = mean_spec_exp, fill = realm)) + 
  geom_histogram(position = position_dodge(), binwidth = 0.05) +
  theme_light() +
  labs(x = "", y = "Frequency", fill = "Realm") + #x = "Mean local spectral exponent"
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none") + coord_flip() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.title.x = element_text(size = 10),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = -0.9, unit = "cm")) +
  facet_wrap(~low_or_high)

map_legend <- data.frame(rasterToPoints(legend_rast)) %>%
  mutate(layer = as.factor(layer)) %>%
  ggplot(., aes(x = x, y = y, fill = layer)) +
  geom_raster() +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none")

## boxplot:
box <- df %>%
  filter(low_or_high %in% c("both", "high")) %>%
  ggplot(., aes(y = mean_spec_exp, x = realm, fill = realm)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none") +
  scale_x_discrete(labels = c("Land", "Ocean")) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3) +
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")) +
  facet_wrap(~low_or_high)

## arrange plots side by side:
lay <- rbind(c(1,1,1,1,2,2),
             c(1,1,1,1,2,2))
besties <- grid.arrange(box, hist, layout_matrix = lay)

ggsave(besties, path = "figures/analysis-workflow", filename = "o-vs-l_mean-local-spec-exp_mf.png", 
       device = "png", width = 6, height = 6)


## correlation between low and high exponent?
land_df_low <- data.frame(rasterToPoints(mean_land_low))
land_df_low$realm = "Land"
land_df_low$low_or_high = "low"
land_df_high <- data.frame(rasterToPoints(mean_land_high))
land_df_high$realm = "Land"
land_df_high$low_or_high = "high"
ocean_df_low <- data.frame(rasterToPoints(mean_ocean_low))
ocean_df_low$realm = "Ocean"
ocean_df_low$low_or_high = "low"
ocean_df_high <- data.frame(rasterToPoints(mean_ocean_high))
ocean_df_high$realm = "Ocean"
ocean_df_high$low_or_high = "high"

## rename columns 
colnames(land_df_low)[3] = colnames(land_df_high)[3] = colnames(ocean_df_low)[3] = 
  colnames(ocean_df_high)[3] = "mean_spec_exp"

map_df <- rbind(ocean_df_low, land_df_low) %>%
  rbind(., ocean_df_high) %>%
  rbind(., land_df_high) %>%
  mutate(group = paste(realm, low_or_high, sep = ""))

map_df %>%
  filter(low_or_high != "both") %>%
  select(-group) %>%
  spread(key = "low_or_high", value = "mean_spec_exp") %>%
  ggplot(., aes(x = low, y = high, colour = realm)) + geom_point(size = 0.1) +
  facet_wrap(~realm) +
  scale_colour_manual(values = pal_realm) +
  geom_smooth( method = "lm", colour = "red")+
  theme_light()


## NOT ADAPTED FOR MF ANALYSIS YET
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 1b. Analyze spectral exponent across tropical vs. temperate latitudes #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specexp
mosaic_specexp_low <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent_low.rds")
mosaic_specexp_high <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent_high.rds")

## split data into tropics and temperate latitudes
tropical <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)

temperate <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y > 23.5 | y < -23.5) %>%
  rasterFromXYZ(.)

## summary stats and figures:
## mean
mean(values(tropical), na.rm=TRUE)
mean(values(temperate), na.rm=TRUE)

## standard deviation
sd(values(tropical), na.rm=TRUE)
sd(values(temperate), na.rm=TRUE)

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

## NOT ADAPTED FOR MF ANALYSIS YET
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
tropical_ocean <- mask(tropical, tropical_land, inverse = TRUE)
plot(tropical_land[[1]])
plot(tropical_ocean[[1]])

temperate_land <- crop(land_mask, temperate)
temperate_land <- mask(temperate, temperate_land)
temperate_ocean <- mask(temperate, temperate_land, inverse = TRUE)
plot(temperate_land[[1]])
plot(temperate_ocean[[1]])

## summary stats and figures:
## mean
mean(values(tropical_land), na.rm=TRUE)
mean(values(temperate_land), na.rm=TRUE)
mean(values(tropical_ocean), na.rm=TRUE)
mean(values(temperate_ocean), na.rm=TRUE)

## standard deviation
sd(values(tropical_land), na.rm=TRUE)
sd(values(temperate_land), na.rm=TRUE)
sd(values(tropical_ocean), na.rm=TRUE)
sd(values(temperate_ocean), na.rm=TRUE)

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

## reorder factor levels
df$lat_and_realm <- factor(df$lat_and_realm, levels = c("Temperate_Land","Tropical_Land","Temperate_Ocean","Tropical_Ocean"), 
                           ordered = TRUE)
hist <- df %>%
  ggplot(., aes(x = mean_spec_exp, fill = lat_and_realm)) + 
  geom_histogram(position = position_dodge(), binwidth = 0.05) +
  theme_light() +
  labs(x = "", y = "Frequency", fill = "Latitude and realm") + #x = "Mean local spectral exponent"
  scale_fill_manual(values = pal_lat_and_realm, 
                    labels = c("Temperate land", "Tropical land", "Temperate ocean", "Tropical ocean")) +
  coord_flip() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.title.x = element_text(size = 10),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = -0.9, unit = "cm")) 

## extract legend 
legend <- get_legend(hist)
hist <- hist + guides(fill = "none")

## map it!
trop_land_df <- data.frame(rasterToPoints(tropical_land))
trop_land_df$lat_and_realm = "Tropical_Land"
temp_land_df <- data.frame(rasterToPoints(temperate_land))
temp_land_df$lat_and_realm = "Temperate_Land"
trop_ocean_df <- data.frame(rasterToPoints(tropical_ocean))
trop_ocean_df$lat_and_realm = "Tropical_Ocean"
temp_ocean_df <- data.frame(rasterToPoints(temperate_ocean))
temp_ocean_df$lat_and_realm = "Temperate_Ocean"
map_df <- rbind(temp_land_df,trop_ocean_df, trop_land_df, temp_ocean_df)

## reorder factor levels
map_df$lat_and_realm <- factor(map_df$lat_and_realm, 
                               levels = c("Temperate_Land","Tropical_Land","Temperate_Ocean","Tropical_Ocean"), 
                           ordered = TRUE)

map_df %>%
  ggplot(., aes(x = x, y = y, fill = lat_and_realm)) + 
  geom_raster() +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = pal_lat_and_realm) +
  labs(fill = "Mean local\nspectral\nexponent") +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) 
  

## boxplot:
box <- df %>%
  ggplot(., aes(y = mean_spec_exp, x = lat_and_realm, fill = lat_and_realm)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  scale_fill_manual(values = pal_lat_and_realm) +
  guides(fill = "none")+
  scale_x_discrete(labels = 
                     c("Temperate land", "Tropical land", "Temperate ocean", "Tropical ocean")) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3) +
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))

## arrange plots side by side:
lay <- rbind(c(1,1,1,1,2,2))
besties <- grid.arrange(box, hist, layout_matrix = lay)

ggsave(besties, path = "figures/analysis-workflow", filename = "o-vs-l-across-lat_mean-local-spec-exp.png", 
       device = "png", width = 8, height = 4)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 2a. Analyze change in spectral exponent on land vs. ocean #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specchange
mosaic_specchange_mf <- readRDS("data-processed/01_tos-tas-mosaic_spectral-change_multifrac.rds")
mosaic_specchange <- readRDS("data-processed/01_tos-tas-mosaic_spectral-change.rds")

land <- extend(land, c(0, 360, -90, 90))

## split data into land and ocean
land_low <- mask(mosaic_specchange_mf$s_estimate_PSD_low, land)
plot(land_low[[1]])
land_high <- mask(mosaic_specchange_mf$s_estimate_PSD_high, land)
plot(land_high[[1]])
land_both <- mask(mosaic_specchange$s_estimate, land)
plot(land_both[[1]])
ocean_low <- mask(mosaic_specchange_mf$s_estimate_PSD_low, land, inverse = TRUE)
plot(ocean_low[[1]])
ocean_high <- mask(mosaic_specchange_mf$s_estimate_PSD_high, land, inverse = TRUE)
plot(ocean_high[[1]])
ocean_both <- mask(mosaic_specchange$s_estimate, land, inverse = TRUE)
plot(ocean_both[[1]])

max(values(land_low[[1]]), na.rm=TRUE)
min(values(land_low[[1]]),  na.rm=TRUE)
mean(values(land_low[[1]]), na.rm=TRUE)
max(values(ocean_low[[1]]), na.rm=TRUE)
min(values(ocean_low[[1]]),  na.rm=TRUE)
mean(values(ocean_low[[1]]), na.rm=TRUE)

## summary stats and figures:
## mean
mean(values(land_low), na.rm=TRUE)
mean(values(ocean_low), na.rm=TRUE)
mean(values(land_high), na.rm=TRUE)
mean(values(ocean_high), na.rm=TRUE)
mean(values(land_both), na.rm=TRUE)
mean(values(ocean_both), na.rm=TRUE)

## standard deviation
sd(values(land_low), na.rm=TRUE)
sd(values(ocean_low), na.rm=TRUE)
sd(values(land_high), na.rm=TRUE)
sd(values(ocean_high), na.rm=TRUE)
sd(values(land_both), na.rm=TRUE)
sd(values(ocean_both), na.rm=TRUE)

## histogram:
low_or_high = append(append(rep("low", length(values(land_low))), 
                            rep("high", length(values(land_high)))),
                     append(rep("low", length(values(ocean_low))), 
                            rep("high", length(values(ocean_high)))))
low_or_high <- append(low_or_high, append(rep("both", length(values(land_both))), 
                                          rep("both", length(values(ocean_both)))))
realm =  append(rep("Land", length(values(land_low)) + length(values(land_high))), 
                rep("Ocean", length(values(ocean_low)) + length(values(ocean_high))))
realm <- append(realm, append(rep("Land", length(values(land_both))), 
                              rep("Ocean", length(values(ocean_both)))))
values <- append(append(values(land_low), values(land_high)),
                 append(values(ocean_low), values(ocean_high)))
values <- append(values, append(-values(land_both), -values(ocean_both)))

df <- data.frame(slope_spec_exp = values,
                 low_or_high = low_or_high,
                 realm = realm)
df <- filter(df, !is.na(slope_spec_exp)) %>%
  mutate(group = paste(realm, low_or_high, sep = ""))

df %>%
  ggplot(., aes(x = slope_spec_exp, group = group, fill = ..x..)) + 
  geom_histogram(position = position_dodge()) +
  theme_light() +
  labs(x = "Slope of spectral exponent",
       y = "Frequency",
       fill = "Slope of\nspectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  facet_wrap(low_or_high~realm, nrow = 3) +
  geom_vline(xintercept = 0) 

## map it!
land_df_low <- data.frame(rasterToPoints(land_low))
land_df_low$realm = "Land"
land_df_low$low_or_high = "low"
land_df_high <- data.frame(rasterToPoints(land_high))
land_df_high$realm = "Land"
land_df_high$low_or_high = "high"
ocean_df_low <- data.frame(rasterToPoints(ocean_low))
ocean_df_low$realm = "Ocean"
ocean_df_low$low_or_high = "low"
ocean_df_high <- data.frame(rasterToPoints(ocean_high))
ocean_df_high$realm = "Ocean"
ocean_df_high$low_or_high = "high"
ocean_df_both <- data.frame(rasterToPoints(ocean_both))
ocean_df_both$x = ocean_df_both$x +180
ocean_df_both$realm = "Ocean"
ocean_df_both$low_or_high = "both" 
land_df_both <- data.frame(rasterToPoints(land_both))
land_df_both$x = land_df_both$x +180
land_df_both$realm = "Land"
land_df_both$low_or_high = "both"

## rename columns 
colnames(land_df_low)[3] = colnames(land_df_high)[3] = colnames(ocean_df_low)[3] = 
  colnames(ocean_df_high)[3] = colnames(land_df_both)[3] = colnames(ocean_df_both)[3] = "slope_spec_exp"

map_df <- rbind(ocean_df_low, land_df_low) %>%
  rbind(., ocean_df_high) %>%
  rbind(., land_df_high) %>%
  rbind(., ocean_df_both) %>%
  rbind(., land_df_both) %>%
  mutate(group = paste(realm, low_or_high, sep = ""))

map_df %>%
  ggplot(., aes(x = x, y = y, fill = slope_spec_exp)) + 
  geom_raster() +
  coord_fixed() +
  theme_void() +
  labs(fill = "Slope of\nspectral\nexponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  facet_wrap(low_or_high~realm, nrow = 3) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) 

hist <- df %>%
  ggplot(., aes(x = slope_spec_exp, fill = realm)) + 
  geom_histogram(position = position_dodge(), binwidth = 0.0002) +
  theme_light() +
  labs(x = "", y = "Frequency", fill = "Realm") + #x = "Slope of spectral exponent"
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none") + coord_flip() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.title.x = element_text(size = 10),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = -0.9, unit = "cm")) +
  facet_wrap(~group, nrow = 4)
  
## boxplot:
box <- df %>%
  ggplot(., aes(y = slope_spec_exp, x = realm, fill = realm)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Slope of spectral exponent", x = "") +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none") +
  scale_x_discrete(labels = c("Land", "Ocean")) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3) +
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")) +
  facet_wrap(~low_or_high, nrow = 1) +
  geom_hline(yintercept = 0) 

## arrange plots side by side:
lay <- rbind(c(1,1,1,1,2,2),
             c(1,1,1,1,2,2))
besties <- grid.arrange(box, hist, layout_matrix = lay)

ggsave(besties, path = "figures/analysis-workflow", filename = "o-vs-l_spec-change_mf.png", 
       device = "png", width = 4, height = 6)


## when high freq. exponent is increasing, is low frequency exponent decreasing?
samesign_low <- land_low 
samesign_low[land_low < 0] <- -1
samesign_low[land_low > 0] <- 1
samesign_high <- land_high 
samesign_high[land_high < 0] <- -1
samesign_high[land_high > 0] <- 1
plot(samesign_low == samesign_high)

samesign_low <- ocean_low 
samesign_low[ocean_low < 0] <- -1
samesign_low[ocean_low > 0] <- 1
samesign_high <- ocean_high 
samesign_high[ocean_high < 0] <- -1
samesign_high[ocean_high > 0] <- 1
plot(samesign_low == samesign_high)

map_df %>%
  filter(low_or_high != "both") %>%
  select(-group) %>%
  spread(key = "low_or_high", value = "slope_spec_exp") %>%
  ggplot(., aes(x = low, y = high, colour = realm)) + geom_point(size = 0.1) +
  facet_wrap(~realm) +
  scale_colour_manual(values = pal_realm) +
  theme_light() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Slope of low-frequency exponent", y = "Slope of high-freqency exponent")



## how different is the global average spectral change across the ocean vs. land?
## use mosaic_specexp
mosaic_specexp_low <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent_low.rds")
mosaic_specexp_high <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent_high.rds")
mosaic_specexp <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent.rds")
mosaic_specexp_low <- extend(mosaic_specexp_low, land)
mosaic_specexp_high <- extend(mosaic_specexp_high, land)

## split data into land and ocean
land_low <- mask(mosaic_specexp_low, land)
plot(land_low[[1]])
ocean_low <- mask(mosaic_specexp_low, land, inverse = TRUE)
plot(ocean_low[[1]])

land_high <- mask(mosaic_specexp_high, land)
plot(land_high[[1]])
ocean_high <- mask(mosaic_specexp_high, land, inverse = TRUE)
plot(ocean_high[[1]])

land_both <- mask(mosaic_specexp, land_mask)
plot(land_both[[1]])
ocean_both <- mask(mosaic_specexp, land_mask, inverse = TRUE)
plot(ocean_both[[1]])

## convert raster layer to data frame
df_land_low <- data.frame(rasterToPoints(land_low))
df_land_low <- gather(df_land_low, key = "window_number", value = "spec_exp", c(3:ncol(df_land_low)))

df_land_high <- data.frame(rasterToPoints(land_high))
df_land_high <- gather(df_land_high, key = "window_number", value = "spec_exp", c(3:ncol(df_land_high)))

df_land_both <- data.frame(rasterToPoints(land_both))
df_land_both <- gather(df_land_both, key = "window_number", value = "spec_exp", c(3:ncol(df_land_both)))
df_land_both$x <- df_land_both$x + 180

land_average_low <- df_land_low %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Land" = mean_spec_exp) 
land_average_low$year = seq(1871, 2081, by = 10)

land_average_high <- df_land_high %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Land" = mean_spec_exp) 
land_average_high$year = seq(1871, 2081, by = 10)

land_average_both <- df_land_both %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Land" = mean_spec_exp) 
land_average_both$year = seq(1871, 2081, by = 10)

df_ocean_low <- data.frame(rasterToPoints(ocean_low))
df_ocean_low <- gather(df_ocean_low, key = "window_number", value = "spec_exp", c(3:ncol(df_ocean_low)))
df_ocean_high <- data.frame(rasterToPoints(ocean_high))
df_ocean_high <- gather(df_ocean_high, key = "window_number", value = "spec_exp", c(3:ncol(df_ocean_high)))
df_ocean_both <- data.frame(rasterToPoints(ocean_both))
df_ocean_both <- gather(df_ocean_both, key = "window_number", value = "spec_exp", c(3:ncol(df_ocean_both)))
df_ocean_both$x <- df_ocean_both$x + 180

ocean_average_low <- df_ocean_low %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Ocean" = mean_spec_exp)
ocean_average_low$year = seq(1871, 2081, by = 10)

ocean_average_high <- df_ocean_high %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Ocean" = mean_spec_exp)
ocean_average_high$year = seq(1871, 2081, by = 10)

ocean_average_both <- df_ocean_both %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Ocean" = mean_spec_exp)
ocean_average_both$year = seq(1871, 2081, by = 10)

average_low <- left_join(land_average_low, ocean_average_low)
average_low <- gather(average_low, key = 'realm', value = 'mean_spec_exp', "Land", 
                  "Ocean") %>%
  filter(!is.na(mean_spec_exp))

average_low %>%
  ggplot(., aes(x = year, y = mean_spec_exp, fill = realm, colour = realm)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") + 
  geom_smooth(method = "lm")  +
  scale_colour_manual(values = pal_realm) +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none", colour = "none")

average_low %>%
  filter(realm == "Ocean") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average_low %>%
  filter(realm == "Land") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average_high <- left_join(land_average_high, ocean_average_high)
average_high <- gather(average_high, key = 'realm', value = 'mean_spec_exp', "Land", 
                      "Ocean") %>%
  filter(!is.na(mean_spec_exp))

average_high %>%
  ggplot(., aes(x = year, y = mean_spec_exp, fill = realm, colour = realm)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") + 
  geom_smooth(method = "lm")  +
  scale_colour_manual(values = pal_realm) +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none", colour = "none")

average_high %>%
  filter(realm == "Ocean") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average_high %>%
  filter(realm == "Land") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average_high$type = "high"


average_both <- left_join(land_average_both, ocean_average_both)
average_both <- gather(average_both, key = 'realm', value = 'mean_spec_exp', "Land", 
                       "Ocean") %>%
  filter(!is.na(mean_spec_exp))

average_both %>%
  ggplot(., aes(x = year, y = mean_spec_exp, fill = realm, colour = realm)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") + 
  geom_smooth(method = "lm")  +
  scale_colour_manual(values = pal_realm) +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none", colour = "none")

average_both %>%
  filter(realm == "Ocean") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average_both %>%
  filter(realm == "Land") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average_both$type = "both"
average_both$mean_spec_exp = -average_both$mean_spec_exp

average_low %>%
  mutate(type = "low") %>%
  rbind(average_high, .) %>%
  rbind(average_both, .) %>%
  ggplot(., aes(x = year, y = mean_spec_exp, fill = realm, colour = realm,
                shape = type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") + 
  geom_smooth(method = "lm")  +
  scale_colour_manual(values = pal_realm) +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none", colour = "none") + facet_wrap(~realm)

saveRDS(average_low, "data-processed/temp_avg-low.rds")

## NOT ADAPTED FOR MF ANALYSIS YET
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
mean(values(tropical[[2]]), na.rm=TRUE)
mean(values(temperate[[2]]), na.rm=TRUE)

## standard deviation
sd(values(tropical[[2]]), na.rm=TRUE)
sd(values(temperate[[2]]), na.rm=TRUE)

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

## NOT ADAPTED FOR MF ANALYSIS YET
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
tropical_ocean <- mask(tropical, tropical_land, inverse = TRUE)
plot(tropical_land[[1]])
plot(tropical_ocean[[1]])

temperate_land <- crop(land_mask, temperate)
temperate_land <- mask(temperate, temperate_land)
temperate_ocean <- mask(temperate, temperate_land, inverse = TRUE)
plot(temperate_land[[1]])
plot(temperate_ocean[[1]])

## summary stats and figures:
## mean
mean(values(tropical_land[[2]]), na.rm=TRUE)
mean(values(temperate_land[[2]]), na.rm=TRUE)
mean(values(tropical_ocean[[2]]), na.rm=TRUE)
mean(values(temperate_ocean[[2]]), na.rm=TRUE)

## standard deviation
sd(values(tropical_land[[2]]), na.rm=TRUE)
sd(values(temperate_land[[2]]), na.rm=TRUE)
sd(values(tropical_ocean[[2]]), na.rm=TRUE)
sd(values(temperate_ocean[[2]]), na.rm=TRUE)

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

## reorder factor levels
df$lat_and_realm <- factor(df$lat_and_realm, 
                           levels = c("Temperate_Land","Tropical_Land","Temperate_Ocean","Tropical_Ocean"), 
                           ordered = TRUE)
hist <- df %>%
  ggplot(., aes(x = slope_spec_exp, fill = lat_and_realm)) + 
  geom_histogram(position = position_dodge(), binwidth = 0.0002) +
  theme_light() +
  labs(x = "", y = "Frequency", fill = "Latitude and realm") + #x = "Slope of spectral exponent"
  scale_fill_manual(values = pal_lat_and_realm, 
                    labels = c("Temperate land", "Tropical land", "Temperate ocean", "Tropical ocean")) +
  coord_flip() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.title.x = element_text(size = 10),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = -0.9, unit = "cm")) 

## extract legend 
legend <- get_legend(hist)
hist <- hist + guides(fill = "none")

## map it!
trop_land_df <- data.frame(rasterToPoints(tropical_land))
trop_land_df$lat_and_realm = "Tropical_Land"
temp_land_df <- data.frame(rasterToPoints(temperate_land))
temp_land_df$lat_and_realm = "Temperate_Land"
trop_ocean_df <- data.frame(rasterToPoints(tropical_ocean))
trop_ocean_df$lat_and_realm = "Tropical_Ocean"
temp_ocean_df <- data.frame(rasterToPoints(temperate_ocean))
temp_ocean_df$lat_and_realm = "Temperate_Ocean"
map_df <- rbind(temp_land_df,trop_ocean_df, trop_land_df, temp_ocean_df)

## reorder factor levels
map_df$lat_and_realm <- factor(map_df$lat_and_realm, 
                               levels = c("Temperate_Land","Tropical_Land","Temperate_Ocean","Tropical_Ocean"), 
                               ordered = TRUE)

map_df %>%
  ggplot(., aes(x = x, y = y, fill = lat_and_realm)) + 
  geom_raster() +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = pal_lat_and_realm) +
  labs(fill = "Slope of\nspectral\nexponent") +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) 

## boxplot:

## boxplot:
box <- df %>%
  ggplot(., aes(y = slope_spec_exp, x = lat_and_realm, fill = lat_and_realm)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Slope of spectral exponent", x = "") +
  scale_fill_manual(values = pal_lat_and_realm) +
  guides(fill = "none")+
  scale_x_discrete(labels = 
                     c("Temperate land", "Tropical land", "Temperate ocean", "Tropical ocean")) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3) +
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))

## arrange plots side by side:
lay <- rbind(c(1,1,1,1,2,2))
besties <- grid.arrange(box, hist, layout_matrix = lay)

ggsave(besties, path = "figures/analysis-workflow", filename = "o-vs-l-across-lat_spec-change.png", 
       device = "png", width = 8, height = 4)

## how different is the global average spectral change across the temperate vs. tropical land vs. ocean?
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
tropical_ocean <- mask(tropical, tropical_land, inverse = TRUE)
plot(tropical_land[[1]])
plot(tropical_ocean[[1]])

temperate_land <- crop(land_mask, temperate)
temperate_land <- mask(temperate, temperate_land)
temperate_ocean <- mask(temperate, temperate_land, inverse = TRUE)
plot(temperate_land[[1]])
plot(temperate_ocean[[1]])

## convert raster layer to data frame
df_tropland <- data.frame(rasterToPoints(tropical_land))
df_templand <- data.frame(rasterToPoints(temperate_land))
df_tropland <- gather(df_tropland, key = "window_number", value = "spec_exp", c(3:ncol(df_tropland)))
df_templand <- gather(df_templand, key = "window_number", value = "spec_exp", c(3:ncol(df_templand)))

tropland_average <- df_tropland %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Tropical_Land" = mean_spec_exp) 
tropland_average$year = seq(1871, 2081, by = 10)

templand_average <- df_templand %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Temperate_Land" = mean_spec_exp) 
templand_average$year = seq(1871, 2081, by = 10)

df_tropocean <- data.frame(rasterToPoints(tropical_ocean))
df_tempocean <- data.frame(rasterToPoints(temperate_ocean))
df_tropocean <- gather(df_tropocean, key = "window_number", value = "spec_exp", c(3:ncol(df_tropocean)))
df_tempocean <- gather(df_tempocean, key = "window_number", value = "spec_exp", c(3:ncol(df_tempocean)))

tropocean_average <- df_tropocean %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Tropical_Ocean" = mean_spec_exp)
tropocean_average$year = seq(1871, 2081, by = 10)

tempocean_average <- df_tempocean %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Temperate_Ocean" = mean_spec_exp)
tempocean_average$year = seq(1871, 2081, by = 10)

average <- left_join(tropland_average, tropocean_average) %>%
  left_join(., templand_average) %>%
  left_join(., tempocean_average)

average <- gather(average, key = 'lat_and_realm', value = 'mean_spec_exp', "Tropical_Land", 
                  "Tropical_Ocean", "Temperate_Land", "Temperate_Ocean") %>%
  filter(!is.na(mean_spec_exp))

# reorder factor
average$lat_and_realm <- factor(average$lat_and_realm, 
                                levels = c("Temperate_Land","Tropical_Land","Temperate_Ocean","Tropical_Ocean"), 
                                ordered = TRUE)

average %>%
  ggplot(., aes(x = year, y = mean_spec_exp, fill = lat_and_realm, colour = lat_and_realm)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") + 
  geom_smooth(method = "lm")  +
  scale_colour_manual(values = pal_lat_and_realm) +
  scale_fill_manual(values = pal_lat_and_realm) +
  guides(fill = "none", colour = "none")

average %>%
  filter(lat_and_realm == "Temperate_Ocean") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average %>%
  filter(lat_and_realm == "Temperate_Land") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average %>%
  filter(lat_and_realm == "Tropical_Ocean") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average %>%
  filter(lat_and_realm == "Tropical_Land") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

