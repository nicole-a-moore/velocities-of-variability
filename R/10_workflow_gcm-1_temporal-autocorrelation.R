## whole workflow for first GCM  !!!
library(tidyverse)
library(raster)
library(PNWColors)
library(gridExtra)
library(cowplot)
select <- dplyr::select

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
s_stack_tas <- readRDS("/Volumes/SundayLab/CMIP5-GCMs/01_CMCC-CESM/s_stack_list.rds")[[6]] 
extent(s_stack_tas) <- c(-179, 179, -89, 89)
s_stack_tos <- readRDS("data-raw/01_CMCC-CESM/s_stack_list_tos.rds")[[6]] 

s_stack_tas <- crop(s_stack_tas, s_stack_tos)
s_stack_tos <- crop(s_stack_tos, s_stack_tas)

## mosaic rasters together 
s_stack_tas <- mask(s_stack_tas, land)
plot(s_stack_tas[[1]])

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

hist <- df %>%
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
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = -0.9, unit = "cm")) 

map_legend <- data.frame(rasterToPoints(legend_rast)) %>%
  mutate(window_1 = as.factor(window_1)) %>%
  ggplot(., aes(x = x, y = y, fill = window_1)) +
  geom_raster() +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none")

## boxplot:
box <- df %>%
  ggplot(., aes(y = mean_spec_exp, x = realm, fill = realm)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none") +
  scale_x_discrete(labels = c("Land", "Ocean")) +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3) +
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))

## arrange plots side by side:
lay <- rbind(c(1,1,1,1,2,2))
besties <- grid.arrange(box, hist, layout_matrix = lay)

ggsave(besties, path = "figures/analysis-workflow", filename = "o-vs-l_mean-local-spec-exp.png", 
       device = "png", width = 8, height = 4)


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
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = -0.9, unit = "cm")) 

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
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"))

## arrange plots side by side:
lay <- rbind(c(1,1,1,1,2,2))
besties <- grid.arrange(box, hist, layout_matrix = lay)

ggsave(besties, path = "figures/analysis-workflow", filename = "o-vs-l_spec-change.png", 
       device = "png", width = 8, height = 4)


## how different is the global average spectral change across the ocean vs. land?
## use mosaic_specexp
mosaic_specexp <- readRDS("data-processed/01_tos-tas-mosaic_spectral-exponent.rds")

## split data into land and ocean
land <- mask(mosaic_specexp, land_mask)
plot(land[[1]])
ocean <- mask(mosaic_specexp, land_mask, inverse = T)
plot(ocean[[1]])

## convert raster layer to data frame
df_land <- data.frame(rasterToPoints(land))
df_land <- gather(df_land, key = "window_number", value = "spec_exp", c(3:ncol(df_land)))

land_average <- df_land %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Land" = mean_spec_exp) 
land_average$year = seq(1871, 2081, by = 10)

df_ocean <- data.frame(rasterToPoints(ocean))
df_ocean <- gather(df_ocean, key = "window_number", value = "spec_exp", c(3:ncol(df_ocean)))

ocean_average <- df_ocean %>%
  group_by(window_number) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_number, mean_spec_exp) %>%
  unique() %>%
  rename("Ocean" = mean_spec_exp)
ocean_average$year = seq(1871, 2081, by = 10)

average <- left_join(land_average, ocean_average)
average <- gather(average, key = 'realm', value = 'mean_spec_exp', "Land", 
                  "Ocean") %>%
  filter(!is.na(mean_spec_exp))

average %>%
  ggplot(., aes(x = year, y = mean_spec_exp, fill = realm, colour = realm)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") + 
  geom_smooth(method = "lm")  +
  scale_colour_manual(values = pal_realm) +
  scale_fill_manual(values = pal_realm) +
  guides(fill = "none", colour = "none")

average %>%
  filter(realm == "Ocean") %>% 
  lm(mean_spec_exp ~ year,
     data = .)

average %>%
  filter(realm == "Land") %>% 
  lm(mean_spec_exp ~ year,
     data = .)



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
tropical_ocean <- mask(tropical, tropical_land, inverse = T)
plot(tropical_land[[1]])
plot(tropical_ocean[[1]])

temperate_land <- crop(land_mask, temperate)
temperate_land <- mask(temperate, temperate_land)
temperate_ocean <- mask(temperate, temperate_land, inverse = T)
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

