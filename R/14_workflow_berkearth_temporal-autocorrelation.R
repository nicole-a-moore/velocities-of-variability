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
s_stack_tas_PSD_low <- readRDS("data-processed/BerkeleyEarth/be_l_stack_list_PSD_low.rds")[[6]] 
s_stack_tas_AWC <- readRDS("data-processed/BerkeleyEarth/be_s_stack_list_AWC.rds")[[6]] 
s_stack_tas_PSD_high <- readRDS("data-processed/BerkeleyEarth/be_s_stack_list_PSD_high.rds")[[6]] 
s_stack_tas_PSD_all <- readRDS("data-processed/BerkeleyEarth/be_s_stack_list_PSD_all.rds")[[6]] 

## crop to same extent
list_tas <- c(s_stack_tas_PSD_low, s_stack_tas_AWC, s_stack_tas_PSD_high, s_stack_tas_PSD_all)

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
s_stack_tas_PSD_all <- list_tas[[4]]

mosaic_specexp_low <- s_stack_tas_PSD_low
mosaic_specexp_AWC <- s_stack_tas_AWC
mosaic_specexp_high <- s_stack_tas_PSD_high
mosaic_specexp_all <- s_stack_tas_PSD_all

saveRDS(mosaic_specexp_low, "data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_low.rds")
saveRDS(mosaic_specexp_high, "data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_high.rds")
saveRDS(mosaic_specexp_AWC, "data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_AWC.rds")
saveRDS(mosaic_specexp_all, "data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_all.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#####        spectral change       #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
### get tas data
path = "data-processed/BerkeleyEarth/" 

## calculate average change in spectral exponent for each time series, across all sliding window widths 
se_filenames <- readRDS(paste(path, "be_se_filenames.rds",  sep = ""))

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

## get rid of ones with large chunks of ts missing
mosaic <- readRDS("data-processed/BerkeleyEarth/missing-data-count.rds")
## make a mask to get rid of ones with more than 10000 missing values
gt_10000 <- mosaic
gt_10000[mosaic >= 10000] <- 1
gt_10000[mosaic < 10000] <- NA
plot(gt_10000)

## get lats and lons 
gt_10000 <- data.frame(rasterToPoints(gt_10000))

## filter spec exp to exclude cells in gt_10000
spec_exp <- mutate(spec_exp, lat_lon = paste(lat, lon))
gt_10000 <- mutate(gt_10000, lat_lon = paste(y, x))

spec_exp <- filter(spec_exp, !lat_lon %in% gt_10000$lat_lon)

tas <- spec_exp %>%
  filter(time_window_width == "10 years") %>%
  select(lon, lat, l_estimate_PSD_low, s_estimate_PSD_low, l_estimate_PSD_high, s_estimate_PSD_high,
         l_estimate_PSD_all, s_estimate_PSD_all, l_estimate_AWC, s_estimate_AWC) %>%
  unique() 

## mosaic together tas and tos 
raster_tas <- rasterFromXYZ(tas)
raster_tas <- extend(raster_tas, c(0, 360, -90, 90))

mosaic_specchange <- raster_tas
plot(mosaic_specchange$s_estimate_PSD_all)

## save:
saveRDS(mosaic_specchange, "data-processed/BerkeleyEarth/01_tas-mosaic_spectral-change_multifrac.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 1a. Analyze spectral exponent  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specexp
mosaic_specexp_low <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_low.rds")
mosaic_specexp_high <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_high.rds")
mosaic_specexp_AWC <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_AWC.rds")
mosaic_specexp <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_all.rds")
land_mask <-readRDS("data-processed/01_tos-tas-land-mask.rds")

legend_rast <- land
legend_rast[is.na(land)] <- 2

## split data into land and ocean
land_low <- mosaic_specexp_low
plot(land_low[[1]])

land_high <- mask(mosaic_specexp_high, land)
plot(land_high[[1]])

land_AWC <- mask(mosaic_specexp_AWC, land)
plot(land_AWC[[1]])

land_both <- mask(mosaic_specexp, land)
plot(land_both[[1]])


## summary stats and figures:
## mean
mean(values(land_low), na.rm = TRUE)
mean(values(land_AWC), na.rm = TRUE)
mean(values(land_high), na.rm = TRUE)
mean(values(land_both), na.rm = TRUE)

## standard deviation
sd(values(land_low), na.rm = TRUE)
sd(values(land_AWC), na.rm = TRUE)
sd(values(land_high), na.rm = TRUE)
sd(values(land_both), na.rm = TRUE)

## histogram:
mean_land_low <- calc(land_low, mean)
mean_land_AWC <- calc(land_AWC, mean)
mean_land_high <- calc(land_high, mean)
mean_land_both <- calc(land_both, mean)

low_or_high = append(rep("low", length(values(mean_land_low))), 
                            rep("high", length(values(mean_land_high))))
low_or_high <- append(low_or_high, rep("both", length(values(mean_land_both))))
low_or_high <- append(low_or_high, rep("AWC", length(values(mean_land_AWC))))
values = append(values(mean_land_low), values(mean_land_high))
values = append(values,  values(mean_land_both))
values = append(values,  values(mean_land_AWC))

df <- data.frame(mean_spec_exp = values,
                 low_or_high = low_or_high)
df <- filter(df, !is.na(mean_spec_exp)) %>%
  mutate(group = paste(low_or_high, sep = ""))

hist <- df %>%
  ggplot(., aes(x = mean_spec_exp, 
                fill = low_or_high, colour = low_or_high)) + 
  geom_histogram(position = position_dodge(), binwidth = 0.05) +
  theme_light() +
  labs(x = "", y = "Frequency", 
       # fill = "Exponent", colour = "Exponent type"
       ) + #x = "Mean local spectral exponent"
  # scale_colour_discrete(labels = c("AWC method", "All frequencies",
  #                                  "High frequencies", "Low frequencies")) +
  # scale_fill_discrete(labels = c("AWC method", "All frequencies",
  #                                  "High frequencies", "Low frequencies")) +
  coord_flip() + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.title.x = element_text(size = 10),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = -0.9, unit = "cm"),
        legend.position = "none") 

## boxplot:
box <- df %>%
  ggplot(., aes(y = mean_spec_exp, x = low_or_high, fill = low_or_high)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  scale_x_discrete(labels = c("AWC method", "All frequencies",
                              "High frequencies", "Low frequencies")) +
  guides(fill = "none") +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3) +
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")) 

## arrange plots side by side:
lay <- rbind(c(1,1,1,1,2,2),
             c(1,1,1,1,2,2))
besties <- grid.arrange(box, hist, layout_matrix = lay)

ggsave(besties, path = "figures/analysis-workflow", filename = "berkeleyearth_o-vs-l_mean-local-spec-exp_mf.png", 
       device = "png", width =8, height = 4)


## correlation between low and high exponent?
land_df_low <- data.frame(rasterToPoints(mean_land_low))
land_df_low$realm = "Land"
land_df_low$low_or_high = "low"
land_df_high <- data.frame(rasterToPoints(mean_land_high))
land_df_high$realm = "Land"
land_df_high$low_or_high = "high"

## rename columns 
colnames(land_df_low)[3] = colnames(land_df_high)[3] = "mean_spec_exp"

map_df <- rbind(land_df_low, land_df_high) %>%
  mutate(group = paste(low_or_high, sep = ""))

map_df %>%
  filter(low_or_high != "both") %>%
  select(-group) %>%
  spread(key = "low_or_high", value = "mean_spec_exp") %>%
  ggplot(., aes(x = low, y = high)) + geom_point(size = 0.1) +
  geom_smooth( method = "lm", colour = "red")+
  theme_light()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 1b. Analyze spectral exponent across tropical vs. temperate latitudes #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specexp
mosaic_specexp_low <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_low.rds")
mosaic_specexp_high <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_high.rds")
mosaic_specexp_AWC <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_AWC.rds")
mosaic_specexp <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-exponent_all.rds")

## split data into tropics and temperate latitudes
tropical_low <- data.frame(rasterToPoints(mosaic_specexp_low)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)
temperate_low <- data.frame(rasterToPoints(mosaic_specexp_low)) %>%
  filter(y > 23.5 | y < -23.5) %>%
  rasterFromXYZ(.)

tropical_AWC <- data.frame(rasterToPoints(mosaic_specexp_AWC)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)
temperate_AWC <- data.frame(rasterToPoints(mosaic_specexp_AWC)) %>%
  filter(y > 23.5 | y < -23.5) %>%
  rasterFromXYZ(.)

tropical_high <- data.frame(rasterToPoints(mosaic_specexp_high)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)
temperate_high <- data.frame(rasterToPoints(mosaic_specexp_high)) %>%
  filter(y > 23.5 | y < -23.5) %>%
  rasterFromXYZ(.)

tropical_both <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y <= 23.5 & y >= -23.5) %>%
  rasterFromXYZ(.)
temperate_both <- data.frame(rasterToPoints(mosaic_specexp)) %>%
  filter(y > 23.5 | y < -23.5) %>%
  rasterFromXYZ(.)

## summary stats and figures:
## mean
mean(values(tropical_low), na.rm = TRUE)
mean(values(temperate_low), na.rm = TRUE)
mean(values(tropical_AWC), na.rm = TRUE)
mean(values(temperate_AWC), na.rm = TRUE)
mean(values(tropical_high), na.rm = TRUE)
mean(values(temperate_high), na.rm = TRUE)
mean(values(tropical_both), na.rm = TRUE)
mean(values(temperate_both), na.rm = TRUE)

## standard deviation
sd(values(tropical_low), na.rm = TRUE)
sd(values(temperate_low), na.rm = TRUE)
sd(values(tropical_AWC), na.rm = TRUE)
sd(values(temperate_AWC), na.rm = TRUE)
sd(values(tropical_high), na.rm = TRUE)
sd(values(temperate_high), na.rm = TRUE)
sd(values(tropical_both), na.rm = TRUE)
sd(values(temperate_both), na.rm = TRUE)

## histogram:
trop_list <- list(tropical_low, tropical_AWC, tropical_high,tropical_both)
temp_list <- list(temperate_low, temperate_AWC, temperate_high,temperate_both)
means_trop <- lapply(trop_list, FUN = calc, mean)
means_temp <- lapply(temp_list, FUN = calc, mean)
vals_trop <- lapply(means_trop, values)
vals_trop <- unlist(vals_trop)
vals_temp <- lapply(means_temp, values)
vals_temp <- unlist(vals_temp)

df <- data.frame(latitudinal_region = rep("tropical", length(vals_trop)),
                 mean_spec_exp = vals_trop,
                 exponent_type = rep(c("Low frequencies", "AWC method", "High frequencies", "All frequencies"),
                                     each = length(values(means_trop[[1]]))))
df <- rbind(df, data.frame(latitudinal_region = rep("temperate", length(vals_temp)),
                           mean_spec_exp = vals_temp,
                           exponent_type = rep(c("Low frequencies", "AWC method", "High frequencies", "All frequencies"),
                                               each = length(values(means_temp[[1]])))))
df <- filter(df, !is.na(mean_spec_exp))

df %>%
  ggplot(., aes(x = mean_spec_exp, fill = latitudinal_region)) + 
  geom_histogram(position = position_dodge(), binwidth = 0.05) +
  theme_light() +
  labs(x = "Mean local spectral exponent",
       y = "Frequency", fill = "Latitudinal region") +
  scale_fill_manual(values = pal_lat) +
  guides(fill = "none") +
  facet_wrap(~exponent_type)

legend_rast <- temperate_low[[5]]
legend_rast[!is.na(tropical_low[[5]])] <- 2
legend_rast[!is.na(temperate_low[[5]])] <- 1

map_legend <- data.frame(rasterToPoints(legend_rast)) %>%
  mutate(window_5 = as.factor(window_5)) %>%
  ggplot(., aes(x = x, y = y, fill = window_5)) +
  geom_raster() +
  coord_fixed() +
  theme_void() +
  scale_fill_manual(values = pal_lat) +
  guides(fill = "none")  

## boxplot:
df %>%
  ggplot(., aes(y = mean_spec_exp, x = latitudinal_region, fill = latitudinal_region)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Mean local spectral exponent", x = "") +
  scale_fill_manual(values = pal_lat) +
  guides(fill = "none") +
  scale_x_discrete(labels = c("Temperate", "Tropical")) + 
  facet_wrap(~exponent_type)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 2a. Analyze change in spectral exponent  #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specchange
mosaic_specchange_mf <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-change_multifrac.rds") 

## split data into exponent types
land_low <- mosaic_specchange_mf$s_estimate_PSD_low
plot(land_low[[1]])
land_AWC <- mosaic_specchange_mf$s_estimate_AWC
plot(land_AWC[[1]])
land_high <- mosaic_specchange_mf$s_estimate_PSD_high
plot(land_high[[1]])
land_both <- mosaic_specchange_mf$s_estimate_PSD_all
plot(land_both[[1]])

## summary stats and figures:
## mean
mean(values(land_low), na.rm = TRUE)
mean(values(land_AWC), na.rm = TRUE)
mean(values(land_high), na.rm = TRUE)
mean(values(land_both), na.rm = TRUE)

## standard deviation
sd(values(land_low), na.rm = TRUE)
sd(values(land_AWC), na.rm = TRUE)
sd(values(land_high), na.rm = TRUE)
sd(values(land_both), na.rm = TRUE)

## histogram:
low_or_high = append(rep("Low frequencies", length(values(land_low))), 
                     rep("High frequencies", length(values(land_high))))
low_or_high <- append(low_or_high, rep("All frequencies", length(values(land_both))))
low_or_high <- append(low_or_high, rep("AWC method", length(values(land_AWC))))
values <- append(values(land_low), values(land_high))
values <- append(values, values(land_both))
values <- append(values, values(land_AWC))

df <- data.frame(slope_spec_exp = values,
                 low_or_high = low_or_high)
df <- filter(df, !is.na(slope_spec_exp)) %>%
  mutate(group = paste(low_or_high, sep = ""))

df %>%
  ggplot(., aes(x = slope_spec_exp, group = group, fill = ..x..)) + 
  geom_histogram(position = position_dodge()) +
  theme_light() +
  labs(x = "Slope of spectral exponent",
       y = "Frequency",
       fill = "Slope of\nspectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  geom_vline(xintercept = 0) + 
  facet_wrap(~low_or_high)

## map it!
land_df_low <- data.frame(rasterToPoints(land_low))
land_df_low$low_or_high = "Low frequencies"
land_df_AWC <- data.frame(rasterToPoints(land_AWC))
land_df_AWC$low_or_high = "AWC method"
land_df_high <- data.frame(rasterToPoints(land_high))
land_df_high$low_or_high = "High frequencies"
land_df_both <- data.frame(rasterToPoints(land_both))
land_df_both$low_or_high = "All frequencies"

## rename columns 
colnames(land_df_low)[3] = colnames(land_df_high)[3]= colnames(land_df_both)[3] = 
  colnames(land_df_AWC)[3] = "slope_spec_exp"

map_df <- rbind(land_df_high, land_df_low) %>%
  rbind(., land_df_both) %>%
  rbind(., land_df_AWC) %>%
  mutate(group = paste(low_or_high, sep = ""))

map_df %>%
  ggplot(., aes(x = x, y = y, fill = slope_spec_exp)) + 
  geom_raster() +
  coord_fixed() +
  theme_void() +
  labs(fill = "Slope of\nspectral\nexponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) + 
  facet_wrap(~low_or_high) 

hist <- df %>%
  filter(low_or_high %in% c("Low frequencies", "High frequencies")) %>%
  ggplot(., aes(x = slope_spec_exp, colour = low_or_high, fill = low_or_high)) + 
  geom_histogram(position = position_dodge(), binwidth = 0.0002) +
  theme_light() +
  labs(x = "", y = "Frequency") + #x = "Slope of spectral exponent"
  guides(fill = "none") + coord_flip() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 6),
        panel.border = element_blank(),
        panel.grid = element_blank(), 
        axis.title.x = element_text(size = 10),
        plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = -0.9, unit = "cm"), 
        legend.position = "none") +
  scale_x_continuous(limits = c(-0.01, 0.01))

## boxplot:
box <- df %>%
  filter(low_or_high %in% c("Low frequencies", "High frequencies")) %>%
  ggplot(., aes(y = slope_spec_exp, x = low_or_high, fill = low_or_high)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Slope of spectral exponent", x = "") +
  guides(fill = "none") +
  stat_summary(fun = mean, geom = "point", shape = 5, size = 3) +
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm")) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(limits = c(-0.01, 0.01))

## arrange plots side by side:
lay <- rbind(c(1,1,1,1,2,2),
             c(1,1,1,1,2,2))
besties <- grid.arrange(box, hist, layout_matrix = lay)

ggsave(besties, path = "figures/analysis-workflow", filename = "berkeleyearth_o-vs-l_spec-change_mf.png", 
       device = "png", width = 8, height = 4)


## when high freq. exponent is increasing, is low frequency exponent decreasing?
samesign_low <- land_low 
samesign_low[land_low < 0] <- -1
samesign_low[land_low > 0] <- 1
samesign_high <- land_high 
samesign_high[land_high < 0] <- -1
samesign_high[land_high > 0] <- 1
plot(samesign_low == samesign_high)

map_df %>%
  filter(!low_or_high %in% c("All frequencies", "AWC method")) %>%
  select(-group) %>%
  spread(key = "low_or_high", value = "slope_spec_exp") %>%
  rename("low" = "Low frequencies", "high" = "High frequencies") %>%
  ggplot(., aes(x = low, y = high)) + geom_point(size = 0.1) +
  theme_light() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Slope of low-frequency exponent", y = "Slope of high-freqency exponent")


##  global average spectral change
## use mosaic_specexp
spec_exp <- read.csv("data-processed/BerkeleyEarth/spec_exp_qc.csv")

spec_exp %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", c(s_spec_exp_PSD_high, s_spec_exp_PSD_all, s_spec_exp_PSD_low,
                                                 s_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  group_by(window_start_year, exp_type) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_start_year, mean_spec_exp) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent", colour = "Exponent type") + 
  scale_colour_discrete(labels = c("AWC method", "All frequencies", "High frequencies", "Low frequencies")) +
  geom_smooth(method = "lm")  

spec_exp %>%
  filter(time_window_width == "10 years") %>%
  gather(key = "exp_type", value = "spec_exp", c(l_spec_exp_PSD_high, l_spec_exp_PSD_all, l_spec_exp_PSD_low,
                                                 l_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  group_by(window_start_year, exp_type) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window
  select(window_start_year, mean_spec_exp) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent", colour = "Exponent type") + 
  scale_colour_discrete(labels = c("AWC method", "All frequencies", "High frequencies", "Low frequencies")) +
  geom_smooth(method = "lm")

spec_exp %>%
  gather(key = "exp_type", value = "spec_exp", c(s_spec_exp_PSD_high, s_spec_exp_PSD_all, s_spec_exp_PSD_low,
                                                 s_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  group_by(window_start_year, exp_type, time_window_width) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_start_year, mean_spec_exp) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent", colour = "Exponent type") + 
  scale_colour_discrete(labels = c("AWC method", "All frequencies", "High frequencies", "Low frequencies")) +
  geom_smooth(method = "lm")  +
  facet_wrap(~time_window_width)

spec_exp %>%
  gather(key = "exp_type", value = "spec_exp", c(l_spec_exp_PSD_high, l_spec_exp_PSD_all, l_spec_exp_PSD_low,
                                                 l_spec_exp_AWC)) %>%
  filter(!is.na(spec_exp)) %>%
  group_by(window_start_year, exp_type, time_window_width) %>% ## group data by window start year
  mutate(mean_spec_exp = mean(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
  select(window_start_year, mean_spec_exp) %>%
  unique() %>%
  ggplot(., aes(x = window_start_year, y = mean_spec_exp, colour = exp_type)) + geom_point() +
  theme_light() +
  labs(x = "Time window start year", y = "Mean spectral exponent") + 
  geom_smooth(method = "lm")  +
  facet_wrap(~time_window_width)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##### 2b. Analyze change in spectral exponent across tropical vs. temperate latitudes #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## use mosaic_specchange
mosaic_specchange_mf <- readRDS("data-processed/BerkeleyEarth/01_tas-mosaic_spectral-change_multifrac.rds")

## split data into exponent types
land_low <- mosaic_specchange_mf$s_estimate_PSD_low
plot(land_low[[1]])
land_AWC <- mosaic_specchange_mf$s_estimate_AWC
plot(land_AWC[[1]])
land_high <- mosaic_specchange_mf$s_estimate_PSD_high
plot(land_high[[1]])
land_both <- mosaic_specchange_mf$s_estimate_PSD_all
plot(land_both[[1]])

## split data into tropics and temperate latitudes
list_mosaic = list(land_low, land_AWC, land_high, land_both)
list_trop <- lapply(list_mosaic, FUN = function(x) {
  data.frame(rasterToPoints(x)) %>%
    filter(y <= 23.5 & y >= -23.5) %>%
    rasterFromXYZ(.)
})
list_temp <- lapply(list_mosaic, FUN = function(x) {
  data.frame(rasterToPoints(x)) %>%
    filter(y > 23.5 | y < -23.5) %>%
    rasterFromXYZ(.)
})

## summary stats and figures:
## mean
trop_means <- lapply(list_trop, FUN = function(x) {
  mean(values(x), na.rm = TRUE)
})
temp_means <- lapply(list_temp, FUN = function(x) {
  mean(values(x), na.rm = TRUE)
})
## standard deviation
trop_sd <- lapply(list_trop, FUN = function(x) {
  sd(values(x), na.rm = TRUE)
})
temp_sd <- lapply(list_temp, FUN = function(x) {
  sd(values(x), na.rm = TRUE)
})

## histogram:
vals_trop <- lapply(list_trop, values)
vals_trop <- unlist(vals_trop)
vals_temp <- lapply(list_temp, values)
vals_temp <- unlist(vals_temp)

df <- data.frame(latitudinal_region = rep("Tropical", length(vals_trop)),
                 slope_spec_exp = vals_trop,
                 exponent_type = rep(c("Low frequencies", "AWC method", "High frequencies",
                                     "All frequencies"), each = length(values(list_trop[[1]]))))
df <- rbind(df, data.frame(latitudinal_region = rep("Temperate", length(vals_temp)),
                           slope_spec_exp = vals_temp,
                           exponent_type = rep(c("Low frequencies", "AWC method", "High frequencies",
                                  "All frequencies"), each = length(values(list_temp[[1]])))))
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
  facet_wrap(~latitudinal_region) +
  geom_vline(xintercept = 0) + facet_wrap(~exponent_type)

## map it!
trop_df <- lapply(list_trop, FUN = function(x) {
  data.frame(rasterToPoints(x))
})
trop_df <- left_join(trop_df[[1]], trop_df[[2]]) %>%
  left_join(., trop_df[[3]]) %>%
  left_join(., trop_df[[4]]) %>%
  gather(key = "low_or_high", value = estimate, -x, -y) %>%
  mutate(latitudinal_region = "Tropical")

temp_df <- lapply(list_temp, FUN = function(x) {
  data.frame(rasterToPoints(x))
})
temp_df <- left_join(temp_df[[1]], temp_df[[2]]) %>%
  left_join(., temp_df[[3]]) %>%
  left_join(., temp_df[[4]]) %>%
  gather(key = "low_or_high", value = estimate, -x, -y) %>%
  mutate(latitudinal_region = "Temperate")

map_df <- rbind(trop_df, temp_df) %>%
  mutate(exponent_type = ifelse(str_detect(low_or_high, "AWC"), "AWC method", 
                                ifelse(str_detect(low_or_high, "low"), "Low frequencies",
                                       ifelse(str_detect(low_or_high, "high"), "High frequencies", 
                                              "All frequencies"))))

map_df %>%
  ggplot(., aes(x = x, y = y, fill = estimate)) +
  geom_polygon(data = countries, col="black", size = 0.01, fill = "grey", alpha = 0.3,
               aes(x=long+180, y=lat, group = group)) +
  geom_raster() +
  coord_fixed() +
  theme_void() +
  labs(fill = "Slope of\nspectral\nexponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  facet_grid(exponent_type~latitudinal_region) +
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) 

## boxplot:
df %>%
  ggplot(., aes(y = slope_spec_exp, x = latitudinal_region, fill = latitudinal_region)) + 
  geom_boxplot() +
  theme_light() +
  labs(y = "Local slope of spectral exponent", x = "") +
  scale_fill_manual(values = pal_lat) +
  guides(fill = "none") +
  facet_wrap(~exponent_type, nrow = 2) 

