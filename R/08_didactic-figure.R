## making a didactic figure 
library(tidyverse)
library(raster)
select <- dplyr::select

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

avg <- spec_exp %>%
  filter(time_window_width == "5 years") %>%
  select(lat, lon, l_estimate, s_estimate, l_p.value, s_p.value) %>%
  unique()

## which location has fastest reddening?
avg$lat[which(avg$l_estimate == max(avg$l_estimate))] #-7.5
avg$lon[which(avg$l_estimate == max(avg$l_estimate))] #5.5

avg$lon <- ifelse(avg$lon >= 180, avg$lon - 358, avg$lon)

countries <- map_data("world")

## plot
l <- avg %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of spectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 


## remove non-significant slopes
l_sig <- avg %>%
  mutate(l_estimate = ifelse(l_p.value < 0.05, l_estimate, NA)) %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of spectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

s <- avg %>%
  ggplot(., aes(x = lon, y = lat, fill = s_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of spectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

s_sig <- avg %>%
  mutate(s_estimate = ifelse(l_p.value < 0.05, s_estimate, NA)) %>%
  ggplot(., aes(x = lon, y = lat, fill = s_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of spectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", alpha = 0.5,
               aes(x=long, y=lat, group = group)) 

ggsave(l, filename = "figures/01_slope-spec-exp_ldetrended.png", height = 3, width = 7, device = "png")
ggsave(s, filename = "figures/01_slope-spec-exp_sdetrended.png", height = 3, width = 7, device = "png")
ggsave(l_sig, filename = "figures/01_slope-spec-exp_ldetrended_sig.png", height = 3, width = 7, device = "png")
ggsave(s_sig, filename = "figures/01_slope-spec-exp_sdetrended_sig.png", height = 3, width = 7, device = "png")


### ### ### ### ### ###
###  make didactic  ###
### ### ### ### ### ###

## crop data to land only
terr <- raster("data-processed/raster_terr_mask.nc") 
avg <- select(avg, lon, lat, everything())
raster_avg <- rasterFromXYZ(avg)
terr <- crop(terr, raster_avg)
raster_avg <- mask(raster_avg, terr)
plot(raster_avg[[1]])

## turn back into data frame for ggplot
data <- data.frame(rasterToPoints(raster_avg))
colnames(data) <- colnames(avg)

data %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of spectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") 


## make plot of local spectral change by showing spectral exponent as a point over time 
data$lat[which(data$l_estimate == min(data$l_estimate))] #81.5
data$lon[which(data$l_estimate == min(data$l_estimate))] #94.5

spec_exp <- readRDS("data-processed/local-spectral-change_lat-81.5_lon-94.5.rds")

## plot out spectral change for 10 year window width on l detrended data
ten <- spec_exp_list[[6]]

spec <- ten[[3]]

spec %>%
  ggplot(aes(x = freq, y = power, group = window_start_year)) + 
  geom_line(aes(alpha = window_start_year)) + scale_alpha(range = c(0.01, 0.3)) +
  scale_y_log10() + scale_x_log10() +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1.5,
            aes(alpha = window_start_year),
            colour = "red") + scale_alpha(range = c(0.01, 0.9)) +
  theme_light() +
  labs(x = "Frequency", y = "Power")

## make plot of loca1l spec1tral change by showing spectral exponent as a changing slope of spectrum over time
exp <- ten[[1]] %>%
  filter(term == "log10(freq)")

exp %>%
  ggplot(aes(x = window_start_year, y = estimate, colour = window_start_year)) + 
  geom_pointrange(aes(ymin = estimate-std.error, ymax = estimate+std.error), colour = "red") +
  theme_light() +
  geom_smooth(data = data.frame(x = exp$window_start_year, 
                                y = exp$estimate), 
              aes(x = x, y = y), method = "lm", alpha = 0.3)





