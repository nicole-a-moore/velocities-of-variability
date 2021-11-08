## making a didactic figure 
library(tidyverse)
library(raster)
library(gridExtra)
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

map <- data %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), legend.key.size = unit(0.5, 'cm'),
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,-10,-10,-10),
        plot.margin = margin(-1,0,0,0),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  labs(x = "", y = "", fill = "Slope of\nspectral\nexponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white",
                       breaks = c(-0.001, -0.0005, 0, 0.0005), 
                     labels =  c("-0.001", "-0.0005", "0", "0.0005")) +
  scale_x_continuous(breaks = c(), labels = c("")) +
  scale_y_continuous(breaks = c(), labels = c("")) +
  annotate("point", y = 81.5, x = -94.5, size = 2, shape = 1, colour = "darkred") 
  


## make plot of local spectral change by showing spectral exponent as a point over time 
data$lat[which(data$l_estimate == min(data$l_estimate))] #81.5
data$lon[which(data$l_estimate == min(data$l_estimate))] #94.5

## extract the colour describing the slope of the spectral exponent at the location of choice
col <- colorRampPalette(c("darkred", "white", "darkblue"))
plot_colours <- col(23918)
positions <- data$l_estimate[order(data$l_estimate)]
which(positions == min(data$l_estimate))

spec_exp <- readRDS("data-processed/local-spectral-change_lat-81.5_lon-94.5.rds")

## plot out spectral change for 10 year window width on l detrended data
ten <- spec_exp[[6]]

spec <- ten[[3]]

## make plot of local spectral change by showing spectral exponent as a changing slope of spectrum over time
exp <- ten[[1]] %>%
  filter(term == "log10(freq)")

spec_change <- exp %>%
  ggplot(aes(x = window_start_year, y = estimate, colour = window_start_year)) + 
  geom_pointrange(aes(ymin = estimate-std.error, ymax = estimate+std.error, 
                      alpha = window_start_year), colour = "darkred") +
  theme_light() +
  geom_smooth(method = "lm", alpha = 0.3, se = F, colour = "black") +
  scale_alpha(range = c(0.01, 0.9)) +
  labs(x = "Time window start year", y = "Spectral exponent",
       alpha = "Time window start year") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        legend.position="top",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8))

## extract legend 
legend <- get_legend(spec_change)
spec_change <- spec_change + guides(alpha = "none")

## filter down to include only 5 sample points using commented line
power_spec <- spec %>%
  filter(window_start_year %in% c(1871, 1921, 1971, 2041, 2071)) %>%
  ggplot(aes(x = freq, y = power, group = window_start_year)) + 
  geom_line(aes(alpha = window_start_year)) + scale_alpha(range = c(0.01, 0.3)) +
  scale_y_log10(breaks = c(0.0000001, 0.00001, 0.001, 0.1),
                labels = c("0.0000001", "0.00001", "0.001", "0.1")) + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1.5,
            aes(alpha = window_start_year),
            colour = "darkred") + scale_alpha(range = c(0.01, 0.9)) +
  theme_light() +
  labs(x = "Frequency", y = "Power") + 
  guides(alpha = "none") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

## combine plots:
lay <- rbind(c(5,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5),
             c(5,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,5),
             c(5,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,5),
             c(5,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,5),
             c(5,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,5),
             c(5,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,5),
             c(5,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,5),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5),
             c(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,5))
empty <- ggplot(spec) + geom_blank() + theme_void()
grid.arrange(power_spec, spec_change, map, legend, empty,
             layout_matrix = lay)

## no but those colours are misleading
spec_change <- exp %>%
  ggplot(aes(x = window_start_year, y = estimate, colour = window_start_year)) + 
  geom_pointrange(aes(ymin = estimate-std.error, ymax = estimate+std.error, 
                      alpha = window_start_year), colour = "black") +
  theme_light() +
  geom_smooth(method = "lm", alpha = 0.3, se = F, colour = "darkred") +
  scale_alpha(range = c(0.01, 0.9)) +
  labs(x = "Time window start year", y = "Spectral exponent",
       alpha = "Time window start year") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        legend.position="top",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8))

## extract legend 
legend <- get_legend(spec_change)
spec_change <- spec_change + guides(alpha = "none")

## filter down to include only 5 sample points using commented line
greys <- colorRampPalette(c("lightgrey", "black"))
greys <- greys(length(unique(spec$window_start_year)))

power_spec <- spec %>%
  filter(window_start_year %in% c(1871, 1921, 1971, 2041, 2071)) %>%
  ggplot(aes(x = freq, y = power, group = as.factor(window_start_year))) + 
  geom_line(aes(colour = factor(window_start_year)), alpha = 0.15) + 
  scale_colour_manual(values = greys) +
  scale_y_log10(breaks = c(0.0000001, 0.00001, 0.001, 0.1),
                labels = c("0.0000001", "0.00001", "0.001", "0.1")) + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1,
            aes(alpha = window_start_year),
            colour = "black") + scale_alpha(range = c(0.01, 0.9)) +
  theme_light() +
  labs(x = "Frequency", y = "Power") + 
  guides(alpha = "none", colour = "none") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

negative_slope_only <- grid.arrange(power_spec, spec_change, map, legend, empty,
             layout_matrix = lay)
ggsave(negative_slope_only, path = "figures/didactic", 
       filename = "spec-exp-didactic_negative-only.png",
       height = 4,
       width = 6)



## this example is fast change - colour accordingly and make example of slow change/opposite slope?
## get max 
data$lat[which(data$l_estimate == max(data$l_estimate))] #-11.5
data$lon[which(data$l_estimate == max(data$l_estimate))] #153.5

spec_exp <- readRDS("data-processed/local-spectral-change_lat--11.5_lon-153.5.rds")

## plot out spectral change for 10 year window width on l detrended data
ten <- spec_exp[[6]]

spec <- ten[[3]]

## make plot of local spectral change by showing spectral exponent as a changing slope of spectrum over time
exp <- ten[[1]] %>%
  filter(term == "log10(freq)")

spec_change_positive <- exp %>%
  ggplot(aes(x = window_start_year, y = estimate, colour = window_start_year)) + 
  geom_pointrange(aes(ymin = estimate-std.error, ymax = estimate+std.error, 
                      alpha = window_start_year), colour = "black") +
  theme_light() +
  geom_smooth(method = "lm", alpha = 0.3, se = F, colour = "darkblue") +
  scale_alpha(range = c(0.01, 0.9)) +
  labs(x = "Time window start year", y = "Spectral exponent",
       alpha = "Time window start year") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        legend.position="top",
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 8)) 

## extract legend 
legend <- get_legend(spec_change_positive)
spec_change_positive <- spec_change_positive + guides(alpha = "none")

## filter down to include only 5 sample points using commented line
greys <- colorRampPalette(c("lightgrey", "black"))
greys <- greys(length(unique(spec$window_start_year)))

power_spec_positive <- spec %>%
  filter(window_start_year %in% c(1871, 1921, 1971, 2041, 2071)) %>%
  ggplot(aes(x = freq, y = power, group = as.factor(window_start_year))) + 
  geom_line(aes(colour = factor(window_start_year)), alpha = 0.15) + 
  scale_colour_manual(values = greys) +
  scale_y_log10(breaks = c(0.0000001, 0.00001, 0.001, 0.1),
                labels = c("0.0000001", "0.00001", "0.001", "0.1")) + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1,
            aes(alpha = window_start_year),
            colour = "black") + scale_alpha(range = c(0.01, 0.9)) +
  theme_light() +
  labs(x = "Frequency", y = "Power") + 
  guides(alpha = "none", colour = "none") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)) 

map = map +  annotate("point", y = -11.5, x = 153.5, size = 2, shape = 1, colour = "darkblue") 


## combine plots:
lay <- rbind(c(5,4,4,4,4,4,4,4,4,5),
             c(5,1,1,2,2,6,6,7,7,5),
             c(5,1,1,2,2,6,6,7,7,5),
             c(5,1,1,2,2,6,6,7,7,5),
             c(5,3,3,3,3,3,3,3,3,5),
             c(5,3,3,3,3,3,3,3,3,5),
             c(5,3,3,3,3,3,3,3,3,5), 
             c(5,3,3,3,3,3,3,3,3,5),
             c(5,3,3,3,3,3,3,3,3,5), 
             c(5,3,3,3,3,3,3,3,3,5))
empty <- ggplot(spec) + geom_blank() + theme_void()
both_directions <- grid.arrange(power_spec, spec_change, map, legend, empty, power_spec_positive, spec_change_positive,
             layout_matrix = lay)
ggsave(both_directions, path = "figures/didactic", 
       filename = "spec-exp-didactic_both-directions.png",
       height = 6,
       width = 10)
