## making a didactic figure 
library(tidyverse)
library(raster)
library(gridExtra)
library(cowplot)
select <- dplyr::select

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######       Setting paths     #######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
path = "/Volumes/SundayLab/CMIP5-GCMs/" 

## create vector of file folders to put data into:
gcm_models <- c("01_CMCC-CESM", "02_CMCC-CM", '03_CMCC-CMS', '04_MPI-ESM-LR', '05_MPI-ESM-MR',
                "06_GFDL-ESM2G", '07_GFDL-CM3', '08_GFDL-ESM2M', '09_HadGEM2-CC', '10_HadGEM2-ES',
                "11_HadGEM2-AO", '12_IPSL-CM5A-LR', '13_IPSL-CM5B-LR', '14_MIROC5', '15_MIROC5-ESM-CHEM',
                '16_MIROC5-ESM', "17_inmcm4", '18_CNRM-CM5', "19_MRI-CGCM3", '20_MRI-ESM1',
                '21_IPSL-CM5A-MR')

folders <- paste(path, gcm_models, "/", sep = "")

path = folders[1]


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######       Making a map     #######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
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

# ## which location has fastest reddening?
# avg$lat[which(avg$l_estimate == max(avg$l_estimate))] #48.5
# avg$lon[which(avg$l_estimate == max(avg$l_estimate))] #154.5

avg$lon <- ifelse(avg$lon >= 180, avg$lon - 358, avg$lon)

countries <- map_data("world")

## add sea surface temperature:
path <- paste("/Users/nikkimoore/Documents/velocities-of-variability/data-raw/", gcm_models[1], "/", sep = "")
se_filenames <- readRDS(paste(path, "se_filenames_tos.rds",  sep = ""))

## combine all spectral exponent csvs into one big dataframe
file = 1
while (file < length(se_filenames) + 1) {
  if (file == 1) {
    spec_exp_tos <- read.csv(se_filenames[file])
  }
  else {
    spec_exp_tos <- rbind(spec_exp_tos, read.csv(se_filenames[file]))
  }
  print(paste("Reading file #", file, "/", length(se_filenames), sep = ""))
  file = file + 1
}

## reorder time_window_width so list elements are in order of increasing time window widths:
spec_exp_tos$time_window_width <- factor(spec_exp_tos$time_window_width, levels = 
                                       c("5 years", "6 years", "7 years", "8 years",
                                         "9 years", "10 years"))

avg_tos <- spec_exp_tos %>%
  filter(time_window_width == "5 years") %>%
  select(lat, lon, l_estimate, s_estimate, l_p.value, s_p.value) %>%
  unique()


## plot
l <- avg %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of spectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0) +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent",
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
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent",
               aes(x=long, y=lat, group = group)) 

s <- avg %>%
  ggplot(., aes(x = lon, y = lat, fill = s_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank()) +
  labs(x = "Longitude", y = "Latitude", fill = "Slope of spectral exponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", 
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
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent", 
               aes(x=long, y=lat, group = group)) 

ggsave(l, filename = "figures/01_slope-spec-exp_ldetrended.png", height = 3, width = 7, device = "png")
ggsave(s, filename = "figures/01_slope-spec-exp_sdetrended.png", height = 3, width = 7, device = "png")
ggsave(l_sig, filename = "figures/01_slope-spec-exp_ldetrended_sig.png", height = 3, width = 7, device = "png")
ggsave(s_sig, filename = "figures/01_slope-spec-exp_sdetrended_sig.png", height = 3, width = 7, device = "png")


## crop data to land only
avg <- select(avg, lon, lat, everything())
avg_both <- mutate(avg, lat_lon = paste(lat, lon, sep = "_")) %>%
  filter(!lat_lon %in% paste(avg_tos$lat, avg_tos$lon, sep = "_")) %>%
  select(-lat_lon)
avg_both <- rbind(avg_both, avg_tos)
## get rid of sea surface temp data that falls outside of air surface temp data
avg_both <- filter(avg_both, lon >= min(avg$lon), lat <= max(avg$lat))
raster_avg <- rasterFromXYZ(avg_both)
plot(raster_avg[[1]])

## turn back into data frame for ggplot
data <- data.frame(rasterToPoints(raster_avg))
colnames(data) <- colnames(avg)

map <- data %>%
  ggplot(., aes(x = lon, y = lat, fill = l_estimate)) + geom_raster() +
  coord_fixed() + 
  theme_minimal() + 
  theme(panel.grid = element_blank(), legend.key.size = unit(0.5, 'cm'),
        # legend.margin = margin(0,0,0,0),
        # legend.box.margin = margin(-10,-10,-10,-10),
        # plot.margin = margin(-1,0,0,0),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.text = element_text(size = 7),
        axis.ticks = element_line(colour = "grey")) +
  labs(x = "", y = "", fill = "Slope of\nspectral\nexponent") +
  scale_fill_gradient2(high = "darkblue", low = "darkred", mid = "#e7d8d3",
                       midpoint = 0, na.value = "white") +
  scale_y_continuous(breaks = c(-90, -60, -30, 0, 30, 60, 90), 
                     labels = c("90°S", "60°S", "30°S", "0°", "30°N", "60°N", "90°N"),
                     expand = c(0,0)) +
  scale_x_continuous(breaks = c(-180, -120, -60, 0, 60, 120, 180), 
                     labels = c("180°E", "120°E", "60°E", "0°", "60°W", "120°W", "180°W"),
                     expand = c(0,0)) +
  annotate("point", y = 60.5, x = 32.5, size = 2, shape = 1, colour = "darkred") +
  geom_polygon(data = countries, col="black", size = 0.1, fill = "transparent",
               aes(x=long, y=lat, group = group)) 


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######   Making panels of decreasing spectral exponent    #######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## make plot of local spectral change by showing spectral exponent as a point over time 
data$lat[which(data$l_estimate == min(data$l_estimate))] #60.5
data$lon[which(data$l_estimate == min(data$l_estimate))] #32.5

## extract the colour describing the slope of the spectral exponent at the location of choice
col <- colorRampPalette(c("darkred", "white", "darkblue"))
plot_colours <- col(23918)
positions <- data$l_estimate[order(data$l_estimate)]
which(positions == min(data$l_estimate))

spec_exp <- readRDS("data-processed/local-spectral-change_lat-60.5_lon-32.5.rds")

## plot out spectral change for 10 year window width on l detrended data
ten <- spec_exp[[6]]

spec <- ten[[3]]

## make plot of local spectral change by showing spectral exponent as a changing slope of spectrum over time
exp <- ten[[1]] %>%
  filter(term == "log10(freq)")

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
        legend.title = element_text(size = 8)) +
  annotate("text", label = "high autocorrelation", x = 1975, y = -1.58, size = 2) +
  annotate("text", label = "low autocorrelation", x = 1975, y = -1.2, size = 2) +
  scale_x_continuous(expand = c(0,10)) +
  scale_y_continuous(limits = c(-1.6, -1.18))

## extract legend 
legend <- get_legend(spec_change)
spec_change <- spec_change + guides(alpha = "none")

ggsave(spec_change, path = "figures/didactic", filename = "spec-change-red.png", device = "png",
       width = 2, height = 1.75)

ggsave(legend, path = "figures/didactic", filename = "legend.png", device = "png",
       width = 4, height = 0.5)


## filter down to include only 5 sample points using commented line
greys <- colorRampPalette(c("lightgrey", "black"))
greys <- greys(length(unique(spec$window_start_year)))

power_spec <- spec %>%
  filter(window_start_year %in% c(1871, 1911, 1971, 2051, 2071)) %>%
  ggplot(aes(x = freq, y = power, group = as.factor(window_start_year))) + 
  geom_line(aes(colour = factor(window_start_year)), alpha = 0.15) + 
  scale_colour_manual(values = greys) +
  scale_y_log10(breaks = c(0.0000001, 0.00001, 0.001, 0.1),
                                                    
                labels = c("0.0000001", "0.00001", "0.001", "0.1")) + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  geom_vline(xintercept = 1/3650, linetype = "dotted") +
  geom_vline(xintercept = 1/365, linetype = "dotted") +
  geom_vline(xintercept = 1/31, linetype = "dotted") +
  annotate(geom = "text", label = "5 years", x = 1/2800, y = 0.000005, size = 2,
           angle = 90) +
  annotate(geom = "text", label = "1 year", x = 1/280, y = 0.000005, size = 2,
           angle = 90) +
  annotate(geom = "text", label = "1 month", x = 1/25, y = 0.000005, size = 2,
           angle = 90) +
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

ggsave(power_spec, path = "figures/didactic", filename = "power-spec-red.png", device = "png",
       width = 2, height = 1.75)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######   Making panels of increasing spectral exponent    #######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## get max 
data$lat[which(data$l_estimate == max(data$l_estimate))] #49.5
data$lon[which(data$l_estimate == max(data$l_estimate))] #154.5

spec_exp <- readRDS("data-processed/local-spectral-change_lat-48.5_lon-154.5.rds")

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
        legend.title = element_text(size = 8)) +
  annotate("text", label = "low autocorrelation", x = 1975, y = -1.02, size = 2) +
  annotate("text", label = "high autocorrelation", x = 1975, y = -1.58, size = 2) +
  scale_y_continuous(limits = c(-1.6, -1)) +
  scale_x_continuous(expand = c(0,10))

## extract legend 
legend <- get_legend(spec_change_positive)
spec_change_positive <- spec_change_positive + guides(alpha = "none")

ggsave(spec_change_positive, path = "figures/didactic", filename = "spec-change-blue.png", device = "png",
       width = 2, height = 1.75)

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
  geom_vline(xintercept = 1/3650, linetype = "dotted", colour = "slategrey") +
  geom_vline(xintercept = 1/365, linetype = "dotted", colour = "slategrey") +
  geom_vline(xintercept = 1/31, linetype = "dotted", colour = "slategrey") +
  annotate(geom = "text", label = "5 years", x = 1/2800, y = 0.0000005, size = 2,
           angle = 90) +
  annotate(geom = "text", label = "1 year", x = 1/280, y = 0.0000005, size = 2,
           angle = 90) +
  annotate(geom = "text", label = "1 month", x = 1/25, y = 0.0000005, size = 2,
           angle = 90)  +
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

ggsave(power_spec_positive, path = "figures/didactic", filename = "power-spec-blue.png", device = "png",
       width = 2.1, height = 1.75)

map = map +  annotate("point", y = 48.5, x = 154.5, size = 2, shape = 1, colour = "darkblue") 

ggsave(map, path = "figures/didactic", filename = "map.png", device = "png",
       width = 8, height = 3)


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
both_directions <- grid.arrange(power_spec, spec_change, map, legend, empty, power_spec_positive, 
                                spec_change_positive,
             layout_matrix = lay)
ggsave(both_directions, path = "figures/didactic", 
       filename = "spec-exp-didactic_both-directions.png",
       height = 6,
       width = 10)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######   Making panel of time series, time series chunk and waves    #######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
### add time series 
library(lubridate)
local_ts <- readRDS("data-processed/local-time-series_lat-60.5_lon-32.5.rds") 
local_ts$date <- as.Date(ymd(local_ts$date), "%Y%m%d")

full_timeseries <- ggplot(local_ts, aes(x = date, y = s_temp)) + 
  geom_line(size = 0.1, aes(alpha = year)) +
  theme_light() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)) +
  scale_x_date(breaks = ymd(paste(seq(1871, 2091, by = 10), "-01-01", sep = "")),
               labels = c(1871, rep("", 2),
                          1901, 1911, rep("", 17), 2091),
               expand = c(0,0.1)) +
  labs(y = "") +
  guides(alpha = "none") +
  annotate("rect", xmin = ymd("1901-01-01"), xmax = ymd("1910-12-31"), 
           ymin = -25, ymax = 25,
           alpha = .1) +
  scale_y_continuous(expand = c(0,0)) 

ggsave(full_timeseries, path = "figures/didactic", filename = "time-series-full.png", device = "png",
       width = 4.5, height = 1)

timeseries <- local_ts %>%
  filter(year %in% 1901:1910) %>%
  ggplot(., aes(x = date, y = s_temp)) + 
  geom_line(size = 0.1) +
  theme_light() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)) +
  scale_x_date(breaks = c(ymd("1901-01-01", "1910-12-31")),
               labels = c("1901", "1911"),
               limits = c(ymd("1901-01-01", "1911-01-15")),
               expand = c(0,10)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "")  +
  annotate("rect", xmin = ymd("1901-01-01"), xmax = ymd("1910-12-31"), 
           ymin = -25, ymax = 25,
           alpha = .1) 

ggsave(timeseries, path = "figures/didactic", filename = "time-series-chunk.png", device = "png",
       width = 3, height = 1)

## generate waves of 1 week, month, year, 5 year frequency 
dates <- local_ts %>%
  filter(year %in% 1901:1910)
dates <- unique(dates$date)

x_seq = seq(0,3650,length.out=3650)

x <- dates
y <- 8*sin(x_seq*(2*pi/(3650/2)))
df <- data.frame(x = x, y = y)

x2 <- dates
y2 <- 4*sin(x_seq*(2*pi/365))
df2 <- data.frame(x = x2, y = y2)

x3 <- dates
y3 <- 2*sin(x_seq*(2*pi/31))
df3 <- data.frame(x = x3, y = y3)

p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + theme_void()
p <- p + geom_line(aes(x = x, y = y), data = df, color = "black",size = 0.75) +
  geom_line(aes(x = x, y = y), data = df2, color = "black", size = 0.75) +
  geom_line(aes(x = x, y = y), data = df3, color = "black", size = 0.35) 

ggsave(p, path = "figures/didactic", filename = "waves.png", device = "png",
       width = 2, height = 1)

local_ts %>%
  filter(year %in% 1901:1910) %>%
  ggplot(., aes(x = date, y = s_temp)) + 
  geom_line(size = 0.1) +
  theme_light() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank()) +
  scale_x_date(breaks = c(ymd("1901-01-01", "1910-12-31")),
               labels = c("1901", "1911")) +
  scale_y_continuous(expand = c(0,0)) +
  labs(y = "")  +
  annotate("rect", xmin = ymd("1901-01-01"), xmax = ymd("1910-12-31"), 
           ymin = -25, ymax = 25,
           alpha = .1) + 
  geom_line(aes(x = x, y = y), data = df, color = "slategrey",size = 1) +
  geom_line(aes(x = x, y = y), data = df2, color = "slategrey", size = 1) +
  geom_line(aes(x = x, y = y), data = df3, color = "slategrey", size = 0.5) 


dates <- local_ts %>%
  filter(year %in% 1901:1910) 
dates <- unique(dates$date)

x_seq = seq(0,3650,length.out=3650)

x <- dates
y <- 8*sin(x_seq*(2*pi/(3650/2))) - 25
df <- data.frame(x = x, y = y)

x2 <- dates
y2 <- 4*sin(x_seq*(2*pi/365)) - 10
df2 <- data.frame(x = x2, y = y2)

x3 <- dates
y3 <- 2*sin(x_seq*(2*pi/31)) 
df3 <- data.frame(x = x3, y = y3)

squiggles <- ggplot(data = data.frame(x = 0), mapping = aes(x = x)) + 
  theme_light() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y.left = element_blank(),
        axis.ticks.x.bottom = element_blank()) +
  geom_line(aes(x = x, y = y), data = df, color = "black",size = 0.25) +
  geom_line(aes(x = x, y = y), data = df2, color = "black", size = 0.25) +
  geom_line(aes(x = x, y = y), data = df3, color = "black", size = 0.25) +
  annotate("text", label = " 5 years", y = mean(y), 
           x = ymd("1912-01-30"), size = 2.5) +
  annotate("text", label = "1 year ", y = mean(y2), 
           x = ymd("1912-01-30"), size = 2.5) +
  annotate("text", label = " 1 month", y = mean(y3), 
           x = ymd("1912-01-30"), size = 2.5) +
  scale_x_date(limits = c(ymd("1900-01-01", "1918-12-31")),
               expand = c(0,0.1))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######   Combining the plots     #######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## time series - full
## time series - subset 
## squiggles
## rest of didactic 

## combine plots:
lay <- rbind(c(11, rep(8,8), 5),
             c(11, rep(8,8), 5),
             c(11, rep(9,9)),
             c(11, rep(9,9)),
             c(5, rep(10,9)),
             c(5, rep(10,9)),
             c(5,4,4,4,4,4,4,4,4,5),
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

axis <- ggplot(data = spec, aes(x = freq, y = power))+
  theme_void() +
  scale_x_continuous(limits = c(0.1,0.26)) +
  annotate("text",  label = "Detrended temperature (°C)", x = 0.25, y = 0, 
           angle = 90, size = 3)

ddtc <- grid.arrange(power_spec, spec_change, map,
            legend, empty, 
            power_spec_positive, spec_change_positive,
            full_timeseries, timeseries, squiggles, axis,
            layout_matrix = lay)
ggsave(ddtc, path = "figures/didactic", 
       filename = "spec-exp-didactic_both-directions_with-time-series.png",
       height = 9,
       width = 9)


  

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
######   Making panel of high versus low autocorrelation     #######
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## read in steep time series
spec_exp_steep <- readRDS("data-processed/local-spectral-change_lat--70.5_lon--141.5.rds")

## plot out spectral change for 10 year window width on l detrended data
ten <- spec_exp_steep[[6]]

spec_steep <- ten[[3]]

one <- filter(spec, window_start_year == 1871)
two <- filter(spec_steep, window_start_year == 1991)

two_examples <- rbind(one, two)
  
two_ex_plot <- two_examples %>%
  ggplot(aes(x = freq, y = power, group = as.factor(window_start_year))) + 
  geom_line(colour = "lightgrey") +  
  scale_y_log10(breaks = c(0.0000001, 0.00001, 0.001, 0.1),
                labels = c("0.0000001", "0.00001", "0.001", "0.1")) + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  geom_vline(xintercept = 1/3650, linetype = "dotted") +
  geom_vline(xintercept = 1/365, linetype = "dotted") +
  geom_vline(xintercept = 1/31, linetype = "dotted") +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1, colour = "black") + 
  theme_light() +
  labs(x = "Frequency", y = "Power") + 
  guides(colour = "none") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

ggsave(two_ex_plot, path = "figures/didactic", filename = "high-vs-low.png", device = "png",
       width = 2.5, height = 2.5)

## make plot with only one
one_ex_plot <- one %>%
  ggplot(aes(x = freq, y = power)) + 
  geom_line() +  
  scale_y_log10(breaks = c(0.0000001, 0.00001, 0.001, 0.1),
                labels = c("0.0000001", "0.00001", "0.001", "0.1")) + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  geom_vline(xintercept = 1/3650, linetype = "dotted") +
  geom_vline(xintercept = 1/365, linetype = "dotted") +
  geom_vline(xintercept = 1/31, linetype = "dotted") +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1, colour = "black") + 
  theme_light() +
  labs(x = "Frequency", y = "Power") + 
  guides(colour = "none") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

ggsave(one_ex_plot, path = "figures/didactic", filename = "red-noise.png", device = "png",
       width = 4.5, height = 3.5)

one_ex_plot <- one %>%
  ggplot(aes(x = freq, y = power)) + 
  geom_line(colour = "transparent") +  
  scale_y_log10(breaks = c(0.0000001, 0.00001, 0.001, 0.1),
                labels = c("0.0000001", "0.00001", "0.001", "0.1")) + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  geom_vline(xintercept = 1/3650, linetype = "dotted") +
  geom_vline(xintercept = 1/365, linetype = "dotted") +
  geom_vline(xintercept = 1/31, linetype = "dotted") +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1, colour = "black") + 
  theme_light() +
  labs(x = "Frequency", y = "Power") + 
  guides(colour = "none") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        rect = element_rect(fill = "transparent"),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))


ggsave(one_ex_plot, path = "figures/didactic", filename = "red-noise_slope.png", device = "png",
       width = 4.5, height = 3.5)

lm(data = one, power ~ freq)

guides = one %>%
  ggplot(aes(x = freq, y = power)) + 
  geom_line(colour = "darkgrey") +  
  scale_y_log10(breaks = c(0.0000001, 0.00001, 0.001, 0.1),
                labels = c("0.0000001", "0.00001", "0.001", "0.1")) + 
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  geom_vline(xintercept = 1/3650, linetype = "dotted") +
  geom_vline(xintercept = 1/365, linetype = "dotted") +
  geom_vline(xintercept = 1/31, linetype = "dotted") +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1, colour = "black") + 
  theme_light() +
  labs(x = "Frequency", y = "Power") + 
  guides(colour = "none") + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        axis.line.y.left = element_line(size=0.5),
        axis.line.x.bottom = element_line(size=0.5),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6)) + 
  geom_abline(intercept = 0.9, slope = 0) + 
  geom_abline(intercept = log(0.07), slope = -1) + 
  geom_abline(intercept = log(0.002), slope = -2) + 
  geom_abline(intercept = log(0.000055), slope = -3) 

ggsave(guides, path = "figures/didactic", filename = "guides.png", device = "png",
       width = 4.5, height = 3.5)
  