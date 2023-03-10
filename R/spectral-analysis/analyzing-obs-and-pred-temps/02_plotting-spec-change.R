## bringing all spectral exponent data together in some plots
## fitting modes to spectral exponent over time 
## and plotting predictions
library(tidyverse)
library(PNWColors)
library(ggExtra)
library(cowplot)

#################################################
###                setting paths               ## 
#################################################
## set 'path' to where you have the GCM files stored on your computer
## for me, they are here:
path = "/Volumes/NIKKI/CMIP5-GCMs/" 

## create vector of file folders to put data into:
gcm_models <- c("01_CMCC-CMS_tas",
                "02_GFDL-CM3_tas",
                "03_GFDL-ESM2G_tas",
                "04_HadGEM2-ES_tas",
                "05_inmcm4_tas",
                "06_IPSL-CM5A-MR_tas",
                "07_MIROC-ESM-CHEM_tas",
                "08_MIROC5_tas",
                "09_MPI-ESM-LR_tas",
                "10_MPI-ESM-MR_tas",
                "11_MRI-CGCM3_tas", 
                "BerkeleyEarth",
                "NOAA-OISST")

folders <- paste(path, gcm_models, "/", sep = "")


## another showing effect of AWC versus PSD, linear versus seasonally detrending on low exponent - 10 year time window


###################################################################
###     calling functions on spectral exponent files             ## 
###################################################################
## make a plot of global gcm trends + berk earth + noaa together 
num = 1
widths = c("5_years", "6_years", "7_years", "8_years", "9_years", "10_years")
while(num <= length(widths)) {
  ## read in and row bind csvs
  all <- c()
  for (i in 1:length(gcm_models)) {
    p = folders[i]
    gcm = gcm_models[i]
    
    cur <- read.csv(paste("data-processed/spectral-change-files/", gcm, "_average-se-over-time_", widths[num], ".csv", sep = ""))
    cur$gcm <- gcm
    
    all <- rbind(all, cur)
  }
  
  ## now plot them
  ## create colour palette:
  pal = c("#F3DE8A", "#3FA34D", "#4F86C6")
  
  if(num == 6) {
    plot <- all %>%
      mutate(group = paste(gcm, exponent_type)) %>%
      mutate(exponent_type = ifelse(exponent_type == "s_PSD_high", 
                                    "High-frequency spectral change (seasonally-detrended)",
                                    ifelse(exponent_type == "l_PSD_high", 
                                           "High-frequency spectral change (linearly-detrended)",
                                           ifelse(exponent_type == "s_PSD_low", 
                                                  "Low-frequency spectral change (seasonally-detrended)",
                                                  ifelse(exponent_type == "l_PSD_low", 
                                                         "Low-frequency spectral change (linearly-detrended)",
                                                         ifelse(exponent_type == "l_AWC", 
                                                                "Low-frequency spectral change\n(linearly-detrended, AWC method)",
                                                                ifelse(exponent_type == "s_AWC", 
                                                                       "Low-frequency spectral change\n(seasonally-detrended, AWC method)",
                                                                       NA))))))) %>%
      mutate(colour = ifelse(gcm == "NOAA-OISST",
                             "Observed sea surface temperature",
                             ifelse(gcm == "BerkeleyEarth", "Observed air surface temperature",
                                    "Model-predicted air surface temperature"))) %>%
      ggplot(., aes(x = year, y = mean_spec_exp, colour = colour, group = group)) +
      theme_light() +
      geom_errorbar(aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                    alpha = 0.75) +
      labs(x = "Time window start year", y = "Mean spectral exponent", colour = "Dataset:") +
      geom_smooth(method = "lm", se = FALSE) +
      geom_point() +
      facet_wrap(~exponent_type, nrow = 3) +
      scale_color_manual(values = pal) +
      theme(legend.position = "bottom")
    
    ggsave(plot, path = "figures/spectral-analysis_GCMs/", 
           filename = paste("all-datasets_", widths[num], ".png", sep = ""), device = "png",
           width = 8, height = 10)
  }
  else {
    plot <- all %>%
      mutate(group = paste(gcm, exponent_type)) %>%
      mutate(exponent_type = ifelse(exponent_type == "s_PSD_high", 
                                    "High-frequency spectral change",
                                    "Low-frequency spectral change")) %>%
      mutate(exponent_type = factor(.$exponent_type, 
                                    levels = c("Low-frequency spectral change", "High-frequency spectral change"),
                                    ordered = TRUE)) %>%
      mutate(colour = ifelse(gcm == "NOAA-OISST",
                             "Observed sea surface temperature",
                             ifelse(gcm == "BerkeleyEarth", "Observed air surface temperature",
                                    "Model-predicted air surface temperature"))) %>%
      ggplot(., aes(x = year, y = mean_spec_exp, colour = colour, group = group)) +
      theme_light() +
      geom_errorbar(aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                    alpha = 0.75) +
      labs(x = "Time window start year", y = "Mean spectral exponent", colour = "Dataset:") +
      geom_smooth(method = "lm", se = FALSE)  +
      geom_point() +
      facet_wrap(~exponent_type) +
      scale_color_manual(values = pal) +
      theme(legend.position = "bottom")
    
    ggsave(plot, path = "figures/spectral-analysis_GCMs/", 
           filename = paste("all-datasets_", widths[num], ".png", sep = ""), device = "png",
           width = 9, height = 5)
  }
  
  if(num ==1) {
    obs <- all %>%
      filter(exponent_type == "s_PSD_low") %>%
      mutate(group = gcm) %>%
      mutate(colour = ifelse(gcm == "NOAA-OISST",
                             "Observed sea surface temperature",
                             ifelse(gcm == "BerkeleyEarth", "Observed air surface temperature",
                                    "Model-predicted air surface temperature"))) %>%
      mutate(color = factor(.$colour, levels = c("Model-predicted air surface temperature", 
                                                 "Observed air surface temperature", 
                                                 "Observed sea surface temperature"),
                            ordered = TRUE)) %>%
      filter(colour %in% c( "Observed sea surface temperature", "Observed air surface temperature"))
    plot_mainfig <- all %>%
      filter(exponent_type == "s_PSD_low") %>%
      mutate(group = gcm) %>%
      mutate(colour = ifelse(gcm == "NOAA-OISST",
                             "Observed sea surface temperature",
                             ifelse(gcm == "BerkeleyEarth", "Observed air surface temperature",
                                    "Model-predicted air surface temperature"))) %>%
      mutate(color = factor(.$colour, levels = c("Model-predicted air surface temperature", 
                                                 "Observed air surface temperature", 
                                                 "Observed sea surface temperature"),
                            ordered = TRUE)) %>%
      ggplot(., aes(x = year, y = mean_spec_exp, colour = colour, group = group)) +
      theme_light()  +
      labs(x = "Year", y = "Mean spectral exponent", colour = "") +
      geom_smooth(method = "lm", se = FALSE) +
      geom_errorbar(aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                    alpha = 0.75) +
      geom_point() +
      geom_errorbar(data = obs, aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                    alpha = 0.75) +
      scale_color_manual(values = pal) +
      theme(panel.grid = element_blank())
  
    ggsave(plot_mainfig, path = "figures/spectral-analysis_GCMs/", 
           filename = paste("all-datasets_", widths[num], "_mainfig.png", sep = ""), device = "png",
           width = 10, height = 3)
    }
  
  num = num + 1
}

## read in and row bind csvs
num = 1
all <- c()
widths = c("5_years", "6_years", "7_years", "8_years", "9_years", "10_years")
while(num <= length(widths)) {

  for (i in 1:length(gcm_models)) {
    p = folders[i]
    gcm = gcm_models[i]
    
    cur <- read.csv(paste("data-processed/spectral-change-files/", gcm, "_average-se-over-time_", widths[num], ".csv", sep = ""))
    cur$gcm <- gcm
    cur$time_window_width = str_replace_all(widths[num], "\\_", " ")
    
    all <- rbind(all, cur)
  }
  
  num = num + 1
}

pal_pnw <- pnw_palette("Starfish", 5)

## now make a plot showing effect of time window width + how high and low frequency compare
plot <- all %>%
  filter(exponent_type %in% c("s_PSD_high", "s_PSD_low")) %>% ## filter to analyses on seasonally detrended data
  mutate(group = paste(gcm, exponent_type)) %>%
  mutate(time_window_width = factor(.$time_window_width, ordered = TRUE, levels = c("5 years", "6 years", "7 years", 
                                                                                  "8 years", "9 years", "10 years"))) %>%
  mutate(gcm = ifelse(gcm %in% c("BerkeleyEarth", "NOAA-OISST"), paste("_", gcm, "_", sep = ""), gcm)) %>%
  mutate(gcm = str_split_fixed(gcm, "\\_", 3)[,2]) %>%
  mutate(gcm = factor(.$gcm, levels = c("BerkeleyEarth", "NOAA-OISST", "CMCC-CMS","GFDL-CM3","GFDL-ESM2G","HadGEM2-ES","inmcm4","IPSL-CM5A-MR",
                                        "MIROC-ESM-CHEM", "MIROC5","MPI-ESM-LR",
                                        "MPI-ESM-MR","MRI-CGCM3"), ordered = TRUE)) %>%
  ggplot(., aes(x = year, y = mean_spec_exp, colour = exponent_type)) +
  theme_light() +
  geom_errorbar(aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                alpha = 0.75) +
  labs(x = "Time window start year", y = "Mean spectral exponent", colour = "") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point() +
  facet_grid(gcm~time_window_width) +
  scale_color_manual(values = pal_pnw[c(2,4)], labels = c("High-frequency spectral exponent", "Low-frequency spectral exponent")) +
  theme(legend.position = "bottom", 
        panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(size = 8), strip.text.x = element_text(size = 8))

ggsave(plot, path = "figures/spectral-analysis_GCMs/", 
       filename = paste("all-datasets_time-window-width-x-high-low.png", sep = ""), device = "png",
       width = 8, height = 12)

## and now make a plot showing effect of measurement method (PSD, AWC) and detrending data on low frequency exponent 
plot <- all %>%
  filter(time_window_width == "10 years") %>% ## filter to 10 year time windows 
  filter(!exponent_type %in% c("s_PSD_high", "l_PSD_high")) %>%  ## get rid of high frequency exponents 
  mutate(seasonal_or_linear = ifelse(str_detect(exponent_type, "s_"), "Seasonally-detrended", "Linearly-detrended"),
         wavelet_or_spectral = ifelse(str_detect(exponent_type, "PSD"), "Spectral analysis", "Wavelet analysis")) %>% ## making grouping variables 
  mutate(gcm = ifelse(gcm %in% c("BerkeleyEarth", "NOAA-OISST"), paste("_", gcm, "_", sep = ""), gcm)) %>%
  mutate(gcm = str_split_fixed(gcm, "\\_", 3)[,2]) %>%
  mutate(gcm = factor(.$gcm, levels = c("BerkeleyEarth", "NOAA-OISST", "CMCC-CMS","GFDL-CM3","GFDL-ESM2G","HadGEM2-ES","inmcm4","IPSL-CM5A-MR",
                                        "MIROC-ESM-CHEM", "MIROC5","MPI-ESM-LR",
                                        "MPI-ESM-MR","MRI-CGCM3"), ordered = TRUE)) %>%
  ggplot(., aes(x = year, y = mean_spec_exp, colour = wavelet_or_spectral)) +
  theme_light() +
  geom_errorbar(aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                alpha = 0.75) +
  labs(x = "Time window start year", y = "Mean spectral exponent", colour = "") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_point() +
  facet_grid(seasonal_or_linear~gcm) +
  scale_color_manual(values = pal_pnw[c(1,3)], 
                     #labels = c("High-frequency spectral exponent", "Low-frequency spectral exponent")
                     ) +
  theme(legend.position = "bottom", 
        panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(size = 8), strip.text.x = element_text(size = 8))

ggsave(plot, path = "figures/spectral-analysis_GCMs/", 
       filename = paste("all-datasets_seasonal-or-linear-x-wave-spec.png", sep = ""), device = "png",
       width = 12, height = 4)





## make a plot of the distribution of local trends for gcms vs obs
num = 1
widths = c("5_years", "6_years", "7_years", "8_years", "9_years", "10_years")
while(num <= length(widths)) {
  ## read in and row bind csvs
  all <- c()
  for (i in 1:length(gcm_models)) {
    p = folders[i]
    gcm = gcm_models[i]
    
    cur <- read.csv(paste("data-processed/spectral-change-files/", gcm, "_spec-change-all.csv", sep = ""))
    cur$gcm <- gcm
    
    all <- rbind(all, cur)
  }
  
  ## now plot them
  ## create colour palette:
  pal = c("#F3DE8A", "#3FA34D", "#4F86C6")
  
  if(num == 6) {
    plot <- all %>%
      filter(time_window_width == str_replace_all(widths[num], "\\_", " ")) %>%
      mutate(group = paste(gcm, exponent_type)) %>%
      mutate(exponent_type = ifelse(exponent_type == "s_PSD_high", 
                                    "High-frequency spectral change\n(seasonally-detrended)",
                                    ifelse(exponent_type == "l_PSD_high", 
                                           "High-frequency spectral change\n(linearly-detrended)",
                                           ifelse(exponent_type == "s_PSD_low", 
                                                  "Low-frequency spectral change\n(seasonally-detrended)",
                                                  ifelse(exponent_type == "l_PSD_low", 
                                                         "Low-frequency spectral change\n(linearly-detrended)",
                                                         ifelse(exponent_type == "l_AWC", 
                                                                "Low-frequency spectral change\n(linearly-detrended, AWC method)",
                                                                ifelse(exponent_type == "s_AWC", 
                                                                       "Low-frequency spectral change\n(seasonally-detrended, AWC method)",
                                                                       NA))))))) %>%
      mutate(colour = ifelse(gcm == "NOAA-OISST",
                             "Observed\nsea surface\ntemperature",
                             ifelse(gcm == "BerkeleyEarth", "Observed\nair surface\ntemperature",
                                    "Model-predicted\nair surface\ntemperature"))) %>%
      ggplot(., aes(x = colour, y = estimate, fill = colour)) + 
      geom_boxplot() +
      theme_light() +
      facet_wrap(~exponent_type) +
      scale_fill_manual(values = pal) +
      theme(legend.position = "none") +
      labs(y = "Local change in\nspectral exponent per year", fill = "", x = "")
    
    ggsave(plot, path = "figures/spectral-analysis_GCMs/", 
           filename = paste("all-datasets_", widths[num], "_boxplot.png", sep = ""), device = "png",
           width = 8, height = 3)
    
  }
  else {
    plot <- all %>%
      filter(time_window_width == str_replace_all(widths[num], "\\_", " ")) %>%
      mutate(group = paste(gcm, exponent_type)) %>%
      mutate(exponent_type = ifelse(exponent_type == "s_PSD_high", 
                                    "High-frequency spectral change",
                                    "Low-frequency spectral change")) %>%
      mutate(exponent_type = factor(.$exponent_type, 
                                    levels = c("Low-frequency spectral change", "High-frequency spectral change"),
                                    ordered = TRUE)) %>%
      mutate(colour = ifelse(gcm == "NOAA-OISST",
                             "Observed\nsea surface\ntemperature",
                             ifelse(gcm == "BerkeleyEarth", "Observed\nair surface\ntemperature",
                                    "Model-predicted\nair surface\ntemperature"))) %>%
      ggplot(., aes(x = colour, y = estimate, fill = colour)) + 
      geom_boxplot() +
      theme_light() +
      facet_wrap(~exponent_type) +
      scale_fill_manual(values = pal) +
      theme(legend.position = "none") +
      labs(y = "Local change in\nspectral exponent per year", fill = "", x = "")
    
    ggsave(plot, path = "figures/spectral-analysis_GCMs/", 
           filename = paste("all-datasets_", widths[num], "_boxplot.png", sep = ""), device = "png",
           width = 8, height = 3)
  
  }
  
  
  ## make special figure for 5-year window
  if(num == 1) {
    plot_box_mainfig <- all %>%
      filter(exponent_type == "s_PSD_low") %>%
      filter(time_window_width == str_replace_all(widths[num], "\\_", " ")) %>%
      mutate(group = gcm) %>%
      mutate(colour = ifelse(gcm == "NOAA-OISST",
                             "Observed\nsea surface\ntemperature",
                             ifelse(gcm == "BerkeleyEarth", "Observed\nair surface\ntemperature",
                                    "Model-predicted\nair surface\ntemperature"))) %>%
      group_by(group) %>%
      mutate(mean = mean(estimate), sd = sd(estimate)) %>%
      ungroup() %>%
      ggplot(., aes(x = colour, y = estimate, fill = colour, group = colour)) +  
      geom_abline(intercept = 0, slope = 0) +
      geom_boxplot(lwd = 0.3) +
      theme_light() +
      scale_fill_manual(values = pal) +
      theme(legend.position = "none", panel.grid = element_blank()) +
      labs(y = "Local change in\nspectral exponent per year", fill = "", x = "")
    
    ggsave(plot_box_mainfig, path = "figures/spectral-analysis_GCMs/", 
           filename = paste("all-datasets_", widths[num], "_boxplot-obs-and-pred.png", sep = ""), device = "png",
           width = 3.5, height = 3)
    
    ## now make one with only observed temperature
    plot_box_obs <- all %>%
      filter(gcm %in% c("BerkeleyEarth", "NOAA-OISST")) %>%
      filter(exponent_type == "s_PSD_low") %>%
      filter(time_window_width == str_replace_all(widths[num], "\\_", " ")) %>%
      mutate(group = gcm) %>%
      mutate(colour = ifelse(gcm == "NOAA-OISST",
                             "Observed\nsea surface\ntemperature",
                             ifelse(gcm == "BerkeleyEarth", "Observed\nair surface\ntemperature", NA))) %>%
      group_by(group) %>%
      mutate(mean = mean(estimate), sd = sd(estimate)) %>%
      ungroup() %>%
      ggplot(., aes(x = colour, y = estimate, fill = colour, group = colour)) +  
      geom_abline(intercept = 0, slope = 0) +
      geom_boxplot(lwd = 0.3) +
      theme_light() +
      scale_fill_manual(values = pal[c(2,3)]) +
      theme(legend.position = "none", panel.grid = element_blank()) +
      labs(y = "Local change in\nspectral exponent per year", fill = "", x = "") +
      theme(panel.background = element_rect(fill='transparent'),
            plot.background = element_rect(fill='transparent', color=NA),
            legend.background = element_rect(fill='transparent'),
            legend.box.background = element_rect(fill='transparent')) +
      scale_y_continuous(limits = c(-0.04, 0.04))
    
    ggsave(plot_box_obs, path = "figures/spectral-analysis_GCMs/", 
           filename = paste("all-datasets_", widths[num], "_boxplot-mainfig.png", sep = ""), device = "png",
           width = 2.5, height = 3, bg = "transparent")
  }
  
  num = num + 1
}

## now make maps :-)
num = 1
widths = c("5_years", "6_years", "7_years", "8_years", "9_years", "10_years")
while(num <= length(widths)) {
  
  ## read in and row bind csvs for all climate datasets
  all <- c()
  for (i in 1:length(gcm_models)) {
    p = folders[i]
    gcm = gcm_models[i]
    
    cur <- read.csv(paste("data-processed/spectral-change-files/", gcm, "_spec-change-all.csv", sep = ""))
    
    cur <- cur %>%
      filter(exponent_type == "s_PSD_low") %>%
      select(lat, lon, estimate, time_window_width) %>% ## select only low seasonal spectral exponent slope 
      filter(time_window_width == str_replace_all(widths[[num]], "\\_", " ")) %>%
      distinct()
    
    cur$gcm <- gcm
    
    all <- rbind(all, cur)
  }
  
  ## now plot them
  all_maps <- all %>%
    mutate(gcm = ifelse(gcm %in% c("BerkeleyEarth", "NOAA-OISST"), paste("_", gcm, "_", sep = ""), gcm)) %>%
    mutate(gcm = str_split_fixed(gcm, "\\_", 3)[,2]) %>%
    mutate(gcm = factor(.$gcm, levels = c("BerkeleyEarth", "NOAA-OISST", "CMCC-CMS","GFDL-CM3","GFDL-ESM2G","HadGEM2-ES","inmcm4","IPSL-CM5A-MR",
                                          "MIROC-ESM-CHEM", "MIROC5","MPI-ESM-LR",
                                          "MPI-ESM-MR","MRI-CGCM3"), ordered = TRUE)) %>%
    ggplot(., aes(x = lon, y = lat, fill = estimate)) +
    geom_raster() +
    coord_fixed() +
    theme_void() +
    facet_wrap(~gcm) +
    scale_fill_gradient2(low = "#7A81CD", high = "#FF5C7A", mid = "#FFFFFF",
                         midpoint = 0, limits = c(-0.04, 0.04), na.value = "grey") +
    labs(fill = "Change in\nspectral exponent") +
    theme(panel.background = element_rect(fill = "#ECECEC", colour = "#ECECEC"))
  
  ggsave(all_maps, path = "figures/spectral-analysis_GCMs/", 
         filename = paste("all-datasets_", widths[num], "_maps.png", sep = ""), device = "png",
         width = 8, height = 4)
  
  ## plot average across gcms 
  all %>%
    mutate(group = ifelse(gcm == "NOAA-OISST",
                           "Observed\nsea surface\ntemperature",
                           ifelse(gcm == "BerkeleyEarth", "Observed\nair surface\ntemperature",
                                  "Model-predicted\nair surface\ntemperature"))) %>%
    group_by(lat, lon, group) %>%
    mutate(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
    ungroup() %>%
    select(lat, lon, mean_estimate, group) %>% 
    distinct() %>%
    ggplot(., aes(x = lon, y = lat, fill = mean_estimate)) +
    geom_raster() +
    coord_fixed() +
    theme_void() +
    facet_wrap(~group) +
    scale_fill_gradient2(low = "#7A81CD", high = "#FF5C7A", mid = "#FFFFFF",
                         midpoint = 0, limits = c(-0.04, 0.04), na.value = "grey") +
    labs(fill = "Change in\nspectral exponent\n per year")
  
  ## now make one that combines Berk and Noaa data 
  map_comb <- all %>%
    mutate(group = ifelse(gcm %in% c("BerkeleyEarth", "NOAA-OISST"),
                          "Observed temperature",
                         "Model-predicted temperature")) %>%
    group_by(lat, lon, group) %>%
    mutate(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
    ungroup() %>%
    select(lat, lon, mean_estimate, group) %>% 
    distinct() %>%
    ggplot(., aes(x = lon, y = lat, fill = mean_estimate)) +
    geom_raster() +
    coord_fixed() +
    theme_void() +
    facet_wrap(~group) +
    scale_fill_gradient2(low = "#7A81CD", high = "#FF5C7A", mid = "#FFFFFF",
                         midpoint = 0, limits = c(-0.04, 0.04), na.value = "grey") +
    labs(fill = "Change in\nspectral exponent\nper year")
  
  ggsave(map_comb, path = "figures/spectral-analysis_GCMs/", 
         filename = paste("all-datasets_", widths[num], "_maps-obs-vs-pred.png", sep = ""), device = "png",
         width = 8, height = 2)
  
  ## now make one with only observed temperature 
  data <- all %>%
    filter(gcm %in% c("BerkeleyEarth", "NOAA-OISST")) %>%
    group_by(lat, lon) %>%
    mutate(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
    ungroup() %>%
    select(lat, lon, mean_estimate) %>% 
    distinct() 
  
  main <- data %>%
    ggplot(., aes(x = lon, y = lat, fill = mean_estimate)) +
    geom_raster() +
    coord_fixed() +
    theme_minimal() +
    # scale_fill_gradient2(low = "#03045E", high = "#8B2635", mid = "#e7d8d3",
    #                      midpoint = 0) +
    # scale_fill_gradient2(low = "#979CD8", high = "#FF859B", mid = "#FFFFFF",
    #                      midpoint = 0, limits = c(-0.04, 0.04), na.value = "grey") +
    scale_fill_gradient2(low = "#7A81CD", high = "#FF5C7A", mid = "#FFFFFF",
                         midpoint = 0, limits = c(-0.04, 0.04), na.value = "grey") +
    labs(fill = "Change in\nspectral exponent\nper year") +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0)),
                       breaks = c(60, 120, 180, 240, 300),
                       labels = c( "-120°", "-60°", "0°", "60°", "120°"))+
    scale_y_continuous(expand = c(0,0),
                       breaks = c( -60, -30, 0, 30, 60),
                       labels = c( "-60°", "-30°", "0°", "30°", "60°")) + 
    labs(x = "Longitude", y = "Latitude") +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 8),
          axis.ticks = element_line(linewidth = 0.25))
    
  
  
  mean_data_lat <- all %>%
    filter(gcm %in% c("BerkeleyEarth", "NOAA-OISST")) %>%
    group_by(lat, gcm) %>%
    mutate(mean_estimate_lat = mean(estimate, na.rm = TRUE), 
           sd_estimate_lat = sd(estimate, na.rm = TRUE)) %>%
    ungroup() %>%
    select(lat, mean_estimate_lat, sd_estimate_lat, gcm) %>% 
    distinct() %>%
    arrange(lat)
  mean_data_lon <- all %>%
    filter(gcm %in% c("BerkeleyEarth", "NOAA-OISST")) %>%
    group_by(lon, gcm) %>%
    mutate(mean_estimate_lon = mean(estimate, na.rm = TRUE),
           sd_estimate_lon = sd(estimate, na.rm = TRUE)) %>%
    ungroup() %>%
    select(lon, mean_estimate_lon, sd_estimate_lon, gcm) %>% 
    distinct() 
  
  xdens <- axis_canvas(main, axis = "x") +
    geom_ribbon(data = mean_data_lon, aes(x = lon, ymax = mean_estimate_lon+sd_estimate_lon,
                                             ymin = mean_estimate_lon-sd_estimate_lon,fill = gcm),
               alpha = 0.3) +
    geom_path(data = mean_data_lon, aes(x = lon, y = mean_estimate_lon,
                                         colour = gcm)) +
    scale_colour_manual(values = c("#3FA34D", "#4F86C6"))  +
    scale_fill_manual(values = c("#3FA34D", "#4F86C6"))  +
    geom_abline(intercept = 0, slope = 0, linewidth = 0.25) +
    theme(panel.grid = element_blank()) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0)))
  
  ydens <- axis_canvas(main, axis = "y", coord_flip = TRUE) +
    geom_ribbon(data = mean_data_lat, aes(x = lat, ymax = mean_estimate_lat+sd_estimate_lat,
                                             ymin = mean_estimate_lat-sd_estimate_lat,
                                             fill = gcm),
                   alpha = 0.3) +
    geom_path(data = mean_data_lat, aes(x = lat, y = mean_estimate_lat,
                                        colour = gcm)) +
    coord_flip() +
    scale_colour_manual(values = c("#3FA34D", "#4F86C6")) +
    scale_fill_manual(values = c("#3FA34D", "#4F86C6"))  +
    geom_abline(intercept = 0, slope = 0, linewidth = 0.25) + 
    theme(panel.grid = element_blank()) 
  
  p1 <- insert_xaxis_grob(main, xdens, grid::unit(.1, "null"), position = "top")
  p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.1, "null"), position = "right")
  ggdraw(p2)
  
  ggsave(p2, path = "figures/spectral-analysis_GCMs/", 
         filename = paste("all-datasets_", widths[num], "_maps-obs-only-with-means.png", sep = ""), device = "png",
         width = 9, height = 6)
  
  ## now make one that asks whether direction of change in observed temperatures matches direction of change in predicted temperatures
  map_mismatch <- all %>%
    mutate(group = ifelse(gcm %in% c("BerkeleyEarth", "NOAA-OISST"),
                          "Observed temperature",
                          "Model-predicted temperature")) %>%
    group_by(lat, lon, group) %>%
    mutate(mean_estimate = mean(estimate, na.rm = TRUE)) %>%
    ungroup() %>%
    select(lat, lon, group, mean_estimate) %>%
    distinct() %>%
    spread(key = group, value = mean_estimate) %>%
    mutate(match = ifelse(`Model-predicted temperature` >= 0 & `Observed temperature` >= 0, 
                          "Match",
                          ifelse(`Model-predicted temperature` < 0 & `Observed temperature` < 0, 
                                 "Match", 
                                 "Mismatch"))) %>%
    filter(!is.na(match)) %>% # %>% group_by(match) %>% tally()
    ggplot(., aes(x = lon, y = lat, fill = match)) +
    geom_raster() +
    coord_fixed() +
    theme_void() +
    labs(fill = "")
 
  
  ggsave(map_mismatch, path = "figures/spectral-analysis_GCMs/", 
         filename = paste("all-datasets_", widths[num], "_maps-mismatch.png", sep = ""), device = "png",
         width = 8, height = 2)
  
  num = num + 1
}







# ## plot to see which of missing pixels are ones with missing data 
# mosaic_noaa <- readRDS(paste( "/Volumes/NIKKI/CMIP5-GCMs/NOAA-OISST/missing-data-count.rds", sep = ""))
# mosaic_berk <- readRDS(paste( "/Volumes/NIKKI/CMIP5-GCMs/BerkeleyEarth/missing-data-count.rds", sep = ""))
# 
# gt_10000 <- mosaic_berk
# gt_10000[mosaic_berk >= 10000] <- 1
# gt_10000[mosaic_berk < 10000] <- NA
# 
# mosaic_berk <- gt_10000
# 
# gt_10000 <- mosaic_noaa
# gt_10000[mosaic_noaa >= 10000] <- 1
# gt_10000[mosaic_noaa < 10000] <- NA
# 
# mosaic_noaa <- gt_10000
# mosaic_noaa <- data.frame(rasterToPoints(mosaic_noaa))
# mosaic_noaa$x <- ifelse(mosaic_noaa$x >= 179.5, mosaic_noaa$x - 179, mosaic_noaa$x + 180)
# mosaic_noaa <- rasterFromXYZ(mosaic_noaa)
# #mosaic_noaa <- extend(mosaic_noaa, mosaic_berk)
# #plot(is.na(mosaic_noaa) & is.na(mosaic_berk))
# 
# 
# ## turn into dfs
# df_berk <- data.frame(rasterToPoints(mosaic_berk))
# df_noaa <- data.frame(rasterToPoints(mosaic_noaa))
# colnames(df_noaa)[3] = "noaa"
# colnames(df_berk)[3] = "berk"
# plotdata <- full_join(df_berk, df_noaa)
# plotdata$is_missing_data = ifelse(!is.na(plotdata$berk) | !is.na(plotdata$noaa),
#                                   "No", 
#                                   "Yes")
# plotdata %>% 
#   ggplot(., aes(x = x, y = y, fill = is_missing_data)) +
#   geom_raster() +
#   coord_fixed() +
#   theme_void() 
# 
# ## turn into dfs
# df_berk <- data.frame(rasterToPoints(readRDS(paste("data-processed/spectral-change-files/BerkeleyEarth_s_stack_list_PSD_low.rds", sep = ""))[[1]]))
# df_noaa <- data.frame(rasterToPoints(readRDS(paste("data-processed/spectral-change-files/NOAA-OISST_s_stack_list_PSD_low.rds", sep = ""))[[1]]))
# colnames(df_noaa)[10] = "noaa"
# colnames(df_berk)[30] = "berk"
# df_noaa = select(df_noaa, x,y,noaa)
# df_berk = select(df_berk, x,y,berk)
# plotdata <- full_join(df_berk, df_noaa)
# plotdata$is_missing_data = ifelse(!is.na(plotdata$berk) | !is.na(plotdata$noaa),
#                                   "No", 
#                                   "Yes")
# plotdata %>% 
#   ggplot(., aes(x = x, y = y, fill = is_missing_data)) +
#   geom_raster() +
#   coord_fixed() +
#   theme_void() 










# ## now, fit model and plot model prediction instead of geom_smooth:
# library(nlme) 
# library(MuMIn)
# 
# ## for each gcm model, fit a model with location as random effect on intercept and slope 
# data <- c()
# mod_list <- list()
# all_predictions <- c()
# all <- c()
# for (i in 1:length(gcm_models)) {
#   
#   if(i == 7) {
#     mod_list[[i]] <- NA
#     
#     i = i+1
#   } 
#   else {
#     p = folders[i]
#     gcm = gcm_models[i]
#     
#     cur <- read.csv(paste("data-processed/spectral-change-files/", gcm, "_all-se-over-time_5_years.csv", sep = ""))
#     
#     cur <- cur %>%
#       select(x, y, s_spec_exp_PSD_low, window_start_year) %>% ## select only low seasonal spectral exponent slope 
#       distinct()
#     
#     cur$gcm <- gcm
#     cur$window_number <- gcm
#     
#     cur %>%
#       ggplot(., aes(x = window_start_year, y = s_spec_exp_PSD_low)) +
#       theme_light() +
#       labs(x = "Time window start year", y = "Spectral exponent") +
#       geom_smooth(method = "lm") +
#       geom_point() 
#     
#     ## make unique location id 
#     cur$loc_id <- paste(cur$lat, cur$lon, sep = "_")
#     
#     all <- rbind(all, cur)
#     
#     ## get rid of nas:
#     cur <- filter(cur, !is.na(s_spec_exp_PSD_low))
#     
#     mod <- lme(s_spec_exp_PSD_low ~ window_start_year, 
#                
#                random = ~window_start_year|loc_id,
#                
#                control = lmeControl(msMaxIter = 1000, msMaxEval = 1000),
#                
#                data = cur)
#     
#     mod_list[[i]] <- mod
#     
#     ## r.squaredGLMM(mod)
#     
#     ## make predictions
#     new_data <- data.frame(window_start_year = seq(from = 1871, to = 2100, by = 0.1))
#     
#     predictions = predict(mod, new_data, level = 0, se.fit = T, re.form = NA)
#     
#     fitted_predictions <- new_data %>%
#       mutate(pred_exponent = predictions$fit,
#              pred_exponent_SE = predictions$se.fit)
#     fitted_predictions$gcm = gcm_models[i]
#     
#     all_predictions = rbind(all_predictions, fitted_predictions)
#     
#   }
#   
# }
# names(mod_list) <- gcm_models
# #saveRDS(gcm_models, "data-processed/spectral-change-files/gcm-models.rds")
# 
# plot_data <- left_join(all, cur)
# #saveRDS(plot_data, "data-processed/spectral-change-files/plot-data.rds")
# 
# 
# ## plot the predictions altogether 
# all_predictions %>%
#   mutate(colour = ifelse(gcm == "NOAA-OISST",
#                          "Observed sea surface temperature",
#                          ifelse(gcm == "BerkeleyEarth", "Observed air surface temperature",
#                                 "Model-predicted air surface temperature"))) %>%
#   ggplot(., aes(x = window_start_year, y = pred_exponent, colour = colour)) +
#   theme_light() +
#   geom_ribbon(aes(ymin = pred_exponent - pred_exponent_SE, ymax = pred_exponent + pred_exponent_SE,
#                   group = gcm), 
#                 alpha = 0.3) +
#   labs(x = "Time window start year", y = "Spectral exponent", colour = "Dataset:") +
#  # geom_smooth(method = "lm")  +
#   geom_point() +
#   scale_color_manual(values = pal[c(1,3,5)]) +
#   theme(legend.position = "bottom")
# 
# #  geom_point(aes(x = window_start_year, y = mean_spec_exp, colour = colour))







