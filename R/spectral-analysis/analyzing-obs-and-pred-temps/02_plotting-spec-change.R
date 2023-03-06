## bringing all spectral exponent data together into one plot 
library(tidyverse)
library(PNWColors)

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

###################################################################
###     calling functions on spectral exponent files             ## 
###################################################################
num = 1
widths = c("5_years", "6_years", "7_years", "8_years", "9_years", "10_years")
while(num <= length(widths)) {
  ## read in and row bind csvs
  all <- c()
  for (i in 1:length(gcm_models)) {
    p = folders[i]
    gcm = gcm_models[i]
    
    cur <- read.csv(paste(p, gcm, "_average-se-over-time_", widths[num], ".csv", sep = ""))
    cur$gcm <- gcm
    
    all <- rbind(all, cur)
  }
  
  ## now plot them
  ## create colour palette:
  pal = pnw_palette("Sunset",5, type = "discrete")
  pnw_palette("Starfish",5, type = "discrete")
  
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
                             ifelse(gcm == "BerkeleyEarth", "Observed air surace temperature",
                                    "Model-predicted air surface temperature"))) %>%
      ggplot(., aes(x = year, y = mean_spec_exp, colour = colour, group = group)) +
      theme_light() +
      geom_errorbar(aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                    alpha = 0.3) +
      labs(x = "Time window start year", y = "Mean spectral exponent", colour = "Dataset:") +
      geom_smooth(method = "lm") +
      geom_point() +
      facet_wrap(~exponent_type, nrow = 3) +
      scale_color_manual(values = pal[c(1,3,5)]) +
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
                             ifelse(gcm == "BerkeleyEarth", "Observed air surace temperature",
                                    "Model-predicted air surface temperature"))) %>%
      ggplot(., aes(x = year, y = mean_spec_exp, colour = colour, group = group)) +
      theme_light() +
      geom_errorbar(aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                    alpha = 0.3) +
      labs(x = "Time window start year", y = "Mean spectral exponent", colour = "Dataset:") +
      geom_smooth(method = "lm")  +
      geom_point() +
      facet_wrap(~exponent_type) +
      scale_color_manual(values = pal[c(1,3,5)]) +
      theme(legend.position = "bottom")
    
    ggsave(plot, path = "figures/spectral-analysis_GCMs/", 
           filename = paste("all-datasets_", widths[num], ".png", sep = ""), device = "png",
           width = 9, height = 5)
  }
  
  num = num + 1
}


## now make maps :-)
num = 1
widths = c("5_years", "6_years", "7_years", "8_years", "9_years", "10_years")
while(num <= length(widths)) {
  ## read in and row bind csvs
  all <- c()
  for (i in 1:length(gcm_models)) {
    p = folders[i]
    gcm = gcm_models[i]
    
    cur <- read.csv(paste(p, gcm, "_average-se-over-time_", widths[num], ".csv", sep = ""))
    cur$gcm <- gcm
    
    all <- rbind(all, cur)
  }
  
  ## now plot them
  all %>%
    spread(gcm, spec_e)
  
  
  
  
  
  num = num + 1
}





