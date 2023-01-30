## calculating local climate velocities 

## first: install the package VoCC from github:
if (!"remotes" %in% rownames(installed.packages())) {
  install.packages("remotes")
  Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
  remotes::install_github("JorGarMol/VoCC")
}

library(VoCC)
library(tidyverse)
library(raster)
select <- dplyr::select


########################################################################
#####    1.  transform spectral exponent data into a rasterStack  ######
########################################################################
## function to convert output from script 04 into 6 rasterStacks (one for each sliding window width, from 5-10 years) that can be used with the gVoCC functions 
create_rasterStack <- function(p, gcm, gcm_num) {
  
  se_filenames <- readRDS(paste(p, "se_filenames.rds",  sep = ""))
  if(!gcm %in% c("BerkeleyEarth", "NOAA-OISST")) {
    se_filenames <- str_replace_all(se_filenames, "CMIP5-GCMs", "/Volumes/NIKKI/CMIP5-GCMs")
  }
  
  ## combine all spectral exponent csvs into one big dataframe
  file = 1
  while (file < length(se_filenames) + 1) {
    if(file.exists(se_filenames[file])) {
      if (file == 1) {
        spec_exp <- read.csv(se_filenames[file])
      }
      else {
        spec_exp <- rbind(spec_exp, read.csv(se_filenames[file]))
      }
    }
    print(paste("Reading file #", file, "/", length(se_filenames), sep = ""))
    file = file + 1
  }
  
  if(!gcm %in% c("BerkeleyEarth")) {
    ## change lon so it matches sea surface temperature 
    spec_exp <- spec_exp %>%
      mutate(lon = ifelse(lon >= 180, lon - 180, lon + 178))
  }
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
              res = 1)
  
  ## for observed temperature data, remove ones missing too many observations 
  if(gcm %in% c("BerkeleyEarth", "NOAA-OISST")) {
    
    ## get rid of ones with large chunks of ts missing
    mosaic <- readRDS(paste(p, "missing-data-count.rds", sep = ""))
    
    if(gcm == "NOAA-OISST") {
      ## change lon so it matches 
      mosaic = data.frame(rasterToPoints(mosaic))
      mosaic = mutate(mosaic, x = ifelse(x >= 180, x - 180, x + 178))
      mosaic <- rasterFromXYZ(mosaic)
      plot(mosaic)
    }
  
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
  }
  
  ## reorder time_window_width so list elements are in order of increasing time window widths:
  spec_exp$time_window_width <- factor(spec_exp$time_window_width, levels = 
                                         c("5 years", "6 years", "7 years", "8 years",
                                           "9 years", "10 years"))
  
  ## split spectral exponent data by window widths:
  ww_split <- split(spec_exp, spec_exp$time_window_width)
  
  ## for each window width:
  i = 1 
  while (i < length(ww_split) + 1) {
    
    ## split by window_start_year (each window step), loop through them and create raster layer for each
    ww <- ww_split[[i]]
    step_split <- split(ww, ww$window_start_year)
    step = length(step_split) ## loop backwards since rasterStacks add layers to beginning
    while (step >= 1) {
      
      ## pull out each different slope estimate
      l_sp_df_PSD_low <- step_split[[step]] %>%
        select(lon, lat, l_spec_exp_PSD_low)
      
      s_sp_df_PSD_low <- step_split[[step]] %>%
        select(lon, lat, s_spec_exp_PSD_low) 
      
      l_sp_df_AWC <- step_split[[step]] %>%
        select(lon, lat, l_spec_exp_AWC)
      
      s_sp_df_AWC <- step_split[[step]] %>%
        select(lon, lat, s_spec_exp_AWC) 
      
      l_sp_df_PSD_high <- step_split[[step]] %>%
        select(lon, lat, l_spec_exp_PSD_high)
      
      s_sp_df_PSD_high <- step_split[[step]] %>%
        select(lon, lat, s_spec_exp_PSD_high) 
      
      l_sp_df_PSD_all <- step_split[[step]] %>%
        select(lon, lat, l_spec_exp_PSD_all)
      
      s_sp_df_PSD_all <- step_split[[step]] %>%
        select(lon, lat, s_spec_exp_PSD_all)
      
      
      ## create raster layer:
      l_layer_PSD_low <- rasterFromXYZ(l_sp_df_PSD_low)
      s_layer_PSD_low <- rasterFromXYZ(s_sp_df_PSD_low)
      l_layer_AWC <- rasterFromXYZ(l_sp_df_AWC)
      s_layer_AWC <- rasterFromXYZ(s_sp_df_AWC)
      l_layer_PSD_high <- rasterFromXYZ(l_sp_df_PSD_high)
      s_layer_PSD_high <- rasterFromXYZ(s_sp_df_PSD_high)
      l_layer_PSD_all <- rasterFromXYZ(l_sp_df_PSD_all)
      s_layer_PSD_all <- rasterFromXYZ(s_sp_df_PSD_all)
      #plot(s_layer_PSD_high)
      
      ## add to temporary rasterstack:
      if (step == length(step_split)) {
        l_temp_stack_PSD_low <- l_layer_PSD_low
        s_temp_stack_PSD_low <- s_layer_PSD_low
        l_temp_stack_AWC <- l_layer_AWC
        s_temp_stack_AWC <- s_layer_AWC
        l_temp_stack_PSD_high <- l_layer_PSD_high
        s_temp_stack_PSD_high <- s_layer_PSD_high
        l_temp_stack_PSD_all <- l_layer_PSD_all
        s_temp_stack_PSD_all <- s_layer_PSD_all
      }
      else {
        ## make sure extents are the same
        ## (won't be for layers with time windows that had all data missing)
        l_layer_PSD_low = extend(l_layer_PSD_low, l_temp_stack_PSD_low)
        s_layer_PSD_low = extend(s_layer_PSD_low, s_temp_stack_PSD_low)
        l_layer_AWC = extend(l_layer_AWC, l_temp_stack_AWC)
        s_layer_AWC = extend(s_layer_AWC, s_temp_stack_AWC)
        l_layer_PSD_high = extend(l_layer_PSD_high, l_temp_stack_PSD_high)
        s_layer_PSD_high = extend(s_layer_PSD_high, s_temp_stack_PSD_high)
        l_layer_PSD_all = extend(l_layer_PSD_all, l_temp_stack_PSD_all)
        s_layer_PSD_all = extend(s_layer_PSD_all, s_temp_stack_PSD_all)
        
        l_temp_stack_PSD_low <- addLayer(l_layer_PSD_low, l_temp_stack_PSD_low)
        s_temp_stack_PSD_low <- addLayer(s_layer_PSD_low, s_temp_stack_PSD_low)
        l_temp_stack_AWC <- addLayer(l_layer_AWC, l_temp_stack_AWC)
        s_temp_stack_AWC <- addLayer(s_layer_AWC, s_temp_stack_AWC)
        l_temp_stack_PSD_high <- addLayer(l_layer_PSD_high, l_temp_stack_PSD_high)
        s_temp_stack_PSD_high <- addLayer(s_layer_PSD_high, s_temp_stack_PSD_high)
        l_temp_stack_PSD_all <- addLayer(l_layer_PSD_all, l_temp_stack_PSD_all)
        s_temp_stack_PSD_all <- addLayer(s_layer_PSD_all, s_temp_stack_PSD_all)
      }
      
      ## move to nested for loop
      step = step - 1
    }
    
    names(l_temp_stack_PSD_low) <- names(l_temp_stack_PSD_high) <- names(s_temp_stack_PSD_low) <- 
      names(s_temp_stack_PSD_high) <- names(l_temp_stack_AWC) <- names(s_temp_stack_AWC) <- 
      paste("window", 1:nlayers(l_temp_stack_PSD_low), sep = "_")
    
    ## save temporary raster stack: 
    if (i == 1) {
      l_stack_list_PSD_low <- list(l_temp_stack_PSD_low)
      s_stack_list_PSD_low <- list(s_temp_stack_PSD_low)
      l_stack_list_AWC <- list(l_temp_stack_AWC)
      s_stack_list_AWC <- list(s_temp_stack_AWC)
      l_stack_list_PSD_high <- list(l_temp_stack_PSD_high)
      s_stack_list_PSD_high <- list(s_temp_stack_PSD_high)
      l_stack_list_PSD_all <- list(l_temp_stack_PSD_all)
      s_stack_list_PSD_all <- list(s_temp_stack_PSD_all)
    }
    else {
      l_stack_list_PSD_low <- append(l_stack_list_PSD_low, l_temp_stack_PSD_low)
      s_stack_list_PSD_low <- append(s_stack_list_PSD_low, s_temp_stack_PSD_low)
      l_stack_list_AWC <- append(l_stack_list_AWC, l_temp_stack_AWC)
      s_stack_list_AWC <- append(s_stack_list_AWC, s_temp_stack_AWC)
      l_stack_list_PSD_high <- append(l_stack_list_PSD_high, l_temp_stack_PSD_high)
      s_stack_list_PSD_high <- append(s_stack_list_PSD_high, s_temp_stack_PSD_high)
      l_stack_list_PSD_all <- append(l_stack_list_PSD_all, l_temp_stack_PSD_all)
      s_stack_list_PSD_all <- append(s_stack_list_PSD_all, s_temp_stack_PSD_all)
    }
    
    ## move to next time window width
    i = i + 1
  }
  
  ## name the list items 
  names(l_stack_list_PSD_low) <- names(l_stack_list_PSD_high) <- names(s_stack_list_PSD_low) <-
    names(s_stack_list_PSD_high) <- names(l_stack_list_AWC) <- names(s_stack_list_AWC) <-
    names(l_stack_list_PSD_all) <- names(s_stack_list_PSD_all) <-
    names(ww_split)
  
  ## save the rasterstack 
  saveRDS(l_stack_list_PSD_low, paste(p, gcm, "_l_stack_list_PSD_low.rds", sep = ""))
  saveRDS(s_stack_list_PSD_low, paste(p, gcm, "_s_stack_list_PSD_low.rds", sep = ""))
  saveRDS(l_stack_list_AWC, paste(p, gcm, "_l_stack_list_AWC.rds", sep = ""))
  saveRDS(s_stack_list_AWC, paste(p, gcm, "_s_stack_list_AWC.rds", sep = ""))
  saveRDS(l_stack_list_PSD_high, paste(p, gcm, "_l_stack_list_PSD_high.rds", sep = ""))
  saveRDS(s_stack_list_PSD_high, paste(p, gcm, "_s_stack_list_PSD_high.rds", sep = ""))
  saveRDS(l_stack_list_PSD_all, paste(p, gcm, "_l_stack_list_PSD_all.rds", sep = ""))
  saveRDS(s_stack_list_PSD_all, paste(p, gcm, "_s_stack_list_PSD_all.rds", sep = ""))
  
  stacks <- list(l_stack_list_PSD_low, s_stack_list_PSD_low, 
                 l_stack_list_AWC, s_stack_list_AWC,
                 l_stack_list_PSD_high, s_stack_list_PSD_high,
                 l_stack_list_PSD_all, s_stack_list_PSD_all)
  
  ## return the 6 lists of rasterStacks
  return(stacks)
}

#########################################################################################
##### 2. calculate average change in spectral exponent for each time series in gcm ######
#########################################################################################
df_average_se_over_time <- function(p, gcm, gcm_num) {
  
  num = 1
  widths = c("5 years", "6 years", "7 years", "8 years", "9 years", "10 years")
  while(num <= length(widths)) {
    ## read in raster stack layers for the time window 
    s_stack_tas_PSD_low <- readRDS(paste(p, gcm, "_s_stack_list_PSD_low.rds", sep = ""))[[num]] 
    s_stack_tas_PSD_high <- readRDS(paste(p, gcm, "_s_stack_list_PSD_high.rds", sep = ""))[[num]] 
    
    if(num == 6) {
      s_stack_tas_AWC <- readRDS(paste(p, gcm, "_s_stack_list_AWC.rds", sep = ""))[[6]] 
      l_stack_tas_PSD_low <- readRDS(paste(p, gcm, "_l_stack_list_PSD_low.rds", sep = ""))[[6]] 
      l_stack_tas_AWC <- readRDS(paste(p, gcm, "_l_stack_list_AWC.rds", sep = ""))[[6]] 
      l_stack_tas_PSD_high <- readRDS(paste(p, gcm, "_l_stack_list_PSD_high.rds", sep = ""))[[6]]
      
      ## make list of sensitivity layers
      list_tas <- c(s_stack_tas_PSD_low, s_stack_tas_PSD_high, s_stack_tas_AWC,
                    l_stack_tas_PSD_low, l_stack_tas_PSD_high, l_stack_tas_AWC)

    }
    else {
      list_tas <- c(s_stack_tas_PSD_low, s_stack_tas_PSD_high)
    }
    
    if(!gcm %in% c("NOAA-OISST")) {
      ## crop to only land
      land <- raster("data-processed/masks/cmip5-land.grd") 
      extent(land) <- c(0, 360, -90, 90)
      land <- crop(land, list_tas[[1]])
      land[land != 1] <- NA
      
      list_tas <- sapply(list_tas, FUN = mask, land)
    }
    
    list_tas <- sapply(list_tas, FUN = stack)
    
    if(gcm == "BerkeleyEarth") {
      years = seq(1880, 2021-num-4, by = num+4)
    }
    else if(gcm == "NOAA-OISST") {
      years = seq(1981, 2022-num-4, by = num+4)
    }
    else {
      years = seq(1871, 2100-num-4, by = num+4)
    }
    
    ## convert raster layer to data frame
    df_low <- data.frame(rasterToPoints(list_tas[[1]]))
    df_low <- gather(df_low, key = "window_number", value = "spec_exp", c(3:ncol(df_low)))
    df_low$exponent_type = "s_PSD_low"
    df_low$time_window_width = widths[num]
   
    df_high <- data.frame(rasterToPoints(list_tas[[2]]))
    df_high <- gather(df_high, key = "window_number", value = "spec_exp", c(3:ncol(df_high)))
    df_high$exponent_type = "s_PSD_high"
    df_high$time_window_width = widths[num]
    
    average_low <- df_low %>%
      group_by(window_number) %>% ## group data by window start year
      mutate(mean_spec_exp = mean(spec_exp),
             sd_spec_exp = sd(spec_exp))  %>% ## calculate average spectral exponent across all locations for each window 
      ungroup() %>%
      select(-x, -y, -spec_exp) %>%
      unique() %>%
      select(-window_number)
    average_low$year = years
    
    average_high <- df_high %>%
      group_by(window_number) %>% ## group data by window start year
      mutate(mean_spec_exp = mean(spec_exp),
             sd_spec_exp = sd(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
      ungroup() %>%
      select(-x, -y, -spec_exp) %>%
      unique() %>%
      select(-window_number)
    average_high$year = years
    
    all = rbind(average_high, average_low)
    all_data = rbind(df_high, df_low)
    
    if(num == 6) {
      df_awc_s <- data.frame(rasterToPoints(list_tas[[3]]))
      df_awc_s <- gather(df_awc_s, key = "window_number", value = "spec_exp", c(3:ncol(df_awc_s)))
      df_awc_s$exponent_type = "s_AWC"
      df_awc_s$time_window_width = widths[num]
      
      average_awc_s <- df_awc_s %>%
        group_by(window_number) %>% ## group data by window start year
        mutate(mean_spec_exp = mean(spec_exp),
               sd_spec_exp = sd(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
        ungroup() %>%
        select(-x, -y, -spec_exp) %>%
        unique() %>%
        select(-window_number)
      average_awc_s$year = years
      
      df_low_l <- data.frame(rasterToPoints(list_tas[[4]]))
      df_low_l <- gather(df_low_l, key = "window_number", value = "spec_exp", c(3:ncol(df_low_l)))
      df_low_l$exponent_type = "l_PSD_low"
      df_low_l$time_window_width = widths[num]
      
      df_high_l <- data.frame(rasterToPoints(list_tas[[5]]))
      df_high_l <- gather(df_high_l, key = "window_number", value = "spec_exp", c(3:ncol(df_high_l)))
      df_high_l$exponent_type = "l_PSD_high"
      df_high_l$time_window_width = widths[num]
      
      df_awc_l <- data.frame(rasterToPoints(list_tas[[6]]))
      df_awc_l <- gather(df_awc_l, key = "window_number", value = "spec_exp", c(3:ncol(df_awc_l)))
      df_awc_l$exponent_type = "l_AWC"
      df_awc_l$time_window_width = widths[num]
      
      average_low_l <- df_low_l %>%
        group_by(window_number) %>% ## group data by window start year
        mutate(mean_spec_exp = mean(spec_exp),
               sd_spec_exp = sd(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
        ungroup() %>%
        select(-x, -y, -spec_exp) %>%
        unique() %>%
        select(-window_number)
      average_low_l$year = years
      
      average_high_l <- df_high_l %>%
        group_by(window_number) %>% ## group data by window start year
        mutate(mean_spec_exp = mean(spec_exp),
               sd_spec_exp = sd(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
        ungroup() %>%
        select(-x, -y, -spec_exp) %>%
        unique() %>%
        select(-window_number)
      average_high_l$year = years
      
      average_awc_l <- df_awc_l %>%
        group_by(window_number) %>% ## group data by window start year
        mutate(mean_spec_exp = mean(spec_exp),
               sd_spec_exp = sd(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
        ungroup() %>%
        select(-x, -y, -spec_exp) %>%
        unique() %>%
        select(-window_number)
      average_awc_l$year = years
      
      all = rbind(all, average_awc_s) %>%
        rbind(., average_low_l) %>%
        rbind(., average_high_l) %>%
        rbind(., average_awc_l)
      
      all_data = rbind(all_data, df_awc_s) %>%
        rbind(., df_low_l) %>%
        rbind(., df_high_l) %>%
        rbind(., df_awc_l)
      
      key = data.frame(window_number = unique(df_low$window_number),
                       year = years)
      all_data = left_join(all_data, key)
      
    }
    
    plot = all %>%
      ggplot(., aes(x = year, y = mean_spec_exp, colour = exponent_type)) + geom_point() +
      theme_light() +
      geom_errorbar(aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                    alpha = 0.3) +
      labs(x = "Time window start year", y = "Mean spectral exponent", colour = "Sensitivity set:") +
      geom_smooth(method = "lm")
    
    width = str_replace_all(widths[num], " ", "_")
    
    ggsave(plot, 
           filename = paste("figures/spectral-analysis_GCMs/", gcm, "_global-se-over-time_", width, ".png", sep = ""), 
           width = 6, height = 4, device = "png")
    
    ## write out
    write.csv(all, paste(p, gcm, "_average-se-over-time_", width, ".csv", sep = ""))
    write.csv(all_data, paste(p, gcm, "_all-se-over-time_", width, ".csv", sep = ""))
    
    num = num + 1
  }
  
  return(NA)
}

raster_avg_change <- function(p, gcm) {
  ## calculate average change in spectral exponent for each time series, across all sliding window widths 
  se_filenames <- readRDS(paste(p, "se_filenames.rds",  sep = ""))
  if(!gcm %in% c("BerkeleyEarth", "NOAA-OISST")) {
    se_filenames <- str_replace_all(se_filenames, "CMIP5-GCMs", "/Volumes/NIKKI/CMIP5-GCMs")
  }
  
  ## combine all spectral exponent csvs into one big dataframe
  file = 1
  while (file < length(se_filenames) + 1) {
    if(file.exists(se_filenames[file])) {
      if (file == 1) {
        spec_exp <- read.csv(se_filenames[file])
      }
      else {
        spec_exp <- rbind(spec_exp, read.csv(se_filenames[file]))
      }
      print(paste("Reading file #", file, "/", length(se_filenames), sep = ""))
    }
    file = file + 1
  }
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
              res = 1)
  
  ## reorder time_window_width so list elements are in order of increasing time window widths:
  spec_exp$time_window_width <- factor(spec_exp$time_window_width, levels = 
                                         c("5 years", "6 years", "7 years", "8 years",
                                           "9 years", "10 years"))
  
  ## loop through and save each window width
  widths = unique(spec_exp$time_window_width)
  num = 1 
  while(num <= length(widths)) {
    
    ## get rid of missing time series BE
    if(gcm %in% c("BerkeleyEarth", "NOAA-OISST")) {
      ## get rid of ones with large chunks of ts missing
      mosaic <- readRDS(paste(p, "missing-data-count.rds", sep = ""))
      
      if(gcm == "NOAA-OISST") {
        ## change lon so it matches 
        mosaic = data.frame(rasterToPoints(mosaic))
        mosaic = mutate(mosaic, x = ifelse(x >= 180, x - 180, x + 178))
        mosaic <- rasterFromXYZ(mosaic)
        plot(mosaic)
      }
      
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
    }
    
    tas <- spec_exp %>%
      filter(time_window_width == widths[num]) %>%
      select(lon, lat, l_estimate_PSD_low, s_estimate_PSD_low, l_estimate_PSD_high, s_estimate_PSD_high,
             l_estimate_AWC, s_estimate_AWC) %>%
      unique() 
    
    if(gcm != "BerkeleyEarth") {
      tas <- mutate(tas, lon = ifelse(lon >= 180, lon - 180, lon + 178))
    }
      
    raster_tas <- rasterFromXYZ(tas)
    raster_tas <- extend(raster_tas, c(0, 360, -90, 90))
    
    ## save:
    width = str_replace_all(widths[num], " ", "_")
    saveRDS(raster_tas, paste(p, gcm, "_spectral-change_multifrac_", width, ".rds", sep = ""))
    
    num = num + 1
  }
  
  return(NA)
}




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
p = folders[gcm_num]
gcm = gcm_models[gcm_num]

stck = create_rasterStack(p = p, gcm = gcm, gcm_num = gcm_num)
df_average_se_over_time(p = p, gcm = gcm, gcm_num = gcm_num)
raster_avg_change(p = p, gcm = gcm)



write.csv(spec_exp, paste(p, gcm, "_spec_exp.csv",  sep = ""), row.names = F)
          
          
          
