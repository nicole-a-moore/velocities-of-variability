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
create_rasterStack <- function(p, gcm) {
  
  se_filenames <- readRDS(paste(p, "se_filenames.rds",  sep = ""))
  se_filenames <- str_replace_all(se_filenames, "CMIP5-GCMs", "/Volumes/NIKKI/CMIP5-GCMs")
  
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
  
  ## change lon so it matches sea surface temperature 
  spec_exp <- spec_exp %>%
    mutate(lon = ifelse(lon >= 180, lon - 180, lon + 178))
  
  r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, 
              crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"),
              res = 1)
  
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
      #plot(l_layer_PSD_high)
      
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
df_avergage_se_over_time <- function(p, gcm, i) {
  
  num = 1
  widths = c("5 years", "6 years", "7 years", "8 years", "9 years", "10 years")
  while(num <= length(widths)) {
    ## read in raster stack layers for the time window 
    s_stack_tas_PSD_low <- readRDS(paste(p, gcm_models[i], "_s_stack_list_PSD_low.rds", sep = ""))[[num]] 
    s_stack_tas_PSD_high <- readRDS(paste(p, gcm_models[i], "_s_stack_list_PSD_high.rds", sep = ""))[[num]] 
    
    if(num == 6) {
      s_stack_tas_AWC <- readRDS(paste(p, gcm_models[i], "_s_stack_list_AWC.rds", sep = ""))[[6]] 
      l_stack_tas_PSD_low <- readRDS(paste(p, gcm_models[i], "_l_stack_list_PSD_low.rds", sep = ""))[[6]] 
      l_stack_tas_AWC <- readRDS(paste(p, gcm_models[i], "_l_stack_list_AWC.rds", sep = ""))[[6]] 
      l_stack_tas_PSD_high <- readRDS(paste(p, gcm_models[i], "_l_stack_list_PSD_high.rds", sep = ""))[[6]]
      
      ## make list of sensitivity layers
      list_tas <- c(s_stack_tas_PSD_low, s_stack_tas_PSD_high, s_stack_tas_AWC,
                    l_stack_tas_PSD_low, l_stack_tas_PSD_high, l_stack_tas_AWC)

    }
    else {
      list_tas <- c(s_stack_tas_PSD_low, s_stack_tas_PSD_high)
    }
    ## crop to only land
    land <- raster("data-processed/masks/cmip5-land.grd") 
    extent(land) <- c(0, 360, -90, 90)
    land <- crop(land, list_tas[[1]])
    land[land != 1] <- NA
    
    list_tas <- sapply(list_tas, FUN = mask, land)
    list_tas <- sapply(list_tas, FUN = stack)
    
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
      select(-x, -y, -window_number, -spec_exp) %>%
      unique() 
    average_low$year = seq(1871, 2100-num-4, by = num+4)
    
    average_high <- df_high %>%
      group_by(window_number) %>% ## group data by window start year
      mutate(mean_spec_exp = mean(spec_exp),
             sd_spec_exp = sd(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
      ungroup() %>%
      select(-x, -y, -window_number, -spec_exp) %>%
      unique() 
    average_high$year = seq(1871, 2100-num-4, by = num+4)
    
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
        select(-x, -y, -window_number, -spec_exp) %>%
        unique() 
      average_awc_s$year = seq(1871, 2100-num-4, by = num+4)
      
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
        select(-x, -y, -window_number, -spec_exp) %>%
        unique() 
      average_low_l$year = seq(1871, 2100-num-4, by = num+4)
      
      average_high_l <- df_high_l %>%
        group_by(window_number) %>% ## group data by window start year
        mutate(mean_spec_exp = mean(spec_exp),
               sd_spec_exp = sd(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
        ungroup() %>%
        select(-x, -y, -window_number, -spec_exp) %>%
        unique() 
      average_high_l$year = seq(1871, 2100-num-4, by = num+4)
      
      average_awc_l <- df_awc_l %>%
        group_by(window_number) %>% ## group data by window start year
        mutate(mean_spec_exp = mean(spec_exp),
               sd_spec_exp = sd(spec_exp)) %>% ## calculate average spectral exponent across all locations for each window 
        ungroup() %>%
        select(-x, -y, -window_number, -spec_exp) %>%
        unique() 
      average_awc_l$year = seq(1871, 2100-num-4, by = num+4)
      
      all = rbind(all, average_awc_s) %>%
        rbind(., average_low_l) %>%
        rbind(., average_high_l) %>%
        rbind(., average_awc_l)
      
      all_data = rbind(all_data, df_awc_s) %>%
        rbind(., df_low_l) %>%
        rbind(., df_high_l) %>%
        rbind(., df_awc_l)
      
      key = data.frame(window_number = unique(df_low$window_number),
                       year = seq(1871, 2100-num-4, by = num+4))
      all_data = left_join(all_data, key)
      
    }
    
    plot = all %>%
      ggplot(., aes(x = year, y = mean_spec_exp, colour = exponent_type)) + geom_point() +
      theme_light() +
      geom_errorbar(aes(ymin = mean_spec_exp - sd_spec_exp, ymax = mean_spec_exp + sd_spec_exp), 
                    alpha = 0.3) +
      labs(x = "Time window start year", y = "Mean spectral exponent", colour = "Sensitivity set:") +
      geom_smooth(method = "lm")
    
    ggsave(plot, 
           filename = paste("figures/spectral-analysis_GCMs/", gcm_models[i], "_global-se-over-time_", width, ".png", sep = ""), 
           width = 6, height = 4, device = "png")
    
    ## write out
    width = str_replace_all(widths[num], " ", "_")
    write.csv(all, paste(p, gcm_models[i], "_average-se-over-time_", width, ".csv", sep = ""))
    write.csv(all_data, paste(p, gcm_models[i], "_all-se-over-time_", width, ".csv", sep = ""))
    
    num = num + 1
  }
  
  return(NA)
}

raster_avg_change <- function(p, gcm) {
  ## calculate average change in spectral exponent for each time series, across all sliding window widths 
  se_filenames <- readRDS(paste(p, "se_filenames.rds",  sep = ""))
  se_filenames <- str_replace_all(se_filenames, "CMIP5-GCMs", "/Volumes/NIKKI/CMIP5-GCMs")
  
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
    tas <- spec_exp %>%
      filter(time_window_width == widths[num]) %>%
      select(lon, lat, l_estimate_PSD_low, s_estimate_PSD_low, l_estimate_PSD_high, s_estimate_PSD_high,
             l_estimate_AWC, s_estimate_AWC) %>%
      unique() %>%
      mutate(lon = ifelse(lon >= 180, lon - 180, lon + 178))
    
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
                "11_MRI-CGCM3_tas")

folders <- paste(path, gcm_models, "/", sep = "")

#################################################
###       make list of GCM file names          ## 
#################################################
CMCC_CMS_hist <- c("tas_day_CMCC-CMS_historical_r1i1p1_18700101-18791231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_18800101-18891231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_18900101-18991231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19000101-19091231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19100101-19191231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19200101-19291231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19300101-19391231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19400101-19491231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19500101-19591231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19600101-19691231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19700101-19791231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19800101-19891231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_19900101-19991231.nc",
                   "tas_day_CMCC-CMS_historical_r1i1p1_20000101-20051231.nc")
CMCC_CMS_rcp85 <- c("tas_day_CMCC-CMS_rcp85_r1i1p1_20060101-20091231.nc",
                    "tas_day_CMCC-CMS_rcp85_r1i1p1_20100101-20191231.nc",
                    "tas_day_CMCC-CMS_rcp85_r1i1p1_20200101-20291231.nc",
                    "tas_day_CMCC-CMS_rcp85_r1i1p1_20300101-20391231.nc",
                    "tas_day_CMCC-CMS_rcp85_r1i1p1_20400101-20491231.nc",
                    "tas_day_CMCC-CMS_rcp85_r1i1p1_20500101-20591231.nc",
                    "tas_day_CMCC-CMS_rcp85_r1i1p1_20600101-20691231.nc",
                    "tas_day_CMCC-CMS_rcp85_r1i1p1_20700101-20791231.nc",
                    "tas_day_CMCC-CMS_rcp85_r1i1p1_20800101-20891231.nc",
                    "tas_day_CMCC-CMS_rcp85_r1i1p1_20900101-21001231.nc")
CMCC_CMS_01 <- list(CMCC_CMS_hist, CMCC_CMS_rcp85)

GFDL_CM3_hist <- c("tas_day_GFDL-CM3_historical_r1i1p1_18700101-18741231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_18750101-18791231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_18800101-18841231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_18850101-18891231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_18900101-18941231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_18950101-18991231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19000101-19041231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19050101-19091231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19100101-19141231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19150101-19191231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19200101-19241231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19250101-19291231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19300101-19341231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19350101-19391231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19400101-19441231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19450101-19491231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19500101-19541231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19550101-19591231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19600101-19641231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19650101-19691231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19700101-19741231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19750101-19791231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19800101-19841231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19850101-19891231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19900101-19941231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_19950101-19991231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_20000101-20041231.nc",
                   "tas_day_GFDL-CM3_historical_r1i1p1_20050101-20051231.nc")
GFDL_CM3_rcp85 <- c("tas_day_GFDL-CM3_rcp85_r1i1p1_20060101-20101231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20110101-20151231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20160101-20201231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20210101-20251231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20260101-20301231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20310101-20351231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20360101-20401231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20410101-20451231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20460101-20501231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20510101-20551231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20560101-20601231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20610101-20651231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20660101-20701231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20710101-20751231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20760101-20801231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20810101-20851231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20860101-20901231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20910101-20951231.nc",
                    "tas_day_GFDL-CM3_rcp85_r1i1p1_20960101-21001231.nc")
GFDL_CM3_02 <- list(GFDL_CM3_hist, GFDL_CM3_rcp85)

GFDL_ESM2G_hist <- c("tas_day_GFDL-ESM2G_historical_r1i1p1_18660101-18701231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_18710101-18751231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_18760101-18801231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_18810101-18851231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_18860101-18901231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_18910101-18951231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_18960101-19001231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19010101-19051231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19060101-19101231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19110101-19151231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19160101-19201231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19210101-19251231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19260101-19301231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19310101-19351231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19360101-19401231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19410101-19451231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19460101-19501231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19510101-19551231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19560101-19601231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19610101-19651231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19660101-19701231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19710101-19751231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19760101-19801231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19810101-19851231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19860101-19901231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19910101-19951231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_19960101-20001231.nc",
                     "tas_day_GFDL-ESM2G_historical_r1i1p1_20010101-20051231.nc")
GFDL_ESM2G_rcp85 <- c("tas_day_GFDL-ESM2G_rcp85_r1i1p1_20060101-20101231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20110101-20151231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20160101-20201231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20210101-20251231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20260101-20301231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20310101-20351231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20360101-20401231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20410101-20451231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20460101-20501231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20510101-20551231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20560101-20601231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20610101-20651231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20660101-20701231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20710101-20751231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20760101-20801231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20810101-20851231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20860101-20901231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20910101-20951231.nc",
                      "tas_day_GFDL-ESM2G_rcp85_r1i1p1_20960101-21001231.nc")
GFDL_ESM2G_03 <- list(GFDL_ESM2G_hist, GFDL_ESM2G_rcp85)

HadGEM2_ES_hist <- c("tas_day_HadGEM2-ES_historical_r1i1p1_18691201-18791130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_18791201-18891130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_18891201-18991130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_18991201-19091130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19091201-19191130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19191201-19291130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19291201-19391130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19391201-19491130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19491201-19591130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19591201-19691130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19691201-19791130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19791201-19891130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19891201-19991130.nc",
                     "tas_day_HadGEM2-ES_historical_r1i1p1_19991201-20051130.nc")
HadGEM2_ES_rcp85 <- c("tas_day_HadGEM2-ES_rcp85_r1i1p1_20051201-20101130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20101201-20151130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20151201-20201130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20201201-20251130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20251201-20301130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20301201-20351130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20351201-20401130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20401201-20451130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20451201-20501130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20501201-20551130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20551201-20601130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20601201-20651130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20651201-20701130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20701201-20751130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20751201-20851130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20851201-20951130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20951201-20991130.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20991201-20991230.nc",
                      "tas_day_HadGEM2-ES_rcp85_r1i1p1_20991201-21091130.nc")
HadGEM2_ES_04 <- list(HadGEM2_ES_hist, HadGEM2_ES_rcp85)

inmcm4_hist <- c("tas_day_inmcm4_historical_r1i1p1_18700101-18791231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_18800101-18891231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_18900101-18991231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19000101-19091231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19100101-19191231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19200101-19291231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19300101-19391231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19400101-19491231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19500101-19591231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19600101-19691231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19700101-19791231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19800101-19891231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_19900101-19991231.nc",
                 "tas_day_inmcm4_historical_r1i1p1_20000101-20051231.nc")
inmcm4_rcp85 <- c("tas_day_inmcm4_rcp85_r1i1p1_20060101-20151231.nc",
                  "tas_day_inmcm4_rcp85_r1i1p1_20160101-20251231.nc",
                  "tas_day_inmcm4_rcp85_r1i1p1_20260101-20351231.nc",
                  "tas_day_inmcm4_rcp85_r1i1p1_20360101-20451231.nc",
                  "tas_day_inmcm4_rcp85_r1i1p1_20460101-20551231.nc",
                  "tas_day_inmcm4_rcp85_r1i1p1_20560101-20651231.nc",
                  "tas_day_inmcm4_rcp85_r1i1p1_20660101-20751231.nc",
                  "tas_day_inmcm4_rcp85_r1i1p1_20760101-20851231.nc",
                  "tas_day_inmcm4_rcp85_r1i1p1_20860101-20951231.nc",
                  "tas_day_inmcm4_rcp85_r1i1p1_20960101-21001231.nc")
inmcm4_05 <- list(inmcm4_hist, inmcm4_rcp85)

IPSL_CM5A_MR_hist <- c("tas_day_IPSL-CM5A-MR_historical_r1i1p1_18500101-18991231.nc",
                       "tas_day_IPSL-CM5A-MR_historical_r1i1p1_19000101-19491231.nc",
                       "tas_day_IPSL-CM5A-MR_historical_r1i1p1_19500101-19991231.nc",
                       "tas_day_IPSL-CM5A-MR_historical_r1i1p1_20000101-20051231.nc")
IPSL_CM5A_MR_rcp85 <- c("tas_day_IPSL-CM5A-MR_rcp85_r1i1p1_20060101-20551231.nc",
                        "tas_day_IPSL-CM5A-MR_rcp85_r1i1p1_20560101-21001231.nc")
IPSL_CM5A_MR_06 <- list(IPSL_CM5A_MR_hist, IPSL_CM5A_MR_rcp85)

MIROC_ESM_CHEM_hist <- c('tas_day_MIROC-ESM-CHEM_historical_r1i1p1_18500101-20051231.nc') 
MIROC_ESM_CHEM_rcp85 <- c("tas_day_MIROC-ESM-CHEM_rcp85_r1i1p1_20060101-21001231.nc")
MIROC_ESM_CHEM_07 <- list(MIROC_ESM_CHEM_hist, MIROC_ESM_CHEM_rcp85)

MIROC5_hist <- c("tas_day_MIROC5_historical_r1i1p1_18700101-18791231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_18800101-18891231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_18900101-18991231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19000101-19091231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19100101-19191231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19200101-19291231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19300101-19391231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19400101-19491231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19500101-19591231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19600101-19691231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19700101-19791231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19800101-19891231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_19900101-19991231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_20000101-20091231.nc")
MIROC5_rcp85 <- c("tas_day_MIROC5_rcp85_r1i1p1_20060101-20091231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_20100101-20191231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_20200101-20291231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_20300101-20391231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_20400101-20491231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_20500101-20591231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_20600101-20691231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_20700101-20791231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_20800101-20891231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_20900101-20991231.nc",
                  "tas_day_MIROC5_rcp85_r1i1p1_21000101-21001231.nc")
MIROC5_08 <- list(MIROC5_hist, MIROC5_rcp85)

MPI_ESM_LR_hist <- c("tas_day_MPI-ESM-LR_historical_r1i1p1_18700101-18791231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_18800101-18891231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_18900101-18991231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19000101-19091231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19100101-19191231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19200101-19291231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19300101-19391231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19400101-19491231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19500101-19591231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19600101-19691231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19700101-19791231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19800101-19891231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_19900101-19991231.nc",
                     "tas_day_MPI-ESM-LR_historical_r1i1p1_20000101-20051231.nc")
MPI_ESM_LR_rcp85 <- c("tas_day_MPI-ESM-LR_rcp85_r1i1p1_20060101-20091231.nc",
                      "tas_day_MPI-ESM-LR_rcp85_r1i1p1_20100101-20191231.nc",
                      "tas_day_MPI-ESM-LR_rcp85_r1i1p1_20200101-20291231.nc",
                      "tas_day_MPI-ESM-LR_rcp85_r1i1p1_20300101-20391231.nc",
                      "tas_day_MPI-ESM-LR_rcp85_r1i1p1_20400101-20491231.nc",
                      "tas_day_MPI-ESM-LR_rcp85_r1i1p1_20500101-20591231.nc",
                      "tas_day_MPI-ESM-LR_rcp85_r1i1p1_20600101-20691231.nc",
                      "tas_day_MPI-ESM-LR_rcp85_r1i1p1_20700101-20791231.nc",
                      "tas_day_MPI-ESM-LR_rcp85_r1i1p1_20800101-20891231.nc",
                      "tas_day_MPI-ESM-LR_rcp85_r1i1p1_20900101-21001231.nc")
MPI_ESM_LR_09 <- list(MPI_ESM_LR_hist, MPI_ESM_LR_rcp85)

MPI_ESM_MR_hist <- c("tas_day_MPI-ESM-MR_historical_r1i1p1_18700101-18791231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_18800101-18891231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_18900101-18991231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19000101-19091231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19100101-19191231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19200101-19291231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19300101-19391231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19400101-19491231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19500101-19591231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19600101-19691231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19700101-19791231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19800101-19891231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_19900101-19991231.nc",
                     "tas_day_MPI-ESM-MR_historical_r1i1p1_20000101-20051231.nc")
MPI_ESM_MR_rcp85 <- c("tas_day_MPI-ESM-MR_rcp85_r1i1p1_20060101-20091231.nc",
                      "tas_day_MPI-ESM-MR_rcp85_r1i1p1_20100101-20191231.nc",
                      "tas_day_MPI-ESM-MR_rcp85_r1i1p1_20200101-20291231.nc",
                      "tas_day_MPI-ESM-MR_rcp85_r1i1p1_20300101-20391231.nc",
                      "tas_day_MPI-ESM-MR_rcp85_r1i1p1_20400101-20491231.nc",
                      "tas_day_MPI-ESM-MR_rcp85_r1i1p1_20500101-20591231.nc",
                      "tas_day_MPI-ESM-MR_rcp85_r1i1p1_20600101-20691231.nc",
                      "tas_day_MPI-ESM-MR_rcp85_r1i1p1_20700101-20791231.nc",
                      "tas_day_MPI-ESM-MR_rcp85_r1i1p1_20800101-20891231.nc",
                      "tas_day_MPI-ESM-MR_rcp85_r1i1p1_20900101-21001231.nc")
MPI_ESM_MR_10 <- list(MPI_ESM_MR_hist, MPI_ESM_MR_rcp85)

MRI_CGCM3_hist <- c("tas_day_MRI-CGCM3_historical_r1i1p1_18700101-18791231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_18800101-18891231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_18900101-18991231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19000101-19091231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19100101-19191231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19200101-19291231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19300101-19391231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19400101-19491231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19500101-19591231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19600101-19691231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19700101-19791231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19800101-19891231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_19900101-19991231.nc",
                    "tas_day_MRI-CGCM3_historical_r1i1p1_20000101-20051231.nc")
MRI_CGCM3_rcp85 <- c("tas_day_MRI-CGCM3_rcp85_r1i1p1_20060101-20151231.nc",
                     "tas_day_MRI-CGCM3_rcp85_r1i1p1_20160101-20251231.nc",
                     "tas_day_MRI-CGCM3_rcp85_r1i1p1_20260101-20351231.nc",
                     "tas_day_MRI-CGCM3_rcp85_r1i1p1_20360101-20451231.nc",
                     "tas_day_MRI-CGCM3_rcp85_r1i1p1_20460101-20551231.nc",
                     "tas_day_MRI-CGCM3_rcp85_r1i1p1_20560101-20651231.nc",
                     "tas_day_MRI-CGCM3_rcp85_r1i1p1_20660101-20751231.nc",
                     "tas_day_MRI-CGCM3_rcp85_r1i1p1_20760101-20851231.nc",
                     "tas_day_MRI-CGCM3_rcp85_r1i1p1_20860101-20951231.nc",
                     "tas_day_MRI-CGCM3_rcp85_r1i1p1_20960101-21001231.nc")
MRI_CGCM3_11 <- list(MRI_CGCM3_hist, MRI_CGCM3_rcp85)

gcm_files <- list(CMCC_CMS_01, GFDL_CM3_02, GFDL_ESM2G_03, HadGEM2_ES_04,
                  inmcm4_05, IPSL_CM5A_MR_06, MIROC_ESM_CHEM_07, MIROC5_08,
                  MPI_ESM_LR_09, MPI_ESM_MR_10, MRI_CGCM3_11)

###################################################################
###     calling functions on spectral exponent files             ## 
###################################################################
p = folders[i]
gcm = gcm_models[i]

spec_exp <- create_rasterStack(p = p, gcm = gcm)
calculate_avg_change(p = p, gcm = gcm)

write.csv(spec_exp, paste(p, gcm, "_spec_exp.csv",  sep = ""), row.names = F)
          
          
          
