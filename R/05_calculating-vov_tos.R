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
create_rasterStack <- function(path, type) {
 
  ## read in spatial chunk file names:
  if (type == "GCM") {
    se_filenames <- readRDS(paste(path, "se_filenames_tos.rds",  sep = ""))
  }
  else if (type == "ERA40") {
    se_filenames <- readRDS("data-processed/ERA-40/era-40_sst_sp_files.rds")
  }

  #se_filenames <- str_replace_all(se_filenames, "data-raw/", "/Volumes/SundayLab/CMIP5-GCMs_tos/")
  
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
      
      ## create raster layer:
      l_layer_PSD_low <- rasterFromXYZ(l_sp_df_PSD_low)
      s_layer_PSD_low <- rasterFromXYZ(s_sp_df_PSD_low)
      l_layer_AWC <- rasterFromXYZ(l_sp_df_AWC)
      s_layer_AWC <- rasterFromXYZ(s_sp_df_AWC)
      l_layer_PSD_high <- rasterFromXYZ(l_sp_df_PSD_high)
      s_layer_PSD_high <- rasterFromXYZ(s_sp_df_PSD_high)
      #plot(l_layer_PSD_high)
      
      ## add to temporary rasterstack:
      if (step == length(step_split)) {
        l_temp_stack_PSD_low <- l_layer_PSD_low
        s_temp_stack_PSD_low <- s_layer_PSD_low
        l_temp_stack_AWC <- l_layer_AWC
        s_temp_stack_AWC <- s_layer_AWC
        l_temp_stack_PSD_high <- l_layer_PSD_high
        s_temp_stack_PSD_high <- s_layer_PSD_high
      }
      else {
        l_temp_stack_PSD_low <- addLayer(l_layer_PSD_low, l_temp_stack_PSD_low)
        s_temp_stack_PSD_low <- addLayer(s_layer_PSD_low, s_temp_stack_PSD_low)
        l_temp_stack_AWC <- addLayer(l_layer_AWC, l_temp_stack_AWC)
        s_temp_stack_AWC <- addLayer(s_layer_AWC, s_temp_stack_AWC)
        l_temp_stack_PSD_high <- addLayer(l_layer_PSD_high, l_temp_stack_PSD_high)
        s_temp_stack_PSD_high <- addLayer(s_layer_PSD_high, s_temp_stack_PSD_high)
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
    }
    else {
      l_stack_list_PSD_low <- append(l_stack_list_PSD_low, l_temp_stack_PSD_low)
      s_stack_list_PSD_low <- append(s_stack_list_PSD_low, s_temp_stack_PSD_low)
      l_stack_list_AWC <- append(l_stack_list_AWC, l_temp_stack_AWC)
      s_stack_list_AWC <- append(s_stack_list_AWC, s_temp_stack_AWC)
      l_stack_list_PSD_high <- append(l_stack_list_PSD_high, l_temp_stack_PSD_high)
      s_stack_list_PSD_high <- append(s_stack_list_PSD_high, s_temp_stack_PSD_high)
    }
    
    ## move to next time window width
    i = i + 1
  }
  
  ## name the list items 
  names(l_stack_list_PSD_low) <- names(l_stack_list_PSD_high) <- names(s_stack_list_PSD_low) <-
    names(s_stack_list_PSD_high) <- names(l_stack_list_AWC) <- names(s_stack_list_PSD_high) <-
    names(ww_split)
  
  ## save the rasterstack 
  saveRDS(l_stack_list_PSD_low, paste(path, type, "_l_stack_list_PSD_low_tos.rds", sep = ""))
  saveRDS(s_stack_list_PSD_low, paste(path, type, "_s_stack_list_PSD_low_tos.rds", sep = ""))
  saveRDS(l_stack_list_AWC, paste(path, type, "_l_stack_list_AWC_tos.rds", sep = ""))
  saveRDS(s_stack_list_AWC, paste(path, type, "_s_stack_list_AWC_tos.rds", sep = ""))
  saveRDS(l_stack_list_PSD_high, paste(path, type, "_l_stack_list_PSD_high_tos.rds", sep = ""))
  saveRDS(s_stack_list_PSD_high, paste(path, type, "s_stack_list_PSD_high_tos.rds", sep = ""))
  
  stacks <- list(l_stack_list_PSD_low, s_stack_list_PSD_low, 
                 l_stack_list_AWC, s_stack_list_AWC,
                 l_stack_list_PSD_high, s_stack_list_PSD_high)
  
  ## return the 6 lists of rasterStacks
  return(stacks)
}

#################################################
###                setting paths               ## 
#################################################
## set 'path' to where you have the GCM files stored on your computer
## for me, they are here:
#path = "/Volumes/SundayLab/CMIP5-GCMs_tos/" ## change me
#path = "data-raw/"
path = "CMIP5-GCMs_tos/"

## create vector of file folders to put data into:
gcm_models <- c("01_CMCC-CMS_tos",
                "02_GFDL-CM3_tos",
                "03_GFDL-ESM2G_tos",
                "04_HadGEM2-ES_tos",
                "05_inmcm4_tos",
                "06_IPSL-CM5A-MR_tos",
                "07_MIROC-ESM-CHEM_tos",
                "08_MIROC5_tos",
                "09_MPI-ESM-LR_tos",
                "10_MPI-ESM-MR_tos",
                "11_MRI-CGCM3_tos")

folders <- paste(path, gcm_models, "/", sep = "")

#################################################
###       make list of GCM file names          ## 
#################################################
CMCC_CMS_hist <- c("tos_day_CMCC-CMS_historical_r1i1p1_18700101-18791231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_18800101-18891231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_18900101-18991231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19000101-19091231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19100101-19191231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19200101-19291231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19300101-19391231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19400101-19491231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19500101-19591231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19600101-19691231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19700101-19791231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19800101-19891231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_19900101-19991231.nc",
                   "tos_day_CMCC-CMS_historical_r1i1p1_20000101-20051231.nc")
CMCC_CMS_rcp85 <- c("tos_day_CMCC-CMS_rcp85_r1i1p1_20060101-20091231.nc",
                    "tos_day_CMCC-CMS_rcp85_r1i1p1_20100101-20191231.nc",
                    "tos_day_CMCC-CMS_rcp85_r1i1p1_20200101-20291231.nc",
                    "tos_day_CMCC-CMS_rcp85_r1i1p1_20300101-20391231.nc",
                    "tos_day_CMCC-CMS_rcp85_r1i1p1_20400101-20491231.nc",
                    "tos_day_CMCC-CMS_rcp85_r1i1p1_20500101-20591231.nc",
                    "tos_day_CMCC-CMS_rcp85_r1i1p1_20600101-20691231.nc",
                    "tos_day_CMCC-CMS_rcp85_r1i1p1_20700101-20791231.nc",
                    "tos_day_CMCC-CMS_rcp85_r1i1p1_20800101-20891231.nc",
                    "tos_day_CMCC-CMS_rcp85_r1i1p1_20900101-21001231.nc")
CMCC_CMS_01 <- list(CMCC_CMS_hist, CMCC_CMS_rcp85)

GFDL_CM3_hist <- c("tos_day_GFDL-CM3_historical_r1i1p1_18700101-18741231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_18750101-18791231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_18800101-18841231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_18850101-18891231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_18900101-18941231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_18950101-18991231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19000101-19041231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19050101-19091231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19100101-19141231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19150101-19191231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19200101-19241231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19250101-19291231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19300101-19341231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19350101-19391231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19400101-19441231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19450101-19491231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19500101-19541231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19550101-19591231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19600101-19641231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19650101-19691231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19700101-19741231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19750101-19791231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19800101-19841231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19850101-19891231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19900101-19941231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_19950101-19991231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_20000101-20041231.nc",
                   "tos_day_GFDL-CM3_historical_r1i1p1_20050101-20051231.nc")
GFDL_CM3_rcp85 <- c("tos_day_GFDL-CM3_rcp85_r1i1p1_20060101-20101231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20110101-20151231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20160101-20201231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20210101-20251231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20260101-20301231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20310101-20351231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20360101-20401231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20410101-20451231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20460101-20501231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20510101-20551231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20560101-20601231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20610101-20651231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20660101-20701231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20710101-20751231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20760101-20801231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20810101-20851231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20860101-20901231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20910101-20951231.nc",
                    "tos_day_GFDL-CM3_rcp85_r1i1p1_20960101-21001231.nc")
GFDL_CM3_02 <- list(GFDL_CM3_hist, GFDL_CM3_rcp85)

GFDL_ESM2G_hist <- c("tos_day_GFDL-ESM2G_historical_r1i1p1_18660101-18701231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_18710101-18751231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_18760101-18801231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_18810101-18851231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_18860101-18901231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_18910101-18951231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_18960101-19001231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19010101-19051231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19060101-19101231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19110101-19151231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19160101-19201231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19210101-19251231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19260101-19301231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19310101-19351231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19360101-19401231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19410101-19451231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19460101-19501231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19510101-19551231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19560101-19601231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19610101-19651231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19660101-19701231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19710101-19751231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19760101-19801231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19810101-19851231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19860101-19901231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19910101-19951231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_19960101-20001231.nc",
                     "tos_day_GFDL-ESM2G_historical_r1i1p1_20010101-20051231.nc")
GFDL_ESM2G_rcp85 <- c("tos_day_GFDL-ESM2G_rcp85_r1i1p1_20060101-20101231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20110101-20151231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20160101-20201231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20210101-20251231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20260101-20301231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20310101-20351231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20360101-20401231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20410101-20451231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20460101-20501231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20510101-20551231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20560101-20601231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20610101-20651231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20660101-20701231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20710101-20751231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20760101-20801231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20810101-20851231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20860101-20901231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20910101-20951231.nc",
                      "tos_day_GFDL-ESM2G_rcp85_r1i1p1_20960101-21001231.nc")
GFDL_ESM2G_03 <- list(GFDL_ESM2G_hist, GFDL_ESM2G_rcp85)

HadGEM2_ES_hist <- c("tos_day_HadGEM2-ES_historical_r1i1p1_18691201-18791130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_18791201-18891130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_18891201-18991130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_18991201-19091130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19091201-19191130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19191201-19291130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19291201-19391130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19391201-19491130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19491201-19591130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19591201-19691130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19691201-19791130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19791201-19891130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19891201-19991130.nc",
                     "tos_day_HadGEM2-ES_historical_r1i1p1_19991201-20051130.nc")
HadGEM2_ES_rcp85 <- c("tos_day_HadGEM2-ES_rcp85_r1i1p1_20051201-20101130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20101201-20151130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20151201-20201130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20201201-20251130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20251201-20301130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20301201-20351130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20351201-20401130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20401201-20451130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20451201-20501130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20501201-20551130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20551201-20601130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20601201-20651130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20651201-20701130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20701201-20751130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20751201-20851130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20851201-20951130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20951201-20991130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20991201-20991230.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20991201-21091130.nc")
HadGEM2_ES_04 <- list(HadGEM2_ES_hist, HadGEM2_ES_rcp85)

inmcm4_hist <- c("tos_day_inmcm4_historical_r1i1p1_18700101-18791231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_18800101-18891231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_18900101-18991231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19000101-19091231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19100101-19191231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19200101-19291231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19300101-19391231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19400101-19491231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19500101-19591231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19600101-19691231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19700101-19791231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19800101-19891231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_19900101-19991231.nc",
                 "tos_day_inmcm4_historical_r1i1p1_20000101-20051231.nc")
inmcm4_rcp85 <- c("tos_day_inmcm4_rcp85_r1i1p1_20060101-20151231.nc",
                  "tos_day_inmcm4_rcp85_r1i1p1_20160101-20251231.nc",
                  "tos_day_inmcm4_rcp85_r1i1p1_20260101-20351231.nc",
                  "tos_day_inmcm4_rcp85_r1i1p1_20360101-20451231.nc",
                  "tos_day_inmcm4_rcp85_r1i1p1_20460101-20551231.nc",
                  "tos_day_inmcm4_rcp85_r1i1p1_20560101-20651231.nc",
                  "tos_day_inmcm4_rcp85_r1i1p1_20660101-20751231.nc",
                  "tos_day_inmcm4_rcp85_r1i1p1_20760101-20851231.nc",
                  "tos_day_inmcm4_rcp85_r1i1p1_20860101-20951231.nc",
                  "tos_day_inmcm4_rcp85_r1i1p1_20960101-21001231.nc")
inmcm4_05 <- list(inmcm4_hist, inmcm4_rcp85)

IPSL_CM5A_MR_hist <- c("tos_day_IPSL-CM5A-MR_historical_r1i1p1_18500101-18991231.nc",
                       "tos_day_IPSL-CM5A-MR_historical_r1i1p1_19000101-19491231.nc",
                       "tos_day_IPSL-CM5A-MR_historical_r1i1p1_19500101-19991231.nc",
                       "tos_day_IPSL-CM5A-MR_historical_r1i1p1_20000101-20051231.nc")
IPSL_CM5A_MR_rcp85 <- c("tos_day_IPSL-CM5A-MR_rcp85_r1i1p1_20060101-20551231.nc",
                        "tos_day_IPSL-CM5A-MR_rcp85_r1i1p1_20560101-21001231.nc")
IPSL_CM5A_MR_06 <- list(IPSL_CM5A_MR_hist, IPSL_CM5A_MR_rcp85)

MIROC_ESM_CHEM_hist <- c('tos_day_MIROC-ESM-CHEM_historical_r1i1p1_18500101-20051231.nc') 
MIROC_ESM_CHEM_rcp85 <- c("tos_day_MIROC-ESM-CHEM_rcp85_r1i1p1_20060101-21001231.nc")
MIROC_ESM_CHEM_07 <- list(MIROC_ESM_CHEM_hist, MIROC_ESM_CHEM_rcp85)

MIROC5_hist <- c("tos_day_MIROC5_historical_r1i1p1_18700101-18791231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_18800101-18891231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_18900101-18991231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19000101-19091231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19100101-19191231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19200101-19291231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19300101-19391231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19400101-19491231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19500101-19591231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19600101-19691231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19700101-19791231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19800101-19891231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_19900101-19991231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_20000101-20091231.nc",
                 "tos_day_MIROC5_historical_r1i1p1_20100101-20121231.nc")
MIROC5_rcp85 <- c("tos_day_MIROC5_rcp85_r1i1p1_20060101-20091231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_20100101-20191231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_20200101-20291231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_20300101-20391231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_20400101-20491231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_20500101-20591231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_20600101-20691231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_20700101-20791231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_20800101-20891231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_20900101-20991231.nc",
                  "tos_day_MIROC5_rcp85_r1i1p1_21000101-21001231.nc")
MIROC5_08 <- list(MIROC5_hist, MIROC5_rcp85)

MPI_ESM_LR_hist <- c("tos_day_MPI-ESM-LR_historical_r1i1p1_18700101-18791231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_18800101-18891231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_18900101-18991231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19000101-19091231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19100101-19191231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19200101-19291231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19300101-19391231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19400101-19491231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19500101-19591231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19600101-19691231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19700101-19791231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19800101-19891231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_19900101-19991231.nc",
                     "tos_day_MPI-ESM-LR_historical_r1i1p1_20000101-20051231.nc")
MPI_ESM_LR_rcp85 <- c("tos_day_MPI-ESM-LR_rcp85_r1i1p1_20060101-20091231.nc",
                      "tos_day_MPI-ESM-LR_rcp85_r1i1p1_20100101-20191231.nc",
                      "tos_day_MPI-ESM-LR_rcp85_r1i1p1_20200101-20291231.nc",
                      "tos_day_MPI-ESM-LR_rcp85_r1i1p1_20300101-20391231.nc",
                      "tos_day_MPI-ESM-LR_rcp85_r1i1p1_20400101-20491231.nc",
                      "tos_day_MPI-ESM-LR_rcp85_r1i1p1_20500101-20591231.nc",
                      "tos_day_MPI-ESM-LR_rcp85_r1i1p1_20600101-20691231.nc",
                      "tos_day_MPI-ESM-LR_rcp85_r1i1p1_20700101-20791231.nc",
                      "tos_day_MPI-ESM-LR_rcp85_r1i1p1_20800101-20891231.nc",
                      "tos_day_MPI-ESM-LR_rcp85_r1i1p1_20900101-21001231.nc")
MPI_ESM_LR_09 <- list(MPI_ESM_LR_hist, MPI_ESM_LR_rcp85)

MPI_ESM_MR_hist <- c("tos_day_MPI-ESM-MR_historical_r1i1p1_18700101-18791231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_18800101-18891231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_18900101-18991231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19000101-19091231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19100101-19191231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19200101-19291231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19300101-19391231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19400101-19491231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19500101-19591231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19600101-19691231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19700101-19791231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19800101-19891231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_19900101-19991231.nc",
                     "tos_day_MPI-ESM-MR_historical_r1i1p1_20000101-20051231.nc")
MPI_ESM_MR_rcp85 <- c("tos_day_MPI-ESM-MR_rcp85_r1i1p1_20060101-20091231.nc",
                      "tos_day_MPI-ESM-MR_rcp85_r1i1p1_20100101-20191231.nc",
                      "tos_day_MPI-ESM-MR_rcp85_r1i1p1_20200101-20291231.nc",
                      "tos_day_MPI-ESM-MR_rcp85_r1i1p1_20300101-20391231.nc",
                      "tos_day_MPI-ESM-MR_rcp85_r1i1p1_20400101-20491231.nc",
                      "tos_day_MPI-ESM-MR_rcp85_r1i1p1_20500101-20591231.nc",
                      "tos_day_MPI-ESM-MR_rcp85_r1i1p1_20600101-20691231.nc",
                      "tos_day_MPI-ESM-MR_rcp85_r1i1p1_20700101-20791231.nc",
                      "tos_day_MPI-ESM-MR_rcp85_r1i1p1_20800101-20891231.nc",
                      "tos_day_MPI-ESM-MR_rcp85_r1i1p1_20900101-21001231.nc")
MPI_ESM_MR_10 <- list(MPI_ESM_MR_hist, MPI_ESM_MR_rcp85)

MRI_CGCM3_hist <- c("tos_day_MRI-CGCM3_historical_r1i1p1_18700101-18791231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_18800101-18891231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_18900101-18991231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19000101-19091231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19100101-19191231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19200101-19291231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19300101-19391231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19400101-19491231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19500101-19591231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19600101-19691231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19700101-19791231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19800101-19891231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_19900101-19991231.nc",
                    "tos_day_MRI-CGCM3_historical_r1i1p1_20000101-20051231.nc")
MRI_CGCM3_rcp85 <- c("tos_day_MRI-CGCM3_rcp85_r1i1p1_20060101-20151231.nc",
                     "tos_day_MRI-CGCM3_rcp85_r1i1p1_20160101-20251231.nc",
                     "tos_day_MRI-CGCM3_rcp85_r1i1p1_20260101-20351231.nc",
                     "tos_day_MRI-CGCM3_rcp85_r1i1p1_20360101-20451231.nc",
                     "tos_day_MRI-CGCM3_rcp85_r1i1p1_20460101-20551231.nc",
                     "tos_day_MRI-CGCM3_rcp85_r1i1p1_20560101-20651231.nc",
                     "tos_day_MRI-CGCM3_rcp85_r1i1p1_20660101-20751231.nc",
                     "tos_day_MRI-CGCM3_rcp85_r1i1p1_20760101-20851231.nc",
                     "tos_day_MRI-CGCM3_rcp85_r1i1p1_20860101-20951231.nc",
                     "tos_day_MRI-CGCM3_rcp85_r1i1p1_20960101-21001231.nc")
MRI_CGCM3_11 <- list(MRI_CGCM3_hist, MRI_CGCM3_rcp85)

gcm_files <- list(CMCC_CMS_01, GFDL_CM3_02, GFDL_ESM2G_03, HadGEM2_ES_04,
                  inmcm4_05, IPSL_CM5A_MR_06, MIROC_ESM_CHEM_07, MIROC5_08,
                  MPI_ESM_LR_09, MPI_ESM_MR_10, MRI_CGCM3_11)

###################################################################
###     calling functions on spectral exponent files             ## 
###################################################################
path = folders[i]
