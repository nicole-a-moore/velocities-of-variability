### reorganize GCMs into spatial chunks, remove for outliers and detrend time series
### R version 3.6.1 
.libPaths(c("~/projects/def-jsunday/nikkim/VoV/packages", .libPaths()))
library(tidyverse)
library(raster)
library(easyNCDF)
library(lubridate)
library(ncdf4)

#################################################
###                   FUNCTIONS                ## 
#################################################
reorganize_GCM <- function(historical_filenames, rcp85_filenames, path, gcm_num) {
  og_names = paste(path, append(historical_filenames, rcp85_filenames), sep = "")
  
  ## get paths and filenames of resampled GCMs
  if(str_detect(path, "MIROC-ESM-CHEM") | str_detect(path, "HadGEM")) {
    filenames <-  paste("resampled_", append(historical_filenames, rcp85_filenames), sep = "")
  }
  else {
    filenames <-  paste("regridded_", append(historical_filenames, rcp85_filenames), sep = "")
  }
  
  paths <- paste(path, filenames, sep = "")
  
  #####       REORGANIZE INTO SPATIAL CHUNKS      #####
  ## to detrend and perform sliding spectral analysis, need the whole time series for each location at once
  ## to satisfy this requirement while avoiding memory exhaustion, reorganize data
  ## go from time chunks spanning all of space to spatial chunks spanning all of time
  ## break into 8 x 60 degree lat x 60 degree lon chunks with data from 1850-2100
  lon_index = 0
  lat_index = 90
  count = 1
  sp_files <- c()
  while (count <= 18) {
    ## get lat and lon bounds to extract in between:
    lon_bound1 <- lon_index 
    lon_bound2 <- lon_index + 60
    lat_bound1 <- lat_index
    lat_bound2 <- lat_index - 60
    
    ## loop through resampled GCM files, extracting data within spatial chunk and appending 
    file = 1
    dates <- c()
    while (file < length(filenames)+1) {
      start <- Sys.time()
      temps <- stack(paths[file])
      
      ## get dates
      time = ncvar_get(nc_open(og_names[file]), "time") 
      
      ## convert 
      ## days since 1870-1-1
      if(gcm_num == 1) {
        dates <- append(dates, dates <- as.Date(time, origin = '1870-01-01'))
      }
      else {
        ## days since 1850-1-1
        ## 07
        dates <- append(dates, dates <- as.Date(time, origin = '1850-01-01'))
      }
      
      ## crop to new extent within bounds of spatial chunk
      extent <- extent(lon_bound1, lon_bound2, lat_bound2, lat_bound1)
      cropped <- crop(temps, extent)
      
      ## add cropped raster to rest from the GCM:
      if(file == 1) {
        spatial_temps <- cropped
      }
      else {
        spatial_temps <- stack(spatial_temps, cropped)
      }

      ## estimating time left
      end <- Sys.time()
      elapsed <- end - start
      expected <- elapsed[[1]] * (length(filenames)-file) 
      print(paste0("remaining: ", round(expected,3), " minutes = ", round(expected/60,3), " hours"))
      file = file + 1
    }

    ## turn into an array:
    temps_df <- as.array(spatial_temps)  
    
    ## figure out which dates to remove:
    b4_1871_af_2101 <- which(as.numeric(str_split_fixed(dates, "\\-",n=2)[,1]) <= 1870 |
                               as.numeric(str_split_fixed(dates, "\\-",n=2)[,1]) >= 2101)
    ly_days <- which(str_detect(dates, "02-29"))
    remove <- c(b4_1871_af_2101, ly_days)
    date_new <- dates[-remove]
    
    ## remove dates:
    temps_df <- temps_df[,,-remove] 
    
    ## loop through cells in spatial chunk, removing outliers and interpolating missing vals in time series
    l_detrended_temps <- s_detrended_temps <- temps_df
    x = 1 ## longitude
    while (x < ncol(temps_df) + 1) {
      y = 1 ## latitude
      while (y < nrow(temps_df)+1) {
        ## get local time series 
        local_ts <- data.frame(time = 1:(length(temps_df[1,1,])), 
                               temp = temps_df[y,x,])
        
        if (!length(which(is.na(local_ts$temp))) == nrow(local_ts)) {
          #####     REMOVING OUTLIERS, INTERPOLATING MISSING TEMPS    #####
          ## temperature remove values > 60C (333.15K)
          local_ts$temp[which(local_ts$temp > 333.15)] <- NA
          
          ## interpolate if temps are missing
          if(length(which(is.na(local_ts$temp))) != 0 & length(which(is.na(local_ts$temp))) 
             != length(local_ts$temp)) {
            local_ts$temp <- na_kalman(local_ts$temp, smooth = TRUE, model = "StructTS")
            temps_df[y,x,] <- local_ts$temp ## save changed time series if it needed changing
          }
          
          print(paste("Done removing outliers & interpolating x ",x,  
                      " y ", y, " of chunk #", count,sep = ""))
          
          #####     LINEARLY AND SEASONALLY DETREND TIME SERIES IN EACH RASTER CELL    #####
          ## create empty objects
          local_ts$md <- str_split_fixed(date_new, "\\-", n=2)[,2]
          
          ts_df <- local_ts %>%
            group_by(md) %>%
            do(mutate(., temp_profile = mean(.$temp))) %>% ## compute temp climatology for each day of year
            ungroup() %>%
            mutate(s_detrended_temp = temp - temp_profile) %>% ## create column representing seasonally detrended
            arrange(., time)
          
          ## run linear regression for grid cell
          l_output <- lm(ts_df, formula = temp ~ time)
          s_output <- lm(ts_df, formula = s_detrended_temp ~ time)
          
          ## extract residuals and add to detrended temps objects:
          l_detrended_temps[y,x,] <- l_output$residuals
          s_detrended_temps[y,x,] <- s_output$residuals
          
          print(paste("Done detrending x ", x,  " y ", y, " of chunk #", count,sep = ""))
        }
        y = y + 1
      }
      x = x + 1
    }
    
    ## save:
    sp_files[count] <- paste(path, "spatial_temps_lon-", lon_bound1,"-", lon_bound2,
                             "_lat-", lat_bound1, "-", lat_bound2,"_tos.nc", sep = "")
    ArrayToNc(temps_df, file_path = paste(path, "spatial_temps_lon-", lon_bound1,"-", lon_bound2,
                                          "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(l_detrended_temps, file_path = paste(path, "l-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,"_tos.nc", sep = ""))
    ArrayToNc(s_detrended_temps, file_path = paste(path, "s-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,"_tos.nc", sep = ""))
    
    ## advance lat and lon indecies to move to next spatial chunk
    if (count %in% c(6, 12)) {
      lat_index <- lat_index - 60
      lon_index <- 0
    }
    else {
      lon_index <- lon_index + 60
    }
    
    ## move to next spatial chunk 
    count = count + 1
  }
  
  saveRDS(date_new, paste(path, "date_new_tos.rds", sep = ""))
  saveRDS(sp_files, paste(path, "sp_files_tos.rds", sep = ""))
  ## returns list of spatial chunk filenames 
  return(sp_files)
}




#################################################
###                setting paths               ## 
#################################################
## set 'path' to where you have the GCM files stored on your computer
## for me, they are here:
#path = "/Volumes/NIKKI/CMIP5-GCMs/" ## change me
#path = "data-raw/"
path = "CMIP5-GCMs/"

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


MIROC_ESM_CHEM_hist <- c("tos_day_MIROC-ESM-CHEM_historical_r1i1p1_18700101-18891231.nc",
                         "tos_day_MIROC-ESM-CHEM_historical_r1i1p1_18900101-19091231.nc",
                         "tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19100101-19291231.nc",
                         "tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19300101-19491231.nc",
                         "tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19500101-19691231.nc",
                         "tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19700101-19891231.nc",
                         "tos_day_MIROC-ESM-CHEM_historical_r1i1p1_19900101-20051231.nc") 
MIROC_ESM_CHEM_rcp85 <- c("tos_day_MIROC-ESM-CHEM_rcp85_r1i1p1_20060101-20251231.nc",
                          "tos_day_MIROC-ESM-CHEM_rcp85_r1i1p1_20260101-20451231.nc",
                          "tos_day_MIROC-ESM-CHEM_rcp85_r1i1p1_20460101-20651231.nc",
                          "tos_day_MIROC-ESM-CHEM_rcp85_r1i1p1_20660101-20851231.nc",
                          "tos_day_MIROC-ESM-CHEM_rcp85_r1i1p1_20860101-21001231.nc")
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


#################################################
###     calling functions on a GCM             ## 
#################################################
## read the command line arguments and call the functions on the corresponding GCM 

command_args <- commandArgs(trailingOnly = TRUE)
i = as.numeric(command_args[1])

element = gcm_files[[i]]
hfn <- element[[1]]
rfn <- element[[2]]
p <- folders[i]
gcm_num = 1

reorganize = reorganize_GCM(historical_filenames = hfn, rcp85_filenames = rfn, path = p, gcm_num = gcm_num)

