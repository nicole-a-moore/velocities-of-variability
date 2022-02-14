### reorganize GCMs into spatial chunks, remove for outliers and detrend time series
### R version 3.6.1 
.libPaths(c("~/projects/def-jsunday/nikkim/VoV/packages", .libPaths()))
library(tidyverse)
library(raster)
library(easyNCDF)

#################################################
###                   FUNCTIONS                ## 
#################################################
reorganize_GCM <- function(historical_filenames, rcp85_filenames, path, dates) {
  
  ## get paths and filenames of resampled GCMs
  filenames <- paste("resampled_", append(historical_filenames, rcp85_filenames), sep = "")
  paths <- paste(path, filenames, sep = "")
  
  ## get dates 
  dates <- readRDS(paste(path, "dates.rds", sep = ""))
  ## figure out which dates to remove:
  b4_1871_af_2101 <- which(as.numeric(str_split_fixed(dates, "\\.",n=2)[,1]) <= 1870 |
                             as.numeric(str_split_fixed(dates, "\\.",n=2)[,1]) >= 2101)
  ly_days <- which(str_detect(dates, "02.29"))
  remove <- c(b4_1871_af_2101, ly_days)
  date_new <- dates[-remove]
  
  
  
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
    while (file < length(filenames)+1) {
      start <- Sys.time()
      temps <- stack(paths[file])
      
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
        
        ggplot(local_ts, aes(x = time, y = temp)) + geom_line() +
          geom_smooth(method = "lm")
        
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
          local_ts$md <- str_split_fixed(date_new, "\\.", n=2)[,2]
          
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
                             "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = "")
    # ArrayToNc(temps_df, file_path = paste(path, "spatial_temps_lon-", lon_bound1,"-", lon_bound2,
    #                                       "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(l_detrended_temps, file_path = paste(path, "l-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(s_detrended_temps, file_path = paste(path, "s-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    
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
  
  saveRDS(date_new, paste(path, "date_new.rds", sep = ""))
  saveRDS(sp_files, paste(path, "sp_files.rds", sep = ""))
  ## returns list of spatial chunk filenames 
  return(sp_files)
}




#################################################
###                setting paths               ## 
#################################################
## set 'path' to where you have the GCM files stored on your computer
## for me, they are here:
#path = "/Volumes/SundayLab/CMIP5-GCMs/" ## change me
path = "CMIP5-GCMs/"

## create vector of file folders to put data into:
gcm_models <- c("01_CMCC-CESM", "02_CMCC-CM", '03_CMCC-CMS', '04_MPI-ESM-LR', '05_MPI-ESM-MR',
                "06_GFDL-ESM2G", '07_GFDL-CM3', '08_GFDL-ESM2M', '09_HadGEM2-CC', '10_HadGEM2-ES',
                "11_HadGEM2-AO", '12_IPSL-CM5A-LR', '13_IPSL-CM5B-LR', '14_MIROC5', '15_MIROC5-ESM-CHEM',
                '16_MIROC5-ESM', "17_inmcm4", '18_CNRM-CM5', "19_MRI-CGCM3", '20_MRI-ESM1',
                '21_IPSL-CM5A-MR')

folders <- paste(path, gcm_models, "/", sep = "")

#################################################
###       make list of GCM file names          ## 
#################################################
CMCC_CESM_hist <- c("tas_day_CMCC-CESM_historical_r1i1p1_18700101-18741231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_18750101-18791231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_18800101-18841231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_18850101-18891231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_18900101-18941231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_18950101-18991231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19000101-19041231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19050101-19091231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19100101-19141231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19150101-19191231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19200101-19241231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19250101-19291231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19300101-19341231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19350101-19391231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19400101-19441231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19450101-19491231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19500101-19541231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19550101-19591231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19600101-19641231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19650101-19691231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19700101-19741231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19750101-19791231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19800101-19841231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19850101-19891231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19900101-19941231.nc",
                    "tas_day_CMCC-CESM_historical_r1i1p1_19950101-19991231.nc")
CMCC_CESM_rcp85 <- c("tas_day_CMCC-CESM_rcp85_r1i1p1_20000101-20041231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20050101-20051231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20060101-20151231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20160101-20251231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20260101-20351231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20360101-20451231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20460101-20551231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20560101-20651231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20660101-20751231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20760101-20851231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20860101-20951231.nc",
                     "tas_day_CMCC-CESM_rcp85_r1i1p1_20960101-21001231.nc")
CMCC_CESM_01 <- list(CMCC_CESM_hist, CMCC_CESM_rcp85)

CMCC_CM_hist <- c("tas_day_CMCC-CM_historical_r1i1p1_18700101-18701231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18710101-18711231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18720101-18721231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18730101-18731231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18740101-18741231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18750101-18751231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18760101-18761231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18770101-18771231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18780101-18781231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18790101-18791231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18800101-18801231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18810101-18811231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18820101-18821231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18830101-18831231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18840101-18841231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18850101-18851231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18860101-18861231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18870101-18871231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18880101-18881231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18890101-18891231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18900101-18901231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18910101-18911231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18920101-18921231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18930101-18931231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18940101-18941231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18950101-18951231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18960101-18961231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18970101-18971231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18980101-18981231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_18990101-18991231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19000101-19001231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19010101-19011231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19020101-19021231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19030101-19031231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19040101-19041231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19050101-19051231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19060101-19061231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19070101-19071231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19080101-19081231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19090101-19091231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19100101-19101231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19110101-19111231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19120101-19121231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19130101-19131231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19140101-19141231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19150101-19151231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19160101-19161231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19170101-19171231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19180101-19181231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19190101-19191231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19200101-19201231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19210101-19211231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19220101-19221231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19230101-19231231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19240101-19241231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19250101-19251231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19260101-19261231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19270101-19271231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19280101-19281231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19290101-19291231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19300101-19301231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19310101-19311231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19320101-19321231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19330101-19331231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19340101-19341231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19350101-19351231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19360101-19361231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19370101-19371231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19380101-19381231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19390101-19391231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19400101-19401231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19410101-19411231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19420101-19421231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19430101-19431231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19440101-19441231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19450101-19451231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19460101-19461231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19470101-19471231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19480101-19481231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19490101-19491231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19500101-19501231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19510101-19511231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19520101-19521231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19530101-19531231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19540101-19541231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19550101-19551231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19560101-19561231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19570101-19571231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19580101-19581231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19590101-19591231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19600101-19601231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19610101-19611231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19620101-19621231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19630101-19631231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19640101-19641231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19650101-19651231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19660101-19661231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19670101-19671231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19680101-19681231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19690101-19691231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19700101-19701231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19710101-19711231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19720101-19721231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19730101-19731231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19740101-19741231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19750101-19751231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19760101-19761231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19770101-19771231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19780101-19781231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19790101-19791231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19800101-19801231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19810101-19811231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19820101-19821231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19830101-19831231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19840101-19841231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19850101-19851231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19860101-19861231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19870101-19871231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19880101-19881231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19890101-19891231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19900101-19901231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19910101-19911231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19920101-19921231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19930101-19931231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19940101-19941231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19950101-19951231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19960101-19961231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19970101-19971231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19980101-19981231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_19990101-19991231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_20000101-20001231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_20010101-20011231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_20020101-20021231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_20030101-20031231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_20040101-20041231.nc",
                  "tas_day_CMCC-CM_historical_r1i1p1_20050101-20051231.nc")
CMCC_CM_rcp85 <- c("tas_day_CMCC-CM_rcp85_r1i1p1_20060101-20061231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20070101-20071231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20080101-20081231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20090101-20091231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20100101-20101231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20110101-20111231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20120101-20121231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20130101-20131231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20140101-20141231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20150101-20151231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20160101-20161231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20170101-20171231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20180101-20181231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20190101-20191231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20200101-20201231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20210101-20211231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20220101-20221231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20230101-20231231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20240101-20241231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20250101-20251231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20260101-20261231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20270101-20271231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20280101-20281231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20290101-20291231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20300101-20301231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20310101-20311231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20320101-20321231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20330101-20331231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20340101-20341231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20350101-20351231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20360101-20361231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20370101-20371231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20380101-20381231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20390101-20391231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20400101-20401231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20410101-20411231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20420101-20421231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20430101-20431231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20440101-20441231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20450101-20451231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20460101-20461231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20470101-20471231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20480101-20481231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20490101-20491231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20500101-20501231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20510101-20511231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20520101-20521231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20530101-20531231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20540101-20541231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20550101-20551231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20560101-20561231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20570101-20571231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20580101-20581231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20590101-20591231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20600101-20601231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20610101-20611231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20620101-20621231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20630101-20631231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20640101-20641231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20650101-20651231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20660101-20661231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20670101-20671231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20680101-20681231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20690101-20691231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20700101-20701231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20710101-20711231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20720101-20721231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20730101-20731231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20740101-20741231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20750101-20751231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20760101-20761231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20770101-20771231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20780101-20781231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20790101-20791231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20800101-20801231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20810101-20811231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20820101-20821231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20830101-20831231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20840101-20841231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20850101-20851231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20860101-20861231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20870101-20871231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20880101-20881231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20890101-20891231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20900101-20901231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20910101-20911231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20920101-20921231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20930101-20931231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20940101-20941231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20950101-20951231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20960101-20961231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20970101-20971231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20980101-20981231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_20990101-20991231.nc",
                   "tas_day_CMCC-CM_rcp85_r1i1p1_21000101-21001231.nc")
CMCC_CM_02 <- list(CMCC_CM_hist, CMCC_CM_rcp85)

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
CMCC_CMS_03 <- list(CMCC_CMS_hist, CMCC_CMS_rcp85)

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
MPI_ESM_LR_04 <- list(MPI_ESM_LR_hist, MPI_ESM_LR_rcp85)

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
MPI_ESM_MR_05 <- list(MPI_ESM_MR_hist, MPI_ESM_MR_rcp85)

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
GFDL_ESM2G_06 <- list(GFDL_ESM2G_hist, GFDL_ESM2G_rcp85)

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
GFDL_CM3_07 <- list(GFDL_CM3_hist, GFDL_CM3_rcp85)

GFDL_ESM2M_hist <- c("tas_day_GFDL-ESM2M_historical_r1i1p1_18660101-18701231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_18710101-18751231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_18760101-18801231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_18810101-18851231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_18860101-18901231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_18910101-18951231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_18960101-19001231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19010101-19051231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19060101-19101231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19110101-19151231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19160101-19201231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19210101-19251231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19260101-19301231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19310101-19351231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19360101-19401231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19410101-19451231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19460101-19501231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19510101-19551231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19560101-19601231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19610101-19651231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19660101-19701231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19710101-19751231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19760101-19801231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19810101-19851231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19860101-19901231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19910101-19951231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_19960101-20001231.nc",
                     "tas_day_GFDL-ESM2M_historical_r1i1p1_20010101-20051231.nc")
GFDL_ESM2M_rcp85 <- c("tas_day_GFDL-ESM2M_rcp85_r1i1p1_20060101-20101231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20110101-20151231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20160101-20201231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20210101-20251231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20260101-20301231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20310101-20351231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20360101-20401231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20410101-20451231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20460101-20501231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20510101-20551231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20560101-20601231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20610101-20651231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20660101-20701231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20710101-20751231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20760101-20801231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20810101-20851231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20860101-20901231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20910101-20951231.nc",
                      "tas_day_GFDL-ESM2M_rcp85_r1i1p1_20960101-21001231.nc")
GFDL_ESM2M_08 <- list(GFDL_ESM2M_hist, GFDL_ESM2M_rcp85)

HadGEM2_CC_hist <- c("tas_day_HadGEM2-CC_historical_r1i1p1_18691201-18741130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_18741201-18791130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_18791201-18841130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_18841201-18891130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_18891201-18941130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_18941201-18991130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_18991201-19041130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19041201-19091130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19091201-19141130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19141201-19191130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19191201-19241130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19241201-19291130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19291201-19341130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19341201-19391130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19391201-19441130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19441201-19491130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19491201-19541130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19541201-19591130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19591201-19641130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19641201-19691130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19691201-19741130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19741201-19791130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19791201-19841130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19841201-19891130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19891201-19941130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19941201-19991130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_19991201-20041130.nc",
                     "tas_day_HadGEM2-CC_historical_r1i1p1_20041201-20051130.nc")
HadGEM2_CC_rcp85 <- c("tas_day_HadGEM2-CC_rcp85_r1i1p1_20051201-20101130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20101201-20151130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20151201-20201130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20201201-20251130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20251201-20301130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20301201-20351130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20351201-20401130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20401201-20451130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20451201-20501130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20501201-20551130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20551201-20601130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20601201-20651130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20651201-20701130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20701201-20751130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20751201-20801130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20801201-20851130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20851201-20901130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20901201-20951130.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_20951201-20991230.nc",
                      "tas_day_HadGEM2-CC_rcp85_r1i1p1_21000101-21001230.nc")
HadGEM2_CC_09 <- list(HadGEM2_CC_hist, HadGEM2_CC_rcp85)

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
HadGEM2_ES_10 <- list(HadGEM2_ES_hist, HadGEM2_ES_rcp85)

HadGEM2_AO_hist <- c("tas_day_HadGEM2-AO_historical_r1i1p1_18600101-20051230.nc")
HadGEM2_AO_rcp85 <- c("tas_day_HadGEM2-AO_rcp85_r1i1p1_20060101-21001230.nc")
HadGEM2_AO_11 <- list(HadGEM2_AO_hist, HadGEM2_AO_rcp85)

IPSL_CM5A_LR_hist <- c("tas_day_IPSL-CM5A-LR_historical_r1i1p1_18500101-19491231.nc",
                       "tas_day_IPSL-CM5A-LR_historical_r1i1p1_19500101-20051231.nc")
IPSL_CM5A_LR_rcp85 <- c("tas_day_IPSL-CM5A-LR_rcp85_r1i1p1_20060101-22051231.nc")
IPSL_CM5A_LR_12 <- list(IPSL_CM5A_LR_hist, IPSL_CM5A_LR_rcp85)

IPSL_CM5B_LR_hist <- c("tas_day_IPSL-CM5B-LR_historical_r1i1p1_18500101-20051231.nc")
IPSL_CM5B_LR_rcp85 <- c("tas_day_IPSL-CM5B-LR_rcp85_r1i1p1_20060101-21001231.nc")
IPSL_CM5B_LR_13 <- list(IPSL_CM5B_LR_hist, IPSL_CM5B_LR_rcp85)

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
                 "tas_day_MIROC5_historical_r1i1p1_20000101-20091231.nc",
                 "tas_day_MIROC5_historical_r1i1p1_20100101-20121231.nc")
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
MIROC5_14 <- list(MIROC5_hist, MIROC5_rcp85)

MIROC5_ESM_CHEM_hist <- c('tas_day_MIROC-ESM-CHEM_historical_r1i1p1_18500101-20051231.nc') 
MIROC5_ESM_CHEM_rcp85 <- c("tas_day_MIROC-ESM-CHEM_rcp85_r1i1p1_20060101-21001231.nc")
MIROC5_ESM_CHEM_15 <- list(MIROC5_ESM_CHEM_hist, MIROC5_ESM_CHEM_rcp85)

MIROC5_ESM_hist <- c('tas_day_MIROC-ESM_historical_r3i1p1_18500101-20051231.nc')
MIROC5_ESM_rcp85 <- c("tas_day_MIROC-ESM_rcp85_r1i1p1_20060101-21001231.nc")
MIROC5_ESM_16 <- list(MIROC5_ESM_hist, MIROC5_ESM_rcp85)

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
inmcm4_17 <- list(inmcm4_hist, inmcm4_rcp85)

CNRM_CM5_hist <- c("tas_day_CNRM-CM5_historical_r10i1p1_18700101-18741231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_18750101-18791231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_18800101-18841231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_18850101-18891231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_18900101-18941231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_18950101-18991231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19000101-19041231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19050101-19091231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19100101-19141231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19150101-19191231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19200101-19241231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19250101-19291231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19300101-19341231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19350101-19391231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19400101-19441231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19450101-19491231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19500101-19541231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19550101-19591231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19600101-19641231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19650101-19691231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19700101-19741231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19750101-19791231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19800101-19841231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19850101-19891231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19900101-19941231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_19950101-19991231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_20000101-20041231.nc",
                   "tas_day_CNRM-CM5_historical_r10i1p1_20050101-20051231.nc")
CNRM_CM5_rcp85 <- c("tas_day_CNRM-CM5_rcp85_r1i1p1_20060101-20101231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20110101-20151231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20160101-20201231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20210101-20251231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20260101-20301231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20310101-20351231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20360101-20401231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20410101-20451231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20460101-20501231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20510101-20551231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20560101-20601231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20610101-20651231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20660101-20701231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20710101-20751231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20760101-20801231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20810101-20851231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20860101-20901231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20910101-20951231.nc",
                    "tas_day_CNRM-CM5_rcp85_r1i1p1_20960101-21001231.nc")
CNRM_CM5_18 <- list(CNRM_CM5_hist, CNRM_CM5_rcp85)

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
MRI_CGCM3_19 <- list(MRI_CGCM3_hist, MRI_CGCM3_rcp85)

MRI_ESM1_hist <- c("tas_day_MRI-ESM1_historical_r1i1p1_18610101-18701231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_18710101-18801231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_18810101-18901231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_18910101-19001231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19010101-19101231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19110101-19201231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19210101-19301231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19310101-19401231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19410101-19501231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19510101-19601231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19610101-19701231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19710101-19801231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19810101-19901231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_19910101-20001231.nc",
                   "tas_day_MRI-ESM1_historical_r1i1p1_20010101-20051231.nc") 
MRI_ESM1_rcp85 <- c("tas_day_MRI-ESM1_rcp85_r1i1p1_20060101-20151231.nc",
                    "tas_day_MRI-ESM1_rcp85_r1i1p1_20160101-20251231.nc",
                    "tas_day_MRI-ESM1_rcp85_r1i1p1_20260101-20351231.nc",
                    "tas_day_MRI-ESM1_rcp85_r1i1p1_20360101-20451231.nc",
                    "tas_day_MRI-ESM1_rcp85_r1i1p1_20460101-20551231.nc",
                    "tas_day_MRI-ESM1_rcp85_r1i1p1_20560101-20651231.nc",
                    "tas_day_MRI-ESM1_rcp85_r1i1p1_20660101-20751231.nc",
                    "tas_day_MRI-ESM1_rcp85_r1i1p1_20760101-20851231.nc",
                    "tas_day_MRI-ESM1_rcp85_r1i1p1_20860101-20951231.nc",
                    "tas_day_MRI-ESM1_rcp85_r1i1p1_20960101-21001231.nc") 
MRI_ESM1_20 <- list(MRI_ESM1_hist, MRI_ESM1_rcp85)

IPSL_CM5A_MR_hist <- c("tas_day_IPSL-CM5A-MR_historical_r1i1p1_18500101-18991231.nc",
                       "tas_day_IPSL-CM5A-MR_historical_r1i1p1_19000101-19491231.nc",
                       "tas_day_IPSL-CM5A-MR_historical_r1i1p1_19500101-19991231.nc",
                       "tas_day_IPSL-CM5A-MR_historical_r1i1p1_20000101-20051231.nc")
IPSL_CM5A_MR_rcp85 <- c("tas_day_IPSL-CM5A-MR_rcp85_r1i1p1_20060101-20551231.nc",
                        "tas_day_IPSL-CM5A-MR_rcp85_r1i1p1_20560101-21001231.nc")
IPSL_CM5A_MR_21 <- list(IPSL_CM5A_MR_hist, IPSL_CM5A_MR_rcp85)


gcm_files <- list(CMCC_CESM_01, 
                  CMCC_CM_02, 
                  CMCC_CMS_03,
                  MPI_ESM_LR_04,
                  MPI_ESM_MR_05,
                  GFDL_ESM2G_06,
                  GFDL_CM3_07,
                  GFDL_ESM2M_08,
                  HadGEM2_CC_09,
                  HadGEM2_ES_10,
                  HadGEM2_AO_11, 
                  IPSL_CM5A_LR_12,
                  IPSL_CM5B_LR_13,
                  MIROC5_14,
                  MIROC5_ESM_CHEM_15,
                  MIROC5_ESM_16,
                  inmcm4_17, 
                  CNRM_CM5_18,
                  MRI_CGCM3_19,
                  MRI_ESM1_20,
                  IPSL_CM5A_MR_21)


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

reorganize = reorganize_GCM(historical_filenames = hfn, rcp85_filenames = rfn, path = p)

