### analyzing era 40 reanalysis of historical observations of air and sea surface temperature 
library(tidyverse)
library(ncdf4)
library(easyNCDF)
library(raster)


#################################################
###            resample observations           ## 
#################################################
r <- raster(ymx = 90, ymn = -90, xmn = 0, xmx = 360, res = 1) 
names_tas <- c()
names_sst <- c()
i=0
while (i < length(filenames)) {
  ## open file, get time variable
  file = paste(path, filenames[length(filenames)-i], sep = "")
  nc <- nc_open(file)
  time <- ncvar_get(nc, "initial_time0") #hours since 1800-01-01 00:00"
  nc_close(nc)
  
  ## rasterize sst and tas:
  sst <- stack(file, varname="SSTK_GDS4_SFC")
  tas <- stack(file, varname="2T_GDS4_SFC")
  
  if (any(str_detect(time, "02/29"))) {
    sst <- sst[[-which(str_detect(time, "02/29"))]] ## get rid of leap year day
    tas <- tas[[-which(str_detect(time, "02/29"))]]
    time <- time[-which(str_detect(time, "02/29"))]
  }

  ## resample to mean daily temp by averaging across samples witin a day (4)
  indices <- rep(1:(length(time)/4),each = 4)
  sst_daily <- stackApply(sst, indices, fun = mean)
  tas_daily <- stackApply(tas, indices, fun = mean)
  
  #####      STANDARDIZE TO 1X1 DEGREE GRID   #####
  ## raster object of desired resolution/extent
  if (i %% 10 != 0) {
    if(i %% 10 == 1){
      rtas <- tas_daily
      rsst <- sst_daily
    }
    else {
      rtas <- stack(rtas, tas_daily)
      rsst <- stack(rsst, sst_daily)
    }
   
  }
  else {
    rtas <- resample(stack(rtas, tas_daily), r, method = 'bilinear', 
                     filename = paste("data-processed/ERA-40/resampled_tas_", 
                                      filenames[length(filenames)-i], sep = ""))
    rsst <- resample(stack(rsst, sst_daily), r, method = 'bilinear',
                     filename = paste("data-processed/ERA-40/resampled_sst_", filenames[length(filenames)-i], sep = ""))
    names_tas <- append(names_tas, paste("data-processed/ERA-40/resampled_tas_", filenames[length(filenames)-i], sep = ""))
    names_sst <- append(names_sst, paste("data-processed/ERA-40/resampled_sst_", filenames[length(filenames)-i], sep = ""))
  }
  
  print(paste0("Done file number: ", i), stdout())
  i = i+2
}

saveRDS(names_sst, 
        "data-processed/ERA-40/sst_filenames.rds")
saveRDS(names_tas, 
        "data-processed/ERA-40/tas_filenames.rds")


## get paths and filenames 
files_sst <- readRDS("data-processed/ERA-40/sst_filenames.rds")
files_tas <- readRDS("data-processed/ERA-40/tas_filenames.rds")

## make date vector
## first date: 09/01/1957
## last date: 08/31/2002
dates <- paste(rep(seq(1957, 2002), each = 365), ".", rep(seq(1:365), 45), sep = "")
dates <- dates[244:(45*365+243)]

reorganize_GCM(filenames = files_sst, type = "sst")
reorganize_GCM(filenames = files_tas, type = "tas")

## detrend and make spatial chunks 
reorganize_GCM <- function(filenames, type) {

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
      temps <- stack(filenames[file])
      
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
      print(paste0("file: ", file))
      file = file + 1
    }
    
    ## turn into an array:
    temps_df <- as.array(spatial_temps)  
    
    ## loop through cells in spatial chunk, removing outliers and interpolating missing vals in time series
    l_detrended_temps <- s_detrended_temps <- temps_df
    x = 1 ## longitude
    while (x < ncol(temps_df) + 1) {
      y = 1 ## latitude
      while (y < nrow(temps_df)+1) {
        ## get local time series 
        local_ts <- data.frame(time = 1:(length(temps_df[1,1,])), 
                               temp = temps_df[y,x,])
        
        # ggplot(local_ts, aes(x = time, y = temp)) + geom_line() +
        #   geom_smooth(method = "lm")
        
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
          local_ts$md <- str_split_fixed(dates, "\\.", n=2)[,2]
          
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
    sp_files[count] <- paste("data-processed/ERA-40/era-40_", type, "_spatial_temps_lon-", lon_bound1,"-", lon_bound2,
                             "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = "")
    # ArrayToNc(temps_df, file_path = paste(path, "spatial_temps_lon-", lon_bound1,"-", lon_bound2,
    #                                       "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(l_detrended_temps, file_path = paste("data-processed/ERA-40/era-40_", type, "_l-detrended_lon-", lon_bound1,"-", lon_bound2,
                                                   "_lat-", lat_bound1, "-", lat_bound2,".nc", sep = ""))
    ArrayToNc(s_detrended_temps, file_path = paste("data-processed/ERA-40/era-40_", type, "_s-detrended_lon-", lon_bound1,"-", lon_bound2,
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
  
  saveRDS(sp_files, paste("data-processed/ERA-40/era-40_", type, "_sp_files.rds", sep = ""))
  ## returns list of spatial chunk filenames 
  return(sp_files)
}


## calculate spectral exponent after getting rid of sea in tas



#################################################
###                 filenames                  ## 
#################################################
path = "data-raw/ERA-40/"

filenames <- c("e4oper.an.sfc.200208.moore543260.nc",
               "e4oper.an.sfc.200207.moore543260.nc",
               "e4oper.an.sfc.200206.moore543260.nc",
               "e4oper.an.sfc.200205.moore543260.nc",
               "e4oper.an.sfc.200204.moore543260.nc",
               "e4oper.an.sfc.200203.moore543260.nc",
               "e4oper.an.sfc.200202.moore543260.nc",
               "e4oper.an.sfc.200201.moore543260.nc",
               "e4oper.an.sfc.200112.moore543260.nc",
               "e4oper.an.sfc.200111.moore543260.nc",
               "e4oper.an.sfc.200110.moore543260.nc",
               "e4oper.an.sfc.200109.moore543260.nc",
               "e4oper.an.sfc.200108.moore543260.nc",
               "e4oper.an.sfc.200107.moore543260.nc",
               "e4oper.an.sfc.200106.moore543260.nc",
               "e4oper.an.sfc.200105.moore543260.nc",
               "e4oper.an.sfc.200104.moore543260.nc",
               "e4oper.an.sfc.200103.moore543260.nc",
               "e4oper.an.sfc.200102.moore543260.nc",
               "e4oper.an.sfc.200101.moore543260.nc",
               "e4oper.an.sfc.200012.moore543260.nc",
               "e4oper.an.sfc.200011.moore543260.nc",
               "e4oper.an.sfc.200010.moore543260.nc",
               "e4oper.an.sfc.200009.moore543260.nc",
               "e4oper.an.sfc.200008.moore543260.nc",
               "e4oper.an.sfc.200007.moore543260.nc",
               "e4oper.an.sfc.200006.moore543260.nc",
               "e4oper.an.sfc.200005.moore543260.nc",
               "e4oper.an.sfc.200004.moore543260.nc",
               "e4oper.an.sfc.200003.moore543260.nc",
               "e4oper.an.sfc.200002.moore543260.nc",
               "e4oper.an.sfc.200001.moore543260.nc",
               "e4oper.an.sfc.199912.moore543260.nc",
               "e4oper.an.sfc.199911.moore543260.nc",
               "e4oper.an.sfc.199910.moore543260.nc",
               "e4oper.an.sfc.199909.moore543260.nc",
               "e4oper.an.sfc.199908.moore543260.nc",
               "e4oper.an.sfc.199907.moore543260.nc",
               "e4oper.an.sfc.199906.moore543260.nc",
               "e4oper.an.sfc.199905.moore543260.nc",
               "e4oper.an.sfc.199904.moore543260.nc",
               "e4oper.an.sfc.199903.moore543260.nc",
               "e4oper.an.sfc.199902.moore543260.nc",
               "e4oper.an.sfc.199901.moore543260.nc",
               "e4oper.an.sfc.199812.moore543260.nc",
               "e4oper.an.sfc.199811.moore543260.nc",
               "e4oper.an.sfc.199810.moore543260.nc",
               "e4oper.an.sfc.199809.moore543260.nc",
               "e4oper.an.sfc.199808.moore543260.nc",
               "e4oper.an.sfc.199807.moore543260.nc",
               "e4oper.an.sfc.199806.moore543260.nc",
               "e4oper.an.sfc.199805.moore543260.nc",
               "e4oper.an.sfc.199804.moore543260.nc",
               "e4oper.an.sfc.199803.moore543260.nc",
               "e4oper.an.sfc.199802.moore543260.nc",
               "e4oper.an.sfc.199801.moore543260.nc",
               "e4oper.an.sfc.199712.moore543260.nc",
               "e4oper.an.sfc.199711.moore543260.nc",
               "e4oper.an.sfc.199710.moore543260.nc",
               "e4oper.an.sfc.199709.moore543260.nc",
               "e4oper.an.sfc.199708.moore543260.nc",
               "e4oper.an.sfc.199707.moore543260.nc",
               "e4oper.an.sfc.199706.moore543260.nc",
               "e4oper.an.sfc.199705.moore543260.nc",
               "e4oper.an.sfc.199704.moore543260.nc",
               "e4oper.an.sfc.199703.moore543260.nc",
               "e4oper.an.sfc.199702.moore543260.nc",
               "e4oper.an.sfc.199701.moore543260.nc",
               "e4oper.an.sfc.199612.moore543260.nc",
               "e4oper.an.sfc.199611.moore543260.nc",
               "e4oper.an.sfc.199610.moore543260.nc",
               "e4oper.an.sfc.199609.moore543260.nc",
               "e4oper.an.sfc.199608.moore543260.nc",
               "e4oper.an.sfc.199607.moore543260.nc",
               "e4oper.an.sfc.199606.moore543260.nc",
               "e4oper.an.sfc.199605.moore543260.nc",
               "e4oper.an.sfc.199604.moore543260.nc",
               "e4oper.an.sfc.199603.moore543260.nc",
               "e4oper.an.sfc.199602.moore543260.nc",
               "e4oper.an.sfc.199601.moore543260.nc",
               "e4oper.an.sfc.199512.moore543260.nc",
               "e4oper.an.sfc.199511.moore543260.nc",
               "e4oper.an.sfc.199510.moore543260.nc",
               "e4oper.an.sfc.199509.moore543260.nc",
               "e4oper.an.sfc.199508.moore543260.nc",
               "e4oper.an.sfc.199507.moore543260.nc",
               "e4oper.an.sfc.199506.moore543260.nc",
               "e4oper.an.sfc.199505.moore543260.nc",
               "e4oper.an.sfc.199504.moore543260.nc",
               "e4oper.an.sfc.199503.moore543260.nc",
               "e4oper.an.sfc.199502.moore543260.nc",
               "e4oper.an.sfc.199501.moore543260.nc",
               "e4oper.an.sfc.199412.moore543260.nc",
               "e4oper.an.sfc.199411.moore543260.nc",
               "e4oper.an.sfc.199410.moore543260.nc",
               "e4oper.an.sfc.199409.moore543260.nc",
               "e4oper.an.sfc.199408.moore543260.nc",
               "e4oper.an.sfc.199407.moore543260.nc",
               "e4oper.an.sfc.199406.moore543260.nc",
               "e4oper.an.sfc.199405.moore543260.nc",
               "e4oper.an.sfc.199404.moore543260.nc",
               "e4oper.an.sfc.199403.moore543260.nc",
               "e4oper.an.sfc.199402.moore543260.nc",
               "e4oper.an.sfc.199401.moore543260.nc",
               "e4oper.an.sfc.199312.moore543260.nc",
               "e4oper.an.sfc.199311.moore543260.nc",
               "e4oper.an.sfc.199310.moore543260.nc",
               "e4oper.an.sfc.199309.moore543260.nc",
               "e4oper.an.sfc.199308.moore543260.nc",
               "e4oper.an.sfc.199307.moore543260.nc",
               "e4oper.an.sfc.199306.moore543260.nc",
               "e4oper.an.sfc.199305.moore543260.nc",
               "e4oper.an.sfc.199304.moore543260.nc",
               "e4oper.an.sfc.199303.moore543260.nc",
               "e4oper.an.sfc.199302.moore543260.nc",
               "e4oper.an.sfc.199301.moore543260.nc",
               "e4oper.an.sfc.199212.moore543260.nc",
               "e4oper.an.sfc.199211.moore543260.nc",
               "e4oper.an.sfc.199210.moore543260.nc",
               "e4oper.an.sfc.199209.moore543260.nc",
               "e4oper.an.sfc.199208.moore543260.nc",
               "e4oper.an.sfc.199207.moore543260.nc",
               "e4oper.an.sfc.199206.moore543260.nc",
               "e4oper.an.sfc.199205.moore543260.nc",
               "e4oper.an.sfc.199204.moore543260.nc",
               "e4oper.an.sfc.199203.moore543260.nc",
               "e4oper.an.sfc.199202.moore543260.nc",
               "e4oper.an.sfc.199201.moore543260.nc",
               "e4oper.an.sfc.199112.moore543260.nc",
               "e4oper.an.sfc.199111.moore543260.nc",
               "e4oper.an.sfc.199110.moore543260.nc",
               "e4oper.an.sfc.199109.moore543260.nc",
               "e4oper.an.sfc.199108.moore543260.nc",
               "e4oper.an.sfc.199107.moore543260.nc",
               "e4oper.an.sfc.199106.moore543260.nc",
               "e4oper.an.sfc.199105.moore543260.nc",
               "e4oper.an.sfc.199104.moore543260.nc",
               "e4oper.an.sfc.199103.moore543260.nc",
               "e4oper.an.sfc.199102.moore543260.nc",
               "e4oper.an.sfc.199101.moore543260.nc",
               "e4oper.an.sfc.199012.moore543260.nc",
               "e4oper.an.sfc.199011.moore543260.nc",
               "e4oper.an.sfc.199010.moore543260.nc",
               "e4oper.an.sfc.199009.moore543260.nc",
               "e4oper.an.sfc.199008.moore543260.nc",
               "e4oper.an.sfc.199007.moore543260.nc",
               "e4oper.an.sfc.199006.moore543260.nc",
               "e4oper.an.sfc.199005.moore543260.nc",
               "e4oper.an.sfc.199004.moore543260.nc",
               "e4oper.an.sfc.199003.moore543260.nc",
               "e4oper.an.sfc.199002.moore543260.nc",
               "e4oper.an.sfc.199001.moore543260.nc",
               "e4oper.an.sfc.198912.moore543260.nc",
               "e4oper.an.sfc.198911.moore543260.nc",
               "e4oper.an.sfc.198910.moore543260.nc",
               "e4oper.an.sfc.198909.moore543260.nc",
               "e4oper.an.sfc.198908.moore543260.nc",
               "e4oper.an.sfc.198907.moore543260.nc",
               "e4oper.an.sfc.198906.moore543260.nc",
               "e4oper.an.sfc.198905.moore543260.nc",
               "e4oper.an.sfc.198904.moore543260.nc",
               "e4oper.an.sfc.198903.moore543260.nc",
               "e4oper.an.sfc.198902.moore543260.nc",
               "e4oper.an.sfc.198901.moore543260.nc",
               "e4oper.an.sfc.198812.moore543260.nc",
               "e4oper.an.sfc.198811.moore543260.nc",
               "e4oper.an.sfc.198810.moore543260.nc",
               "e4oper.an.sfc.198809.moore543260.nc",
               "e4oper.an.sfc.198808.moore543260.nc",
               "e4oper.an.sfc.198807.moore543260.nc",
               "e4oper.an.sfc.198806.moore543260.nc",
               "e4oper.an.sfc.198805.moore543260.nc",
               "e4oper.an.sfc.198804.moore543260.nc",
               "e4oper.an.sfc.198803.moore543260.nc",
               "e4oper.an.sfc.198802.moore543260.nc",
               "e4oper.an.sfc.198801.moore543260.nc",
               "e4oper.an.sfc.198712.moore543260.nc",
               "e4oper.an.sfc.198711.moore543260.nc",
               "e4oper.an.sfc.198710.moore543260.nc",
               "e4oper.an.sfc.198709.moore543260.nc",
               "e4oper.an.sfc.198708.moore543260.nc",
               "e4oper.an.sfc.198707.moore543260.nc",
               "e4oper.an.sfc.198706.moore543260.nc",
               "e4oper.an.sfc.198705.moore543260.nc",
               "e4oper.an.sfc.198704.moore543260.nc",
               "e4oper.an.sfc.198703.moore543260.nc",
               "e4oper.an.sfc.198702.moore543260.nc",
               "e4oper.an.sfc.198701.moore543260.nc",
               "e4oper.an.sfc.198612.moore543260.nc",
               "e4oper.an.sfc.198611.moore543260.nc",
               "e4oper.an.sfc.198610.moore543260.nc",
               "e4oper.an.sfc.198609.moore543260.nc",
               "e4oper.an.sfc.198608.moore543260.nc",
               "e4oper.an.sfc.198607.moore543260.nc",
               "e4oper.an.sfc.198606.moore543260.nc",
               "e4oper.an.sfc.198605.moore543260.nc",
               "e4oper.an.sfc.198604.moore543260.nc",
               "e4oper.an.sfc.198603.moore543260.nc",
               "e4oper.an.sfc.198602.moore543260.nc",
               "e4oper.an.sfc.198601.moore543260.nc",
               "e4oper.an.sfc.198512.moore543260.nc",
               "e4oper.an.sfc.198511.moore543260.nc",
               "e4oper.an.sfc.198510.moore543260.nc",
               "e4oper.an.sfc.198509.moore543260.nc",
               "e4oper.an.sfc.198508.moore543260.nc",
               "e4oper.an.sfc.198507.moore543260.nc",
               "e4oper.an.sfc.198506.moore543260.nc",
               "e4oper.an.sfc.198505.moore543260.nc",
               "e4oper.an.sfc.198504.moore543260.nc",
               "e4oper.an.sfc.198503.moore543260.nc",
               "e4oper.an.sfc.198502.moore543260.nc",
               "e4oper.an.sfc.198501.moore543260.nc",
               "e4oper.an.sfc.198412.moore543260.nc",
               "e4oper.an.sfc.198411.moore543260.nc",
               "e4oper.an.sfc.198410.moore543260.nc",
               "e4oper.an.sfc.198409.moore543260.nc",
               "e4oper.an.sfc.198408.moore543260.nc",
               "e4oper.an.sfc.198407.moore543260.nc",
               "e4oper.an.sfc.198406.moore543260.nc",
               "e4oper.an.sfc.198405.moore543260.nc",
               "e4oper.an.sfc.198404.moore543260.nc",
               "e4oper.an.sfc.198403.moore543260.nc",
               "e4oper.an.sfc.198402.moore543260.nc",
               "e4oper.an.sfc.198401.moore543260.nc",
               "e4oper.an.sfc.198312.moore543260.nc",
               "e4oper.an.sfc.198311.moore543260.nc",
               "e4oper.an.sfc.198310.moore543260.nc",
               "e4oper.an.sfc.198309.moore543260.nc",
               "e4oper.an.sfc.198308.moore543260.nc",
               "e4oper.an.sfc.198307.moore543260.nc",
               "e4oper.an.sfc.198306.moore543260.nc",
               "e4oper.an.sfc.198305.moore543260.nc",
               "e4oper.an.sfc.198304.moore543260.nc",
               "e4oper.an.sfc.198303.moore543260.nc",
               "e4oper.an.sfc.198302.moore543260.nc",
               "e4oper.an.sfc.198301.moore543260.nc",
               "e4oper.an.sfc.198212.moore543260.nc",
               "e4oper.an.sfc.198211.moore543260.nc",
               "e4oper.an.sfc.198210.moore543260.nc",
               "e4oper.an.sfc.198209.moore543260.nc",
               "e4oper.an.sfc.198208.moore543260.nc",
               "e4oper.an.sfc.198207.moore543260.nc",
               "e4oper.an.sfc.198206.moore543260.nc",
               "e4oper.an.sfc.198205.moore543260.nc",
               "e4oper.an.sfc.198204.moore543260.nc",
               "e4oper.an.sfc.198203.moore543260.nc",
               "e4oper.an.sfc.198202.moore543260.nc",
               "e4oper.an.sfc.198201.moore543260.nc",
               "e4oper.an.sfc.198112.moore543260.nc",
               "e4oper.an.sfc.198111.moore543260.nc",
               "e4oper.an.sfc.198110.moore543260.nc",
               "e4oper.an.sfc.198109.moore543260.nc",
               "e4oper.an.sfc.198108.moore543260.nc",
               "e4oper.an.sfc.198107.moore543260.nc",
               "e4oper.an.sfc.198106.moore543260.nc",
               "e4oper.an.sfc.198105.moore543260.nc",
               "e4oper.an.sfc.198104.moore543260.nc",
               "e4oper.an.sfc.198103.moore543260.nc",
               "e4oper.an.sfc.198102.moore543260.nc",
               "e4oper.an.sfc.198101.moore543260.nc",
               "e4oper.an.sfc.198012.moore543260.nc",
               "e4oper.an.sfc.198011.moore543260.nc",
               "e4oper.an.sfc.198010.moore543260.nc",
               "e4oper.an.sfc.198009.moore543260.nc",
               "e4oper.an.sfc.198008.moore543260.nc",
               "e4oper.an.sfc.198007.moore543260.nc",
               "e4oper.an.sfc.198006.moore543260.nc",
               "e4oper.an.sfc.198005.moore543260.nc",
               "e4oper.an.sfc.198004.moore543260.nc",
               "e4oper.an.sfc.198003.moore543260.nc",
               "e4oper.an.sfc.198002.moore543260.nc",
               "e4oper.an.sfc.198001.moore543260.nc",
               "e4oper.an.sfc.197912.moore543260.nc",
               "e4oper.an.sfc.197911.moore543260.nc",
               "e4oper.an.sfc.197910.moore543260.nc",
               "e4oper.an.sfc.197909.moore543260.nc",
               "e4oper.an.sfc.197908.moore543260.nc",
               "e4oper.an.sfc.197907.moore543260.nc",
               "e4oper.an.sfc.197906.moore543260.nc",
               "e4oper.an.sfc.197905.moore543260.nc",
               "e4oper.an.sfc.197904.moore543260.nc",
               "e4oper.an.sfc.197903.moore543260.nc",
               "e4oper.an.sfc.197902.moore543260.nc",
               "e4oper.an.sfc.197901.moore543260.nc",
               "e4oper.an.sfc.197812.moore543260.nc",
               "e4oper.an.sfc.197811.moore543260.nc",
               "e4oper.an.sfc.197810.moore543260.nc",
               "e4oper.an.sfc.197809.moore543260.nc",
               "e4oper.an.sfc.197808.moore543260.nc",
               "e4oper.an.sfc.197807.moore543260.nc",
               "e4oper.an.sfc.197806.moore543260.nc",
               "e4oper.an.sfc.197805.moore543260.nc",
               "e4oper.an.sfc.197804.moore543260.nc",
               "e4oper.an.sfc.197803.moore543260.nc",
               "e4oper.an.sfc.197802.moore543260.nc",
               "e4oper.an.sfc.197801.moore543260.nc",
               "e4oper.an.sfc.197712.moore543260.nc",
               "e4oper.an.sfc.197711.moore543260.nc",
               "e4oper.an.sfc.197710.moore543260.nc",
               "e4oper.an.sfc.197709.moore543260.nc",
               "e4oper.an.sfc.197708.moore543260.nc",
               "e4oper.an.sfc.197707.moore543260.nc",
               "e4oper.an.sfc.197706.moore543260.nc",
               "e4oper.an.sfc.197705.moore543260.nc",
               "e4oper.an.sfc.197704.moore543260.nc",
               "e4oper.an.sfc.197703.moore543260.nc",
               "e4oper.an.sfc.197702.moore543260.nc",
               "e4oper.an.sfc.197701.moore543260.nc",
               "e4oper.an.sfc.197612.moore543260.nc",
               "e4oper.an.sfc.197611.moore543260.nc",
               "e4oper.an.sfc.197610.moore543260.nc",
               "e4oper.an.sfc.197609.moore543260.nc",
               "e4oper.an.sfc.197608.moore543260.nc",
               "e4oper.an.sfc.197607.moore543260.nc",
               "e4oper.an.sfc.197606.moore543260.nc",
               "e4oper.an.sfc.197605.moore543260.nc",
               "e4oper.an.sfc.197604.moore543260.nc",
               "e4oper.an.sfc.197603.moore543260.nc",
               "e4oper.an.sfc.197602.moore543260.nc",
               "e4oper.an.sfc.197601.moore543260.nc",
               "e4oper.an.sfc.197512.moore543260.nc",
               "e4oper.an.sfc.197511.moore543260.nc",
               "e4oper.an.sfc.197510.moore543260.nc",
               "e4oper.an.sfc.197509.moore543260.nc",
               "e4oper.an.sfc.197508.moore543260.nc",
               "e4oper.an.sfc.197507.moore543260.nc",
               "e4oper.an.sfc.197506.moore543260.nc",
               "e4oper.an.sfc.197505.moore543260.nc",
               "e4oper.an.sfc.197504.moore543260.nc",
               "e4oper.an.sfc.197503.moore543260.nc",
               "e4oper.an.sfc.197502.moore543260.nc",
               "e4oper.an.sfc.197501.moore543260.nc",
               "e4oper.an.sfc.197412.moore543260.nc",
               "e4oper.an.sfc.197411.moore543260.nc",
               "e4oper.an.sfc.197410.moore543260.nc",
               "e4oper.an.sfc.197409.moore543260.nc",
               "e4oper.an.sfc.197408.moore543260.nc",
               "e4oper.an.sfc.197407.moore543260.nc",
               "e4oper.an.sfc.197406.moore543260.nc",
               "e4oper.an.sfc.197405.moore543260.nc",
               "e4oper.an.sfc.197404.moore543260.nc",
               "e4oper.an.sfc.197403.moore543260.nc",
               "e4oper.an.sfc.197402.moore543260.nc",
               "e4oper.an.sfc.197401.moore543260.nc",
               "e4oper.an.sfc.197312.moore543260.nc",
               "e4oper.an.sfc.197311.moore543260.nc",
               "e4oper.an.sfc.197310.moore543260.nc",
               "e4oper.an.sfc.197309.moore543260.nc",
               "e4oper.an.sfc.197308.moore543260.nc",
               "e4oper.an.sfc.197307.moore543260.nc",
               "e4oper.an.sfc.197306.moore543260.nc",
               "e4oper.an.sfc.197305.moore543260.nc",
               "e4oper.an.sfc.197304.moore543260.nc",
               "e4oper.an.sfc.197303.moore543260.nc",
               "e4oper.an.sfc.197302.moore543260.nc",
               "e4oper.an.sfc.197301.moore543260.nc",
               "e4oper.an.sfc.197212.moore543260.nc",
               "e4oper.an.sfc.197211.moore543260.nc",
               "e4oper.an.sfc.197210.moore543260.nc",
               "e4oper.an.sfc.197209.moore543260.nc",
               "e4oper.an.sfc.197208.moore543260.nc",
               "e4oper.an.sfc.197207.moore543260.nc",
               "e4oper.an.sfc.197206.moore543260.nc",
               "e4oper.an.sfc.197205.moore543260.nc",
               "e4oper.an.sfc.197204.moore543260.nc",
               "e4oper.an.sfc.197203.moore543260.nc",
               "e4oper.an.sfc.197202.moore543260.nc",
               "e4oper.an.sfc.197201.moore543260.nc",
               "e4oper.an.sfc.197112.moore543260.nc",
               "e4oper.an.sfc.197111.moore543260.nc",
               "e4oper.an.sfc.197110.moore543260.nc",
               "e4oper.an.sfc.197109.moore543260.nc",
               "e4oper.an.sfc.197108.moore543260.nc",
               "e4oper.an.sfc.197107.moore543260.nc",
               "e4oper.an.sfc.197106.moore543260.nc",
               "e4oper.an.sfc.197105.moore543260.nc",
               "e4oper.an.sfc.197104.moore543260.nc",
               "e4oper.an.sfc.197103.moore543260.nc",
               "e4oper.an.sfc.197102.moore543260.nc",
               "e4oper.an.sfc.197101.moore543260.nc",
               "e4oper.an.sfc.197012.moore543260.nc",
               "e4oper.an.sfc.197011.moore543260.nc",
               "e4oper.an.sfc.197010.moore543260.nc",
               "e4oper.an.sfc.197009.moore543260.nc",
               "e4oper.an.sfc.197008.moore543260.nc",
               "e4oper.an.sfc.197007.moore543260.nc",
               "e4oper.an.sfc.197006.moore543260.nc",
               "e4oper.an.sfc.197005.moore543260.nc",
               "e4oper.an.sfc.197004.moore543260.nc",
               "e4oper.an.sfc.197003.moore543260.nc",
               "e4oper.an.sfc.197002.moore543260.nc",
               "e4oper.an.sfc.197001.moore543260.nc",
               "e4oper.an.sfc.196912.moore543260.nc",
               "e4oper.an.sfc.196911.moore543260.nc",
               "e4oper.an.sfc.196910.moore543260.nc",
               "e4oper.an.sfc.196909.moore543260.nc",
               "e4oper.an.sfc.196908.moore543260.nc",
               "e4oper.an.sfc.196907.moore543260.nc",
               "e4oper.an.sfc.196906.moore543260.nc",
               "e4oper.an.sfc.196905.moore543260.nc",
               "e4oper.an.sfc.196904.moore543260.nc",
               "e4oper.an.sfc.196903.moore543260.nc",
               "e4oper.an.sfc.196902.moore543260.nc",
               "e4oper.an.sfc.196901.moore543260.nc",
               "e4oper.an.sfc.196812.moore543260.nc",
               "e4oper.an.sfc.196811.moore543260.nc",
               "e4oper.an.sfc.196810.moore543260.nc",
               "e4oper.an.sfc.196809.moore543260.nc",
               "e4oper.an.sfc.196808.moore543260.nc",
               "e4oper.an.sfc.196807.moore543260.nc",
               "e4oper.an.sfc.196806.moore543260.nc",
               "e4oper.an.sfc.196805.moore543260.nc",
               "e4oper.an.sfc.196804.moore543260.nc",
               "e4oper.an.sfc.196803.moore543260.nc",
               "e4oper.an.sfc.196802.moore543260.nc",
               "e4oper.an.sfc.196801.moore543260.nc",
               "e4oper.an.sfc.196712.moore543260.nc",
               "e4oper.an.sfc.196711.moore543260.nc",
               "e4oper.an.sfc.196710.moore543260.nc",
               "e4oper.an.sfc.196709.moore543260.nc",
               "e4oper.an.sfc.196708.moore543260.nc",
               "e4oper.an.sfc.196707.moore543260.nc",
               "e4oper.an.sfc.196706.moore543260.nc",
               "e4oper.an.sfc.196705.moore543260.nc",
               "e4oper.an.sfc.196704.moore543260.nc",
               "e4oper.an.sfc.196703.moore543260.nc",
               "e4oper.an.sfc.196702.moore543260.nc",
               "e4oper.an.sfc.196701.moore543260.nc",
               "e4oper.an.sfc.196612.moore543260.nc",
               "e4oper.an.sfc.196611.moore543260.nc",
               "e4oper.an.sfc.196610.moore543260.nc",
               "e4oper.an.sfc.196609.moore543260.nc",
               "e4oper.an.sfc.196608.moore543260.nc",
               "e4oper.an.sfc.196607.moore543260.nc",
               "e4oper.an.sfc.196606.moore543260.nc",
               "e4oper.an.sfc.196605.moore543260.nc",
               "e4oper.an.sfc.196604.moore543260.nc",
               "e4oper.an.sfc.196603.moore543260.nc",
               "e4oper.an.sfc.196602.moore543260.nc",
               "e4oper.an.sfc.196601.moore543260.nc",
               "e4oper.an.sfc.196512.moore543260.nc",
               "e4oper.an.sfc.196511.moore543260.nc",
               "e4oper.an.sfc.196510.moore543260.nc",
               "e4oper.an.sfc.196509.moore543260.nc",
               "e4oper.an.sfc.196508.moore543260.nc",
               "e4oper.an.sfc.196507.moore543260.nc",
               "e4oper.an.sfc.196506.moore543260.nc",
               "e4oper.an.sfc.196505.moore543260.nc",
               "e4oper.an.sfc.196504.moore543260.nc",
               "e4oper.an.sfc.196503.moore543260.nc",
               "e4oper.an.sfc.196502.moore543260.nc",
               "e4oper.an.sfc.196501.moore543260.nc",
               "e4oper.an.sfc.196412.moore543260.nc",
               "e4oper.an.sfc.196411.moore543260.nc",
               "e4oper.an.sfc.196410.moore543260.nc",
               "e4oper.an.sfc.196409.moore543260.nc",
               "e4oper.an.sfc.196408.moore543260.nc",
               "e4oper.an.sfc.196407.moore543260.nc",
               "e4oper.an.sfc.196406.moore543260.nc",
               "e4oper.an.sfc.196405.moore543260.nc",
               "e4oper.an.sfc.196404.moore543260.nc",
               "e4oper.an.sfc.196403.moore543260.nc",
               "e4oper.an.sfc.196402.moore543260.nc",
               "e4oper.an.sfc.196401.moore543260.nc",
               "e4oper.an.sfc.196312.moore543260.nc",
               "e4oper.an.sfc.196311.moore543260.nc",
               "e4oper.an.sfc.196310.moore543260.nc",
               "e4oper.an.sfc.196309.moore543260.nc",
               "e4oper.an.sfc.196308.moore543260.nc",
               "e4oper.an.sfc.196307.moore543260.nc",
               "e4oper.an.sfc.196306.moore543260.nc",
               "e4oper.an.sfc.196305.moore543260.nc",
               "e4oper.an.sfc.196304.moore543260.nc",
               "e4oper.an.sfc.196303.moore543260.nc",
               "e4oper.an.sfc.196302.moore543260.nc",
               "e4oper.an.sfc.196301.moore543260.nc",
               "e4oper.an.sfc.196212.moore543260.nc",
               "e4oper.an.sfc.196211.moore543260.nc",
               "e4oper.an.sfc.196210.moore543260.nc",
               "e4oper.an.sfc.196209.moore543260.nc",
               "e4oper.an.sfc.196208.moore543260.nc",
               "e4oper.an.sfc.196207.moore543260.nc",
               "e4oper.an.sfc.196206.moore543260.nc",
               "e4oper.an.sfc.196205.moore543260.nc",
               "e4oper.an.sfc.196204.moore543260.nc",
               "e4oper.an.sfc.196203.moore543260.nc",
               "e4oper.an.sfc.196202.moore543260.nc",
               "e4oper.an.sfc.196201.moore543260.nc",
               "e4oper.an.sfc.196112.moore543260.nc",
               "e4oper.an.sfc.196111.moore543260.nc",
               "e4oper.an.sfc.196110.moore543260.nc",
               "e4oper.an.sfc.196109.moore543260.nc",
               "e4oper.an.sfc.196108.moore543260.nc",
               "e4oper.an.sfc.196107.moore543260.nc",
               "e4oper.an.sfc.196106.moore543260.nc",
               "e4oper.an.sfc.196105.moore543260.nc",
               "e4oper.an.sfc.196104.moore543260.nc",
               "e4oper.an.sfc.196103.moore543260.nc",
               "e4oper.an.sfc.196102.moore543260.nc",
               "e4oper.an.sfc.196101.moore543260.nc",
               "e4oper.an.sfc.196012.moore543260.nc",
               "e4oper.an.sfc.196011.moore543260.nc",
               "e4oper.an.sfc.196010.moore543260.nc",
               "e4oper.an.sfc.196009.moore543260.nc",
               "e4oper.an.sfc.196008.moore543260.nc",
               "e4oper.an.sfc.196007.moore543260.nc",
               "e4oper.an.sfc.196006.moore543260.nc",
               "e4oper.an.sfc.196005.moore543260.nc",
               "e4oper.an.sfc.196004.moore543260.nc",
               "e4oper.an.sfc.196003.moore543260.nc",
               "e4oper.an.sfc.196002.moore543260.nc",
               "e4oper.an.sfc.196001.moore543260.nc",
               "e4oper.an.sfc.195912.moore543260.nc",
               "e4oper.an.sfc.195911.moore543260.nc",
               "e4oper.an.sfc.195910.moore543260.nc",
               "e4oper.an.sfc.195909.moore543260.nc",
               "e4oper.an.sfc.195908.moore543260.nc",
               "e4oper.an.sfc.195907.moore543260.nc",
               "e4oper.an.sfc.195906.moore543260.nc",
               "e4oper.an.sfc.195905.moore543260.nc",
               "e4oper.an.sfc.195904.moore543260.nc",
               "e4oper.an.sfc.195903.moore543260.nc",
               "e4oper.an.sfc.195902.moore543260.nc",
               "e4oper.an.sfc.195901.moore543260.nc",
               "e4oper.an.sfc.195812.moore543260.nc",
               "e4oper.an.sfc.195811.moore543260.nc",
               "e4oper.an.sfc.195810.moore543260.nc",
               "e4oper.an.sfc.195809.moore543260.nc",
               "e4oper.an.sfc.195808.moore543260.nc",
               "e4oper.an.sfc.195807.moore543260.nc",
               "e4oper.an.sfc.195806.moore543260.nc",
               "e4oper.an.sfc.195805.moore543260.nc",
               "e4oper.an.sfc.195804.moore543260.nc",
               "e4oper.an.sfc.195803.moore543260.nc",
               "e4oper.an.sfc.195802.moore543260.nc",
               "e4oper.an.sfc.195801.moore543260.nc",
               "e4oper.an.sfc.195712.moore543260.nc",
               "e4oper.an.sfc.195711.moore543260.nc",
               "e4oper.an.sfc.195710.moore543260.nc",
               "e4oper.an.sfc.195709.moore543260.nc")
