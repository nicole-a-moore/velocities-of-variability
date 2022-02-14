library(ncdf4)
library(tidyverse)
library(raster)

## most GCMs of tos were resampled in NCL because they were curvilinear
## however, HadGem GCMs were not, so they were resampled in R

#################################################
###                   FUNCTIONS                ## 
#################################################
resample_GCM <- function(historical_filenames, rcp85_filenames, path) {
  
  #####      EXTRACTING DATA FROM FILES      #####
  ## create vector of file names:
  filenames <- append(historical_filenames, rcp85_filenames)
  paths <- paste(path, filenames, sep = "")
  ## raster object of desired resolution/extent
  r <- raster(ymx = 90, ymn = -90, xmn = 0, xmx = 360, res = 1) 
  
  ## for each file in the set of GCM files:
  file = 1
  dates <- c()
  while (file < (length(filenames)+1)) { 
    
    ## create raster stack to store data 
    og_temps = brick(paths[file], stopIfNotEqualSpaced = FALSE)
    dates <- append(dates, names(og_temps))
    
    #####      STANDARDIZE TO 1X1 DEGREE GRID   #####
    ## resample temperatures so all GCMs are on a 1 degree x 1 degree grid of the same extent
    new_temps <- resample(og_temps, r, method = 'bilinear',
                          filename = paste(path, "resampled_", filenames[file], sep = ""))
    
    print(paste0("Done file number: ", file), stdout())
    file = file + 1
  }
  dates <- str_replace_all(dates, "X","")
  saveRDS(dates, paste(path, "dates.rds", sep = ""))
  return(dates)
}


#################################################
###      setting paths, calling function       ## 
#################################################
path = "data-raw/09_HadGEM2-CC/"
historical_filenames = c("tos_day_HadGEM2-CC_historical_r1i1p1_18691201-18741130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_18741201-18791130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_18791201-18841130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_18841201-18891130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_18891201-18941130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_18941201-18991130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_18991201-19041130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19041201-19091130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19091201-19141130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19141201-19191130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19191201-19241130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19241201-19291130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19291201-19341130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19341201-19391130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19391201-19441130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19441201-19491130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19491201-19541130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19541201-19591130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19591201-19641130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19641201-19691130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19691201-19741130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19741201-19791130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19791201-19841130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19841201-19891130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19891201-19941130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19941201-19991130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_19991201-20041130.nc",
                         "tos_day_HadGEM2-CC_historical_r1i1p1_20041201-20051130.nc")
                       

rcp85_filenames = c("tos_day_HadGEM2-CC_rcp85_r1i1p1_20051201-20101130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20101201-20151130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20151201-20201130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20201201-20251130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20251201-20301130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20301201-20351130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20351201-20401130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20401201-20451130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20451201-20501130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20501201-20551130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20551201-20601130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20601201-20651130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20651201-20701130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20701201-20751130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20751201-20801130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20801201-20851130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20851201-20901130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20901201-20951130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_20951201-20991130.nc",
                    "tos_day_HadGEM2-CC_rcp85_r1i1p1_21000101-21001230.nc")

resample_GCM(historical_filenames, rcp85_filenames, path)

path = "data-raw/10_HadGEM2-ES/"

historical_filenames = c("tos_day_HadGEM2-ES_rcp85_r1i1p1_20991201-21091130.nc",
                         "tos_day_HadGEM2-ES_historical_r1i1p1_18691201-18791130.nc",
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

rcp85_filenames <- c("tos_day_HadGEM2-ES_rcp85_r1i1p1_20051201-20101130.nc",
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
                     "tos_day_HadGEM2-ES_rcp85_r1i1p1_20851201-20991130.nc")

resample_GCM(historical_filenames, rcp85_filenames, path)

