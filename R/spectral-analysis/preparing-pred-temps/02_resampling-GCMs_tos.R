library(ncdf4)
library(tidyverse)
library(raster)

## most GCMs of tos were resampled in NCL because they were curvilinear
## however, HadGem and MIROC-ESM-CHEM were not, so they were resampled in R

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
                          filename = paste(path, "resampled_", filenames[file], sep = ""),
                          overwrite = TRUE)
    
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
path = "/Volumes/NIKKI/CMIP5-GCMs/"

## create vector of file folders to put data into
gcm_models <- c("04_HadGEM2-ES_tos",
                "07_MIROC-ESM-CHEM_tos")

folders <- paste(path, gcm_models, "/", sep = "")

#################################################
###       make list of GCM file names          ## 
#################################################
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
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20851201-20991130.nc",
                      "tos_day_HadGEM2-ES_rcp85_r1i1p1_20991201-21091130.nc")
HadGEM2_ES_04 <- list(HadGEM2_ES_hist, HadGEM2_ES_rcp85)

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

gcm_files <- list(HadGEM2_ES_04, MIROC_ESM_CHEM_07)

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

d = resample_GCM(historical_filenames = hfn, rcp85_filenames = rfn, path = p)
