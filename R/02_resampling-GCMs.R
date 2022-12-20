### resample organize GCMs to 1x1 degree grid
### R version 3.6.1 
.libPaths(c("~/projects/def-jsunday/nikkim/VoV/packages", .libPaths()))
library(tidyverse)
library(raster)
library(stringr)

#################################################
###                   FUNCTIONS                ## 
#################################################
resample_GCM <- function(historical_filenames, rcp85_filenames, path, file) {
  
  #####      EXTRACTING DATA FROM FILES      #####
  ## create vector of file names:
  filenames <- append(historical_filenames, rcp85_filenames)
  paths <- paste(path, filenames, sep = "")
  ## raster object of desired resolution/extent
  r <- raster(ymx = 90, ymn = -90, xmn = 0, xmx = 360, res = 1) 
  
  ## for each file in the set of GCM files:
  file = file
  dates <- c()
  start <- Sys.time()
  model_3 <- str_detect(path, "03")
  while (file < (length(filenames)+1)) { 
    
    ## create raster stack to store data 
    if(model_3) {
      ## these two annoying models have an unevenly spaced grid. alas, we must get around this...
      og_temps = brick(paths[file], stopIfNotEqualSpaced = FALSE)
    }
    else {
      og_temps <- stack(paths[file])
    }
    
    ## if GCM has overlap between historical and rcp data, crop data at 2005-21-31 / 2006-01-01
    if(historical_filenames[file] == "tas_day_MIROC5_historical_r1i1p1_20000101-20091231.nc") {
      lastday <- which(str_replace_all(names(og_temps), "X","") == "2005.12.31")
      og_temps <- og_temps[[1:lastday]]
      
      ## save dates
      dates <- append(dates, names(og_temps)[1:lastday])
    }
    else {
      ## save dates
      dates <- append(dates, names(og_temps))
    }
    
    #####      STANDARDIZE TO 1X1 DEGREE GRID   #####
    ## resample temperatures so all GCMs are on a 1 degree x 1 degree grid of the same extent
    new_temps <- resample(og_temps, r, method = 'bilinear',
                            filename = paste(path, "test_resampled_", filenames[file], sep = ""),
                          overwrite = TRUE)
    
    print(paste0("Done file number: ", file), stdout())
    file = file + 1
  }
  end <- Sys.time()
  elapsed <- end - start
  print(paste0("This GCM took ", round(elapsed,3), " minutes and", round(elapsed/60,3), " hours to run."), stdout())
  dates <- str_replace_all(dates, "X","")
  saveRDS(dates, paste(path, "dates.rds", sep = ""))
  return(dates)
}

### how to read back in the data from NC:
# nc = nc_open(paste(path, "resampled_", filenames[file], sep = ""))
# temps = ncvar_get(nc, "variable")
# nc_close(nc)
# 
# r = raster(nrow = 180, ncol = 360, res = 1,
#            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
# extent(r) <- c(0,360,-90,90)
# values(r) <- t(temps[1:360,1:180,1])
# plot(r)

#################################################
###                setting paths               ## 
#################################################
## set 'path' to where you have the GCM files stored on your computer
## for me, they are here:
#path = "/Volumes/NIKKI/CMIP5-GCMs/" ## change me
path = "CMIP5-GCMs/"

## create vector of file folders to put data into
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


#################################################
###     calling functions on a GCM             ## 
#################################################
## read the command line arguments and call the functions on the corresponding GCM 

command_args <- commandArgs(trailingOnly = TRUE)
i = as.numeric(command_args[1])
file = as.numeric(command_args[2])

element = gcm_files[[i]]
hfn <- element[[1]]
rfn <- element[[2]]
p <- folders[i]
  
d = resample_GCM(historical_filenames = hfn, rcp85_filenames = rfn, path = p, file = file)
  



