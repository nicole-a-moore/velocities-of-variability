## checking out sample GCM - "CMCC_CESM"
library(tidyverse)
library(ncdf4)
library(abind)
library(sp)
library(sf)
library(raster)
library(broom)
library(evobiR)
library(gridExtra)
library(lubridate)
select <- dplyr::select

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
############
gcm = 2
element = gcm_files[[gcm]]
hfn <- element[[1]]
rfn <- element[[2]]
p <- folders[gcm]

filenames <- append(hfn, rfn)
paths <- paste(p, filenames, sep = "")


i = 1
while (i < length(cmcc_cesm_historical)+length(cmcc_cesm_rcp85)+1) {
  
  ncfile <- nc_open(paths[i])
  
  ## for first file, extract latitude, longitude, air temp and time 
  if (i == 1) {
    lat <- ncvar_get(ncfile, "lat_bnds") 
    long <- ncvar_get(ncfile, "lon_bnds") 
    tas <- ncvar_get(ncfile, "tas") ## air surface temperature in an array of [long, lat, time]
    time <- ncvar_get(ncfile, "time_bnds")
  }
  ## for other files, append air temp and time onto existing variables 
  else {
    tas <- abind(tas, ncvar_get(ncfile, "tas"), along = 3) 
    time_modified <- ncvar_get(ncfile, "time_bnds")
    ## modify time variable to link time series together 
    time_modified[1,] <-  time_modified[1,] + time[1,ncol(time)]+1
    time_modified[2,] <-  time_modified[2,] + time[1,ncol(time)]+1
    time <- cbind(time, time_modified)
  }
  
  ## close the file
  nc_close(ncfile)
  
  ## move to next file
  i = i + 1
}

dimnames(tas) <- NULL

## check that time series have 84,371 days
length(tas[1,1,]) ## it does

## make lat and long coords represent centre of grid cells:
lat <- (lat[1,] + lat[2,]) / 2
long <- (long[1,] + long[2,]) / 2 

## reorganize data so North America is left of Europe when plotted:
long[which(long >= 180)] = long[which(long >= 180)] - 360

## try plotting air temps at one time point:
tp <-  expand.grid(long, lat)
colnames(tp) <- c("longitude", "latitude") 
tp$temp <- as.vector(tas[,,1])

## plot in base R:
r = raster(nrows = 48, ncols = 96, xmn=-180, xmx=176.25, ymn=-87.65043, ymx= 87.65043, 
           crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

raster <- rasterize(tp[,1:2], r, tp[,3])
plot(raster, asp = 1)

## plot in base ggplot:
gg1 <- tp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = temp)) + 
  ggtitle("Mean air surface temperature on Jan 1, 1850") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Temperature (K)")

# plot in C
gg2 <- tp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = temp - 273.15)) + 
  ggtitle("Mean air surface temperature on Jan 1, 1850") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Temperature (C)")

## see how variable temperature is in each location
tas_sd <- tas[,,1]
i=1
while (i < nrow(tas) +1 ) {
  z=1
  while (z < ncol(tas) +1) {
    tas_sd[i,z] <- sd(tas[i,z,])
    z = z +1 
  }
  i = i +1 
}

tp <-  expand.grid(long, lat)
colnames(tp) <- c("longitude", "latitude") 
tp$temp <- as.vector(tas_sd[,])

## plot in base R:
r = raster(nrows = 48, ncols = 96, xmn=-180, xmx=176.25, ymn=-87.65043, ymx= 87.65043, 
           crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

raster <- rasterize(tp[,1:2], r, tp[,3])
plot(raster, asp = 1)

## plot in base ggplot:
gg1 <- tp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = temp)) + 
  ggtitle("SD air surface temperature") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Temperature (C)")

## for now, don't worry about standardizing grid size (will have to standardize across GCMs, though)
## just try FFT 

## create function for making dates since 1850-01-01:
as.date <- function(x, origin = getOption("date.origin")){
  origin <- ifelse(is.null(origin), "1970-01-01", origin)
  as.Date(x, origin)
}
options(date.origin = "1849-12-31")

## 1. Linearly detrend temp time series for each spatial cell separately:
l_detrended_tas <- tas
x = 1 ## represents longitude index
while (x < length(long)+1) {
  y = 1 ## represents latitude index
  while (y < length(lat)+1) {
    local_ts <- l_detrended_tas[x, y, ] ## get the local time series 
    ts_df <- data.frame(time = 1:length(local_ts), ## add simple time column representing days from 1850-01-01
                        temp = local_ts, 
                        date = as.date(1:length(local_ts))) ## add a date column  
    
    ## run linear regression for grid cell 
    output <- lm(ts_df, formula = temp ~ time)
    
    ## extract residuals and add to l_detrended_tas:
    l_detrended_tas[x, y, ] <- output$residuals
    
    y = y + 1
  }
  x = x + 1
}


## save detrended time series object:
##saveRDS(l_detrended_tas, "/Volumes/ADATA HV620/GCM-test/l-detrended-tas.rds")
l_detrended_tas <- readRDS("/Volumes/ADATA HV620/GCM-test/l-detrended-tas.rds")

## try plotting one time point of the new set of detrended temp time series:
detrended_tp <-  expand.grid(long, lat)
colnames(detrended_tp) <- c("longitude", "latitude") 
detrended_tp$temp <- as.vector(l_detrended_tas[,,1])

## plot in base R:
raster <- rasterize(detrended_tp[,1:2], r, detrended_tp[,3])
plot(raster, asp = 1)

## plot in ggplot:
gg2 <- detrended_tp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = temp)) + 
  ggtitle("Linearly detrended air surface temperature (Jan 1, 1850)") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Temperature residual (K)")


## now try with seasonally detrended data
## 2. Seasonally detrend temp time series for each spatial cell separately:
s_detrended_tas <- tas
x = 1 ## represents longitude index
while (x < length(long)+1) {
  y = 1 ## represents latitude index
  while (y < length(lat)+1) {
    local_ts <- s_detrended_tas[x, y, ] ## get the local time series 
    
    ts_df <- data.frame(time = 1:length(local_ts), ## add simple time column representing days from 1850-01-01
                        temp = local_ts, 
                        date = as.date(1:length(local_ts))) ## add a date column
    
    ## compute the within-year temperature profile at each geographical location and subtract from temp
    ts_df <- ts_df %>%
      mutate(year = str_split_fixed(.$date, pattern = "-", n = 2)[,1]) %>% ## add a year column
      group_by(year) %>% ## group by year
      do(mutate(., julian_date = seq(1:length(.$year)))) %>% ## add julian date column
      ungroup() %>%
      group_by(julian_date) %>%
      do(mutate(., temp_profile = mean(.$temp))) %>% ## compute mean temp for each day of year across all years
      ungroup() %>%
      mutate(s_detrended_temp = temp - temp_profile)
    
    ## run linear regression for grid cell 
    output <- lm(ts_df, formula = s_detrended_temp ~ time)
    
    ## extract residuals and add to detrended_tas:
    s_detrended_tas[x, y, ] <- output$residuals
    
    print(paste("On x = ",x, ", y = ", y, sep = ""))
    y = y + 1
  }
  x = x + 1
}


## save detrended time series object:
##saveRDS(s_detrended_tas, "/Volumes/ADATA HV620/GCM-test/s-detrended-tas.rds")
s_detrended_tas <- readRDS("/Volumes/ADATA HV620/GCM-test/s-detrended-tas.rds")


## 3. Compute periodogram using FFT across all years (to start) at each location 
spec2D <- array(dim = c(48, 96)) ## 2D array to store spec exp in 
x = 1
while (x < length(long)+1) {
  y = 1 
  while (y < length(lat)+1) {
    local_ts <- l_detrended_tas[x, y, ] ## get the local detrended time series 
    
    ## calculate periodogram using FFT
    L <- length(local_ts)
    
    # Fourier transform the series: 
    dft <- fft(local_ts)/L
    amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
    amp <- amp[1:(L/2)]	## remove second half of amplitudes (negative half)
    freq <- 1:(L/2)/L ## frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
    
    ## create periodogram data by squaring amplitude of FFT output
    spectral <- data.frame(freq = freq, power = amp^2)
    
    # ## plot spectrum:
    # spectral %>% 
    #   ggplot(aes(x = freq, y = power)) + geom_line() +
    #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")
    
    ## get estimate of spectral exponent over whole time series:
    model_output <- lm(spectral, formula = log10(power) ~ log10(freq)) %>%
      tidy(.) %>%
      filter(term == "log10(freq)") %>%
      mutate(lat = lat[y]) %>%
      mutate(long = long[x]) %>%
      mutate(time_interval = "whole_ts")
    
    ## store average spectral exponent in 2D array represting a static picture of projected variability:
    spec2D[y, x] <- model_output$estimate
    
    ## store:
    if(x == 1 & y == 1) {
      spec_exp <- model_output
    }
    else {
      spec_exp <- rbind(spec_exp, model_output)
    }
    
    y = y + 1
  }
  x = x + 1
}

dimnames(spec2D) <- NULL

## try plotting spectral exponents across space:
sp <-  expand.grid(lat, long)
colnames(sp) <- c("latitude", "longitude") 
sp$spec <- as.vector(spec2D)
sp <- sp %>%
  select(longitude, latitude, spec)

## plot in base R:
raster <- rasterize(sp[,1:2], r, sp[,3])
plot(raster, asp = 1)

## plot in ggplot:
gg3 <- sp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = spec)) + 
  ggtitle("Average spectral exponent over years 1850-2100") + coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "")


## Notes:
##  - dates each file spans is in title of file
##  - bnds means bounds, rows of lat and lon describe extent of each grid square 
##  - must check for gaps in time series and see how they deal with leap years!!!

## idea to find gaps in the time series:
## which(!(time[,2] - time[,1]) == 1)
## i think this will work?


ggsave(path = "figures/", filename = "air-temps-raw_Jan-01-2000.png", gg1, height = 6, width = 9)
ggsave(path = "figures/", filename = "air-temps-detrended_Jan-01-2000.png", gg2, height = 6, width = 9)
ggsave(path = "figures/", filename = "static-spec-exp.png", gg3, height = 6, width = 9)


## function to calculate spectral exponent over time series window
spectral_exponent_calculator <- function(ts_window) {
  l <- length(ts_window)
  
  # Fourier transform the time series window: 
  dft <- fft(ts_window)/l
  amp <- sqrt(Im(dft[-1])^2 + Re(dft[-1])^2) ## get rid of first term (represents DC component - y axis shift)
  amp <- amp[1:(l/2)]	## remove second half of amplitudes (negative half)
  freq <- 1:(l/2)/l ## frequency = period(1 day, 2 days, 3 days.... L/2 days) / length of time series 
  
  ## create periodogram data by squaring amplitude of FFT output
  spectral <- data.frame(freq = freq, power = amp^2)
  
  # ## plot spectrum:
  # spectral %>% 
  #   ggplot(aes(x = freq, y = power)) + geom_line() +
  #   scale_y_log10() + scale_x_log10() + geom_smooth(method = "lm")
  
  ## get estimate of spectral exponent over time series window:
  model_output <- lm(spectral, formula = log10(power) ~ log10(freq)) %>%
    tidy(.) %>%
    filter(term == "log10(freq)")
  
  return(model_output$estimate)
}

## function to calculate spectral exponent over sliding time series windows of varying widths (5-10 years) with a time step of one year 
## returns a list containing: matrix of change in spectral exponents using each window wdith, data frame of spectral exponents within each sliding window across all locations and for all window widths 
sliding_window_spec_exp <- function(detrended_ts) {
  spec3D <- array(dim = c(48, 96, 6)) ## 3D array to store spec exp in where third dimension represents years spectral exponent was calculated over (5,6,7,8,9,10)
  x = 1
  while (x < length(long)+1) {
    y = 1 
    while (y < length(lat)+1) {
      local_ts <- detrended_tas[x, y, ] ## get the local detrended time series 
      
      #########################################
      ##        SENSITIVITY ANALYSIS:        ##
      #########################################
      ## calculate spectral exponent using FFT over n year windows
      ## store spectral exponents and calculate slope
      n = 5
      while (n < 11) {
        ## find width of window in days (365.25*n)
        window_width <- round(365*n, digits = 0)
        
        ## calcualte spectral exponent in each window of n years, moving one window width at a time
        exp <- SlidingWindow(FUN = spectral_exponent_calculator, data = local_ts, step = 365, 
                             window = window_width)
        
        spec_exp <- exp %>%
          as.data.frame() %>%
          mutate(time_window_start = seq(from = 1, by = 365, 
                                         to = (floor(length(local_ts)-window_width)))) %>%
          rename(spec_exp = ".") %>%
          mutate(lat = lat[y]) %>%
          mutate(long = long[x]) %>%
          mutate(time_window_width = paste(n, "years")) %>%
          mutate(time_step = "1y")
        
        # ## plot:
        # spec_exp %>%
        #   ggplot(., aes (x = time_window_start, y = spec_exp)) + geom_point() +
        #   geom_smooth(method = "lm") + labs(x = "Time window start index",
        #                                     y = "Spectral exponent over window")
        
        ## regress spectral exponent and extract slope representing change in spectral exponent over time
        model_output <- lm(spec_exp, formula = spec_exp ~ time_window_start) %>%
          tidy(.) %>%
          filter(term == "time_window_start")
        
        spec_exp <- cbind(spec_exp, model_output)
        
        ## store in database 
        if(x == 1 & y == 1 & n == 5) {
          all_spec_exp <- spec_exp 
        }
        else {
          all_spec_exp <- rbind(all_spec_exp, spec_exp)
        }
        
        ## store in array
        spec3D[y, x, n-4] <- model_output$estimate
        
        ## advance to next time window width
        print(paste("On x = ",x, ", y = ", y, ", window width = ", n, " years.", sep = ""))
        n = n + 1
      }
      
      
      ## advance to next longitude
      y = y + 1
    }
    ## advance to next latitude
    x = x + 1
  }
  
  dimnames(spec3D) <- NULL
  
  return(list(spec3D, all_spec_exp))
}



## 4. calculate spectral exponent in each sliding window of varying widths for both linearly detrended data and seasonally detrended data:

l_spec_exp <- sliding_window_spec_exp(l_detrended_tas)
##saveRDS(l_spec_exp[[1]], "./data-processed/spec3d_year-time-step_l.rds")
write.csv(l_spec_exp[[2]], "./data-processed/cmcc-cesm_spectral-exponent-sliding-window_l.csv",
          row.names = FALSE)

s_spec_exp <- sliding_window_spec_exp(s_detrended_tas)
##saveRDS(s_spec_exp[[1]], "./data-processed/spec3d_year-time-step_s.rds")
write.csv(s_spec_exp[[2]], "./data-processed/cmcc-cesm_spectral-exponent-sliding-window_s.csv",
          row.names = FALSE)


spec3D <- s_spec_exp[[1]]

## try plotting:
n5 <- spec3D[,,1] ## extract change in spectral exponent over time windows with 5 year width

sp <-  expand.grid(lat, long)
colnames(sp) <- c("latitude", "longitude") 
sp$spec <- as.vector(n5)
sp <- sp %>%
  select(longitude, latitude, spec)

gg4 <- sp %>%
  ggplot() +
  geom_raster(aes(x = longitude, y = latitude, fill = spec)) + 
  ggtitle("Change in spectral exponent over years 1850-2100 (5 year window, 1 year step)") + 
  coord_fixed() +
  xlab("Longitude") +
  ylab("Latitude") +
  labs(fill = "Slope of spectral exponent")

## woohoo!!!!!
## before moving on to rest of the GCMs, lets try and calculate some velocities of variation! 




##### exploring sea surface temperature files ######
library(ncdf4)
library(raster)
filename <- "/Volumes/Nikki 6TB/velocities-of-variability/data-raw/01_CMCC-CESM/tos_day_CMCC-CESM_historical_r1i1p1_18680101-18731231.nc"
nc <- nc_open(filename)
tos = ncvar_get(nc, "tos")
lat = ncvar_get(nc, "lat")
lon = ncvar_get(nc, "lon")
lat_bnds = ncvar_get(nc, "lat_vertices")
lon_bnds = ncvar_get(nc, "lon_vertices")
time = ncvar_get(nc, "time")
nc_close(nc)

r <- raster(tos[,,1])
plot(r)
r



## get a feel for the regularity of the grid
## try plotting one layer of temps
element = 1
temps <- list()
for (i in 1:182) {
  for (j in 1:149) {
    temps[[element]] <- c(lon[i,j], lat[i,j], tos[i,j,1])
    element = element + 1
  }
}

df_temps <- data.frame(do.call(rbind, temps))
colnames(df_temps) <- c("lon", "lat", "temp")
df_temps <- filter(df_temps, !is.na(temp))

ggplot(df_temps, aes(y = lat, x = lon, colour = temp)) + geom_point(size = 0.0001) + 
  theme_minimal() + coord_fixed()

sp <- SpatialPointsDataFrame(coords = df_temps, data = df_temps, 
                             proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# ## save spatial points for use in QGIS 
# library(rgdal)
# writeOGR(obj=sp, dsn=".",
#          layer = "tos",
#          driver="ESRI Shapefile", 
#          overwrite = TRUE)


## make polygons for each grid cell
element = 1
i = 1
while (i < 183) {
  j = 1
  while (j < 150) {
    if (!is.na(tos[i,j,1])) {
      mat <- matrix(ncol = 2, c(c(lon_bnds[c(1:4,1),i,j],c(lat_bnds[c(1:4,1),i,j]))))
      p <- Polygon(mat)
      ps = Polygons(list(p),1)
      sps = SpatialPolygons(list(ps))
      #plot(sps)
      
      if (element == 1) {
        # create data frame
        df <- data.frame(lat = lat[i,j], lon = lon[i,j], tos = tos[i,j,1])
        sp_df <- SpatialPolygonsDataFrame(sps, data = df)
        element = element + 1
      }
      else {
        # add shape to df
        df <- data.frame(lat = lat[i,j], lon = lon[i,j], tos = tos[i,j,1])
        sp_df <- rbind(sp_df, SpatialPolygonsDataFrame(sps, data = df))
      }
    }
    j = j+1
  }
  i = i+1
}


library(ncdf4)
library(raster)
filename <- "/Volumes/Nikki 6TB/velocities-of-variability/data-raw/01_CMCC-CESM/regridded_tos_day_CMCC-CESM_historical_r1i1p1_18680101-18731231.nc"
nc <- nc_open(filename)
tos = ncvar_get(nc, "tos")
lat = ncvar_get(nc, "lat")
lon = ncvar_get(nc, "lon")
time = ncvar_get(nc, "time")
nc_close(nc)

r <- stack("/Volumes/Nikki 6TB/velocities-of-variability/data-raw/01_CMCC-CESM/regridded_tos_day_CMCC-CESM_historical_r1i1p1_18680101-18731231.nc")
plot(r[[1000]])
r

r = raster(nrows = 180, ncols = 360, xmn=-180, xmx=180, ymn=-90, ymx=90, 
           crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

values(r) <- as.vector(tos)


### making a land/ocean mask from the SST grids:
r1 <- stack("/Volumes/Nikki 6TB/velocities-of-variability/data-raw/01_CMCC-CESM/regridded_tos_day_CMCC-CESM_historical_r1i1p1_18680101-18731231.nc")[[1]]
r3 <- stack("/Users/nikkimoore/Desktop/regridded_tos_day_CMCC-CMS_historical_r1i1p1_18700101-18791231.nc")[[1]]

r4 <- stack("/Users/nikkimoore/Desktop/regridded_tos_day_MPI-ESM-LR_historical_r2i1p1_18700101-18791231.nc")[[1]]
r5 <- stack("/Users/nikkimoore/Desktop/regridded_tos_day_MPI-ESM-MR_historical_r1i1p1_18700101-18701231.nc")[[1]]
r6 <- stack("/Users/nikkimoore/Desktop/regridded_tos_day_GFDL-ESM2G_historical_r1i1p1_18660101-18701231.nc")[[1]]
r7 <- stack("/Users/nikkimoore/Desktop/regridded_tos_day_GFDL-CM3_historical_r1i1p1_18700101-18741231.nc")[[1]]

plot(r1)
plot(r3)
plot(r4)
plot(r5)
plot(r6)
plot(r7)

r <- mean(r1,r3,r4,r5,r6,r7, na.rm=T)
plot(r)
## these are the places we will definitely have sst estimates for

writeRaster(r, "data-processed/masks/cmip5-ocean.grd", 
            overwrite = TRUE, format = "raster")

land <- r
land[is.na(r)] <- 1
land[r == 1] <- NA
plot(land)

writeRaster(land, "data-processed/masks/cmip5-land.grd", 
            overwrite = TRUE, format = "raster")

## create list of ocean grid cell coordinates
ocean <- data.frame(rasterToPoints(r)[,1:2])
colnames(ocean) <- c("lon", "lat")

## transform coordinates from -180, 180 to 0, 360 to match temperature data
ocean$lon <- ifelse(ocean$lon <= 0, ocean$lon + 360, ocean$lon)

write.csv(ocean, "data-processed/masks/cmip5-ocean-coords.csv", row.names = F)

