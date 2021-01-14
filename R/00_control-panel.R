## control panel: a place to execute all of the functions for each of the 21 GCMs
library(tidyverse)



## IMPORTANT: 
## set 'path' to where you have the GCM files stored on your computer
## for me, they are here:
path = "/Volumes/" ## change me

## create vector of file folders to put data into:
nums <- c(ifelse(nchar(seq(1:4)) == 2, seq(1:21), paste("0", seq(1:21), sep = "")))

gcm_models <- c("CMCC-CESM", "CMCC-CM", "CMCC-CMS", "inmcm4")

folders <- c(paste(nums, gcm_models, sep = "_")) %>%
  paste(path, "ADATA HV620/CMIP5-GCMs/", ., "/", sep = "")

## make list of filenames:

## 1870-01-01:
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

CMCC_CESM <- list(CMCC_CESM_hist, CMCC_CESM_rcp85)

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

CMCC_CM <- list(CMCC_CM_hist, CMCC_CM_rcp85)

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

CMCC_CMS <- list(CMCC_CMS_hist, CMCC_CMS_rcp85)

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

inmcm4 <- list(inmcm4_hist, inmcm4_rcp85)

gcm_files <- list(CMCC_CESM, CMCC_CM, CMCC_CMS, inmcm4)

## 1869-12-01:
HadGEM_AO_hist <- c("tas_day_HadGEM2-AO_historical_r1i1p1_18600101-20051230.nc")
HadGEM_AO_rcp85 <- c("tas_day_HadGEM2-AO_rcp85_r1i1p1_20060101-21001230.nc")



###############################################
###           READING IN FUNCTIONS          ###
###############################################
source("R/02_GCM-file-functions.R")
source("R/03_spectral-exponent-functions.R")
source("R/04_velocity-functions.R")

## loop through each GCM, calling functions that extract, organize, and detrend data before calculating 
## spectral exponent and spectral velocity vectors

## WARNING: 
##    - each GCM will take multiple days to complete (estimate: )
##    - to make sure your computer doesn't run out of memory, have at least 4GB memory free

gcm <- 1
while (gcm < length(gcm_models) + 1) {
  print(paste("Now starting GCM number ", gcm, ": ", gcm_models[gcm], "...", sep = ""))
  
  ###############################################
  ###       EXRACTING AND ORGANIZING DATA     ###
  ###############################################
  ## extract and organize the CMCC-CESM data into spatial chunks:
  files <- gcm_files[[gcm]]
  e_and_o <- extract_and_organize_memory_conscious(path = folders[gcm], historical_filenames = files[[1]], 
                                  rcp85_filenames = files[[2]])
  
  ###############################################
  ###               DETRENDING DATA           ###
  ###############################################
  ## detrend the CMCC-CESM data:
  d_tas <- detrend_tas(e_and_o)
  
  ###############################################
  ###     CALCULATING SPECTRAL EXPONENTS      ###
  ###############################################
  ## calculate spectral exponent over windows of varying widths:
  se <- sliding_window_spec_exp(d_tas)
  
  
  ###############################################
  ###       REORGANIZING INTO RASTERSTACKS    ###
  ###############################################
  ## reformat into rasterStacks for use in VoCC functions:
  stacks <- create_rasterStack(se)
  
  
  ###############################################
  ###       CALCULATING VELOCITY VECTORS      ###
  ###############################################
  ## calculate long term trend in spectral exponent:
  tt <- calculate_tempTrend(stacks)
  
  ## calculate spatial gradient in spectral exponent:
  sg <- calculate_spatGrad(stacks)
  
  ## calculate local velocity vectors:
  vocc <- calculate_VoCC(tt = tt, sg = sg)
  
  
  ###############################################
  ###           SAVING VELOCITY VECTORS       ###
  ###############################################
  filename <- paste("data-processed/VoCCs/", nums[gcm], gcm_models[gcm], ".rds", sep = "")
  saveRDS(vocc, filename)
  
  print(paste("Woohoo! Finally finished GCM: ", gcms[gcm], sep = ""))
  gcm = gcm + 1
}





