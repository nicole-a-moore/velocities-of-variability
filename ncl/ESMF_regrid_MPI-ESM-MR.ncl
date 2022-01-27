;======================================================================
; ESMF_regrid_6.ncl
;
; Concepts illustrated:
;   - Interpolating from one grid to another using ESMF_regrid
;   - Interpolating data from a CMIP5 grid to a 1X1 degree rectilinear grid
;======================================================================
; This example is identical to ESMF_all_6.ncl, except it does the
; regridding in one call to "ESMF_regrid".  See ESMF_wgts_6.ncl
; for a faster example of regridding using an existing weights file.
;======================================================================
; This example uses the ESMF application "ESMF_RegridWeightGen" to 
; generate the weights.
;
; For more information about ESMF:
;
;        http://www.earthsystemmodeling.org/
;======================================================================
; This script regrids a CMIP5 grid to a 1.0 degree world grid and
; plots sea water potential temperature on the new grid.
;
; It uses SCRIP for both the CMIP5 and 1.0 degree world grid.
;======================================================================
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;
; This file still has to be loaded manually
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin
;--- make a list of all filenames:

	filepath  = "/Users/nikkimoore/Documents/velocities-of-variability/data-raw/05_MPI-ESM-MR/"
	filenames = ((/"tos_day_MPI-ESM-MR_historical_r1i1p1_18700101-18701231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18710101-18711231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18720101-18721231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18730101-18731231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18740101-18741231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18750101-18751231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18760101-18761231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18770101-18771231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18780101-18781231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18790101-18791231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18800101-18801231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18810101-18811231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18820101-18821231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18830101-18831231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18840101-18841231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18850101-18851231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18860101-18861231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18870101-18871231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18880101-18881231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18890101-18891231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18900101-18901231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18910101-18911231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18920101-18921231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18930101-18931231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18940101-18941231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18950101-18951231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18960101-18961231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18970101-18971231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18980101-18981231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_18990101-18991231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19000101-19001231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19010101-19011231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19020101-19021231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19030101-19031231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19040101-19041231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19050101-19051231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19060101-19061231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19070101-19071231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19080101-19081231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19090101-19091231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19100101-19101231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19110101-19111231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19120101-19121231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19130101-19131231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19140101-19141231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19150101-19151231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19160101-19161231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19170101-19171231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19180101-19181231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19190101-19191231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19200101-19201231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19210101-19211231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19220101-19221231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19230101-19231231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19240101-19241231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19250101-19251231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19260101-19261231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19270101-19271231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19280101-19281231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19290101-19291231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19300101-19301231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19310101-19311231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19320101-19321231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19330101-19331231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19340101-19341231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19350101-19351231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19360101-19361231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19370101-19371231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19380101-19381231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19390101-19391231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19400101-19401231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19410101-19411231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19420101-19421231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19430101-19431231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19440101-19441231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19450101-19451231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19460101-19461231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19470101-19471231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19480101-19481231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19490101-19491231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19500101-19501231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19510101-19511231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19520101-19521231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19530101-19531231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19540101-19541231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19550101-19551231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19560101-19561231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19570101-19571231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19580101-19581231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19590101-19591231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19600101-19601231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19610101-19611231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19620101-19621231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19630101-19631231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19640101-19641231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19650101-19651231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19660101-19661231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19670101-19671231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19680101-19681231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19690101-19691231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19700101-19701231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19710101-19711231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19720101-19721231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19730101-19731231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19740101-19741231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19750101-19751231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19760101-19761231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19770101-19771231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19780101-19781231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19790101-19791231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19800101-19801231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19810101-19811231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19820101-19821231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19830101-19831231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19840101-19841231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19850101-19851231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19860101-19861231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19870101-19871231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19880101-19881231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19890101-19891231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19900101-19901231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19910101-19911231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19920101-19921231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19930101-19931231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19940101-19941231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19950101-19951231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19960101-19961231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19970101-19971231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19980101-19981231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_19990101-19991231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_20000101-20001231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_20010101-20011231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_20020101-20021231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_20030101-20031231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_20040101-20041231.nc","tos_day_MPI-ESM-MR_historical_r1i1p1_20050101-20051231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20060101-20061231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20070101-20071231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20080101-20081231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20090101-20091231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20100101-20101231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20110101-20111231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20120101-20121231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20130101-20131231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20140101-20141231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20150101-20151231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20160101-20161231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20170101-20171231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20180101-20181231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20190101-20191231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20200101-20201231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20210101-20211231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20220101-20221231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20230101-20231231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20240101-20241231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20250101-20251231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20260101-20261231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20270101-20271231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20280101-20281231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20290101-20291231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20300101-20301231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20310101-20311231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20320101-20321231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20330101-20331231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20340101-20341231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20350101-20351231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20360101-20361231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20370101-20371231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20380101-20381231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20390101-20391231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20400101-20401231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20410101-20411231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20420101-20421231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20430101-20431231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20440101-20441231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20450101-20451231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20460101-20461231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20470101-20471231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20480101-20481231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20490101-20491231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20500101-20501231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20510101-20511231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20520101-20521231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20530101-20531231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20540101-20541231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20550101-20551231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20560101-20561231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20570101-20571231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20580101-20581231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20590101-20591231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20600101-20601231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20610101-20611231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20620101-20621231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20630101-20631231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20640101-20641231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20650101-20651231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20660101-20661231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20670101-20671231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20680101-20681231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20690101-20691231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20700101-20701231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20710101-20711231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20720101-20721231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20730101-20731231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20740101-20741231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20750101-20751231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20760101-20761231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20770101-20771231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20780101-20781231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20790101-20791231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20800101-20801231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20810101-20811231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20820101-20821231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20830101-20831231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20840101-20841231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20850101-20851231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20860101-20861231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20870101-20871231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20880101-20881231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20890101-20891231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20900101-20901231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20910101-20911231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20920101-20921231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20930101-20931231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20940101-20941231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20950101-20951231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20960101-20961231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20970101-20971231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20980101-20981231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_20990101-20991231.nc","tos_day_MPI-ESM-MR_rcp85_r1i1p1_21000101-21001231.nc"/))

;--- for each file in filenames, regrid:
	do n=0,dimsizes(filenames)-1
	
	;---Interpolation methods
    methods      = ("bilinear")

;---Input file
    srcFileName  = filepath + filenames(n)
    wgtFile      = filepath + "CMIP5_2_World_" + methods + ".nc"

;---Get data and lat/lon grid from CMIP5 Grid
    sfile        = addfile(srcFileName,"r")
    tos       = sfile->tos
    tos@lat2d = sfile->lat
    tos@lon2d = sfile->lon

    Opt                  = True
    Opt@SrcFileName      = filepath + "CMIP5_SCRIP.nc"       ; source file name
    Opt@DstFileName      = filepath + "World1deg_SCRIP.nc"   ; destination file name
    Opt@ForceOverwrite   = True

    Opt@SrcGridCornerLat = sfile->lat_vertices    ; corners are necessary
    Opt@SrcGridCornerLon = sfile->lon_vertices    ; for "conserve" method
    Opt@SrcMask2D        = where(.not.ismissing(tos(0,:,:)),1,0)
    
    Opt@DstGridType      = "1x1"              ; Destination grid
    Opt@DstTitle         = "World Grid 1-degree Resolution"
    Opt@DstLLCorner      = (/-89.5d,   -179.5d /)
    Opt@DstURCorner      = (/ 89.5d, 179.5d /) 

    ;;Opt@PrintTimings   = True
    ;;Opt@Debug          = True

;----------------------------------------------------------------------
; Setup for graphics
;----------------------------------------------------------------------
    wks = gsn_open_wks("png","ESMF_regrid")     ; send graphics to PNG file

;---Resources to share between both plots
    res                     = True              ; Plot modes desired.

    res@gsnDraw             = False             ; Will panel later
    res@gsnFrame            = False             ; Will panel later

    res@gsnMaximize         = True              ; Maximize plot

    res@cnFillOn            = True              ; color plot desired
    res@cnFillPalette       = "rainbow"         ; set color map
    res@cnLinesOn           = False             ; turn off contour lines
    res@cnLineLabelsOn      = False             ; turn off contour labels
    res@cnFillMode          = "RasterFill"      ; turn raster on      
  
    res@cnLevelSelectionMode = "ExplicitLevels"
    res@cnLevels             = ispan(270,300,2)
  
    res@mpFillOn            = False

    res@trGridType         = "TriangularMesh"  ; allow missing coordinates
    res@gsnAddCyclic       = False

    res@pmLabelBarWidthF   = 0.7
    res@lbLabelBarOn       = False   ; Will do this in panel

    res@gsnAddCyclic       = False

;---Resources for paneling
    pres                  = True
    pres@gsnMaximize      = True
    pres@gsnPanelLabelBar = True
    pres@lbLabelFontHeightF = 0.01

;----------------------------------------------------------------------
; Loop across each method and generate interpolation weights for 
; CMIP5 Grid to World Grid    
;----------------------------------------------------------------------
    plot_regrid = new(dimsizes(methods),graphic)

    do i=0,dimsizes(methods)-1
      print("Generating interpolation weights from CMIP5 to")
      print("World 1 degree grid using the " + methods(i) + " method.")

      Opt@WgtFileName  = wgtFile(i)
      Opt@InterpMethod = methods(i)

;----------------------------------------------------------------------
; Interpolate data from CMIP5 to World 1-degree grid.
;----------------------------------------------------------------------

      tos_regrid = ESMF_regrid(tos,Opt)
      printVarSummary(tos_regrid)
      
      filename = filepath + "regridded_" + filenames(n)
      
      ;-- first, delete if exists
      if(isfilepresent(filename)) then
		system("/bin/rm " + filename)
	  end if
       
	  ofile = addfile(filename,"c")          ;-- open file for creation
  	  ofile->tos  = tos_regrid  

;---Clean up before next time in loop.
      delete(tos_regrid)
    end do
	
	delete(tos)
	end do
end