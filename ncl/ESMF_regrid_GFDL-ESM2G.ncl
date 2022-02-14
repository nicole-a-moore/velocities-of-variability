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

	filepath  = "/Users/nikkimoore/Documents/velocities-of-variability/data-raw/08_GFDL-ESM2M/"
	filenames = ((/"tos_day_GFDL-ESM2G_historical_r1i1p1_18660101-18701231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_18710101-18751231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_18760101-18801231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_18810101-18851231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_18860101-18901231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_18910101-18951231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_18960101-19001231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19010101-19051231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19060101-19101231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19110101-19151231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19160101-19201231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19210101-19251231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19260101-19301231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19310101-19351231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19360101-19401231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19410101-19451231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19460101-19501231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19510101-19551231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19560101-19601231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19610101-19651231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19660101-19701231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19710101-19751231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19760101-19801231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19810101-19851231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19860101-19901231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19910101-19951231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_19960101-20001231.nc","tos_day_GFDL-ESM2G_historical_r1i1p1_20010101-20051231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20060101-20101231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20110101-20151231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20160101-20201231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20210101-20251231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20260101-20301231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20310101-20351231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20360101-20401231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20410101-20451231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20460101-20501231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20510101-20551231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20560101-20601231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20610101-20651231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20660101-20701231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20710101-20751231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20760101-20801231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20810101-20851231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20860101-20901231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20910101-20951231.nc","tos_day_GFDL-ESM2G_rcp85_r1i1p1_20960101-21001231.nc"/))

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