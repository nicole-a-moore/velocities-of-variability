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

	filepath  = "/Users/nikkimoore/Documents/velocities-of-variability/data-raw/04_MPI-ESM-LR/"
	filenames = ((/"tos_day_MPI-ESM-LR_historical_r2i1p1_18600101-18691231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_18700101-18791231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_18800101-18891231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_18900101-18991231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19000101-19091231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19100101-19191231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19200101-19291231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19300101-19391231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19400101-19491231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19500101-19591231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19600101-19691231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19700101-19791231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19800101-19891231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_19900101-19991231.nc","tos_day_MPI-ESM-LR_historical_r2i1p1_20000101-20051231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20060101-20091231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20100101-20191231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20200101-20291231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20300101-20391231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20400101-20491231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20500101-20591231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20600101-20691231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20700101-20791231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20800101-20891231.nc","tos_day_MPI-ESM-LR_rcp85_r2i1p1_20900101-21001231.nc"/))

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