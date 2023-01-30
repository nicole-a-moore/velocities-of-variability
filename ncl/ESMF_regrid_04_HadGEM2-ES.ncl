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

	filepath  = "/Volumes/NIKKI/CMIP5-GCMs/04_HadGEM2-ES_tos/"
	filenames = ((/"tos_day_HadGEM2-ES_historical_r1i1p1_18691201-18791130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_18791201-18891130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_18891201-18991130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_18991201-19091130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19091201-19191130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19191201-19291130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19291201-19391130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19391201-19491130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19491201-19591130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19591201-19691130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19691201-19791130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19791201-19891130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19891201-19991130.nc","tos_day_HadGEM2-ES_historical_r1i1p1_19991201-20051130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20051201-20101130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20101201-20151130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20151201-20201130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20201201-20251130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20251201-20301130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20301201-20351130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20351201-20401130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20401201-20451130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20451201-20501130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20501201-20551130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20551201-20601130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20601201-20651130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20651201-20701130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20701201-20751130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20751201-20851130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20851201-20991130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_20991201-21091130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21091201-21191130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21191201-21291130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21291201-21391130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21391201-21491130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21491201-21591130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21591201-21691130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21691201-21791130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21791201-21891130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21891201-21991130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_21991201-22091130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22091201-22191130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22191201-22291130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22291201-22391130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22391201-22491130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22491201-22591130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22591201-22691130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22691201-22791130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22791201-22891130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22891201-22991130.nc","tos_day_HadGEM2-ES_rcp85_r1i1p1_22991201-22991230.nc"/))

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

    Opt@SrcGridCornerLat = sfile->lat_bnds    ; corners are necessary
    Opt@SrcGridCornerLon = sfile->lon_bnds    ; for "conserve" method
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