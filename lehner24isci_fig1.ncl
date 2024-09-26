load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
; ================================================;
begin

; ====================================================
; lehner24isci_fig1.ncl
; ====================================================
; Code to produce maps in Fig. 1 in Lehner (2024)
;
; For code to produce the Hawkins&Sutton-type plots in Fig. 1,
; please consult https://github.com/flehner/lehner20esd_code/blob/main/other_figs/lehner20esd_figs.m
;
; Citation:
; Lehner, F. (2024):
; Climate model large ensembles as test beds for applied compound event research
; iScience, DOI: TBD
;
; Notes on code:
; - Only does calculations and plotting, not pre-processing
; - Requires pre-processed data, which is not provided due to data volume:
;   - gridded monthly temperature and precipitation fields from CMIP6 models (CMIP6), conservatively regridded to 2.5° horizontal resolution
; ====================================================

  print_clock("start script")

  pathin  = "/project/ccr/flehner/cmip6-ng/"
  pathout = "/project/ccr/flehner/reservoir_model/cmip6-ng/"

  recalc1 = "yes" ; calc water year mean and anomalies and percentile thresholds for tas and pr
  recalc2 = "yes" ; calc hot-dry frequency (hdf)

  models  = (/"ACCESS-CM2","ACCESS-ESM1-5","AWI-CM-1-1-MR","BCC-CSM2-MR","CAMS-CSM1-0","CESM2","CESM2-WACCM","CNRM-CM6-1","CNRM-CM6-1-HR","CNRM-ESM2-1","CanESM5","CanESM5-CanOE","EC-Earth3","EC-Earth3-Veg","FGOALS-f3-L","FGOALS-g3","GFDL-ESM4","GISS-E2-1-G","INM-CM4-8","INM-CM5-0","IPSL-CM6A-LR","MCM-UA-1-0","MIROC-ES2L","MIROC6","MPI-ESM1-2-HR","MPI-ESM1-2-LR","MRI-ESM2-0","NorESM2-MM","UKESM1-0-LL"/)

  vars      = (/"tas","pr"/)

  ; -- fix parameters for file reading and concatenating --
  histstart   = 185001
  histende    = 201412
  sspstart    = 201501
  sspende     = 210012
  start       = round(histstart/100,3)
  ende        = round(sspende/100,3)
  ; -- choose epochs for difference maps --
  refstart    = 1995
  refende     = 2014
  istart      = refstart-start
  iende       = refende-start
  futstart    = 2051
  futende     = 2100
  jstart      = futstart-start
  jende       = futende-start
  ; -- percentile threshold for compound event
  p1          = 80; percentile threshold for temperature
  p2          = 20; percentile threshold for precipitation

  region      = "north_america_lambert" ; north_america north_america_lambert
  if region .eq. "north_america"
    minlat      = 22
    maxlat      = 74
    minlon      = -170
    maxlon      = -70
  end if
  if region .eq. "north_america_lambert"
    minlat      = 18
    maxlat      = 63
    minlon      = -125
    maxlon      = -72
    p1f         = 33.0
    p2f         = 45.0
    mf          = -105.
  end if

  ; -- create matrix
  ntime     = round(sspende/100,3)-round(histstart/100,3)+1
  nlat      = 72
  nlon      = 144

  scen = "ssp585" ; ssp126 ssp245 ssp370 ssp585

  if recalc1 .eq. "yes"
    var1_wy   = new((/ntime,dimsizes(models),nlat,nlon/),float)
    var2_wy   = var1_wy
    var1_wy_anom = var1_wy
    var2_wy_anom = var1_wy
    do m = 0,dimsizes(models)-1
      print_clock("model = "+models(m) +" | "+m+"/"+dimsizes(models))
      vari = "tas"
      print_clock("  "+vari)
      ifile = systemfunc("ls "+pathin+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p*f*_g025.nc")
      a     = addfile(ifile(0),"r")
      tmp0  = a->$vari$
      tim   = a->time
      cld   = cd_calendar(tim, -1)
      i_start = ind(cld .eq. histstart)
      i_ende = ind(cld .eq. sspende)
      ; -- cut to common length
      tmp1 = tmp0(i_start:i_ende,:,:)
      ; -- calc water year (wy)
      do i = 1,ntime-1
        var1_wy(i,m,:,:) = tofloat( dim_avg_n_Wrap(tmp1(((i-1)*12)+9:(i*12)+8,:,:),0) )
      end do
      copy_VarCoords(tmp1(0,:,:),var1_wy(0,0,:,:))
      var1_wy!0 = "time"
      ; -- calc anomaly to reference period
      clim0 = dim_avg_n_Wrap(var1_wy(refstart-start:refende-start,m,:,:),0)
      var1_wy_anom(:,m,:,:) = var1_wy(:,m,:,:) - conform_dims(dimsizes(var1_wy(:,m,:,:)),clim0,(/1,2/))
      copy_VarCoords(tmp1(0,:,:),var1_wy_anom(0,0,:,:))
      var1_wy_anom!0 = "time"
      ; -- save individual model data
      system("mkdir -p "+pathout+"/"+vari)
      ofile   = pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc"
      system("rm -rf "+ofile)
      a         = addfile(ofile,"c")
      a->tas_wy = rm_single_dims(var1_wy(:,m,:,:))
      a->tas_wy_anom = rm_single_dims(var1_wy_anom(:,m,:,:))
      delete([/tmp0,tmp1,tim,cld,i_start,i_ende,ifile/])
      ; -- use CDO to calculate percentiles outside of NCL (NCL too slow)
      system("rm -rf "+pathout+"/"+vari+"/tmp1.nc")
      system("rm -rf "+pathout+"/"+vari+"/tmp2.nc")
      system("rm -rf "+pathout+"/"+vari+"/tmp3.nc")
      system("rm -rf "+pathout+"/"+vari+"/tmp4.nc")
      system("cdo -seltimestep,"+istart+"/"+iende+" "+pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc "+pathout+"/"+vari+"/tmp1.nc")
      system("cdo -timmin "+pathout+"/"+vari+"/tmp1.nc "+pathout+"/"+vari+"/tmp2.nc")
      system("cdo -seltimestep,"+istart+"/"+iende+" "+pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc "+pathout+"/"+vari+"/tmp3.nc")
      system("cdo -timmax "+pathout+"/"+vari+"/tmp3.nc "+pathout+"/"+vari+"/tmp4.nc")
      system("cdo timpctl,"+p1+" "+pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc "+pathout+"/"+vari+"/tmp2.nc "+pathout+"/"+vari+"/tmp4.nc "+pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy_p"+p1+".nc")

      vari = "pr"
      print_clock("  "+vari)
      ifile = systemfunc("ls "+pathin+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p*f*_g025.nc")
      a     = addfile(ifile(0),"r")
      tmp0  = a->$vari$
      tim   = a->time
      cld   = cd_calendar(tim, -1)
      i_start = ind(cld .eq. histstart)
      i_ende = ind(cld .eq. sspende)
      ; -- cut to common length
      tmp1 = tmp0(i_start:i_ende,:,:)
      ; -- calc water year (wy)
      do i = 1,ntime-1
        var2_wy(i,m,:,:) = tofloat( dim_avg_n_Wrap(tmp1(((i-1)*12)+9:(i*12)+8,:,:),0) )
      end do
      copy_VarCoords(tmp1(0,:,:),var2_wy(0,0,:,:))
      var2_wy!0 = "time"
      ; -- calc anomaly to reference period
      clim0 = dim_avg_n_Wrap(var2_wy(refstart-start:refende-start,m,:,:),0)
      var2_wy_anom(:,m,:,:) = (( var2_wy(:,m,:,:) / conform_dims(dimsizes(var2_wy(:,m,:,:)),clim0,(/1,2/)) ) -1 ) *100; pr expressed as percent anomaly
      copy_VarCoords(tmp1(0,:,:),var2_wy_anom(0,0,:,:))
      var2_wy_anom!0 = "time"
      ; -- save individual model data
      system("mkdir -p "+pathout+"/"+vari)
      ofile   = pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc"
      system("rm -rf "+ofile)
      a         = addfile(ofile,"c")
      a->pr_wy = rm_single_dims(var2_wy(:,m,:,:))
      a->pr_wy_anom = rm_single_dims(var2_wy_anom(:,m,:,:))
      delete([/tmp0,tmp1,tim,cld,i_start,i_ende,ifile/])
      ; -- use CDO to calculate percentiles outside of NCL (NCL too slow)
      system("rm -rf "+pathout+"/"+vari+"/tmp1.nc")
      system("rm -rf "+pathout+"/"+vari+"/tmp2.nc")
      system("rm -rf "+pathout+"/"+vari+"/tmp3.nc")
      system("rm -rf "+pathout+"/"+vari+"/tmp4.nc")
      system("cdo -seltimestep,"+istart+"/"+iende+" "+pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc "+pathout+"/"+vari+"/tmp1.nc")
      system("cdo -timmin "+pathout+"/"+vari+"/tmp1.nc "+pathout+"/"+vari+"/tmp2.nc")
      system("cdo -seltimestep,"+istart+"/"+iende+" "+pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc "+pathout+"/"+vari+"/tmp3.nc")
      system("cdo -timmax "+pathout+"/"+vari+"/tmp3.nc "+pathout+"/"+vari+"/tmp4.nc")
      system("cdo timpctl,"+p2+" "+pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc "+pathout+"/"+vari+"/tmp2.nc "+pathout+"/"+vari+"/tmp4.nc "+pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy_p"+p2+".nc")

    end do ; end of model loop

    ; -- save data
    print("save data...")
    ; -- tas
    vari = "tas"
    system("mkdir -p "+pathout+"/"+vari)
    ofile   = pathout+"/"+vari+"/"+vari+"_mon_cmip6-ng_"+scen+"_r1i1p1f1_g025_wy.nc"
    system("rm -rf "+ofile)
    a         = addfile(ofile,"c")
    a->tas_wy = var1_wy
    a->tas_wy_anom = var1_wy_anom
    ; -- pr
    vari = "pr"
    system("mkdir -p "+pathout+"/"+vari)
    ofile   = pathout+"/"+vari+"/"+vari+"_mon_cmip6-ng_"+scen+"_r1i1p1f1_g025_wy.nc"
    system("rm -rf "+ofile)
    a         = addfile(ofile,"c")
    a->pr_wy = var2_wy
    a->pr_wy_anom = var2_wy_anom
  end if


  if recalc2 .eq. "yes"
    hdf_wy   = new((/ntime,dimsizes(models),nlat,nlon/),float)
    ; -- calc HDf --
    do m = 0,dimsizes(models)-1
      print_clock("model = "+models(m) +" | "+m+"/"+dimsizes(models))
      vari = "hdf"
      print_clock("  "+vari)
      ifile = systemfunc("ls "+pathout+"/tas/tas_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc")
      a     = addfile(ifile,"r")
      tas   = a->tas_wy
      ifile = systemfunc("ls "+pathout+"/tas/tas_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy_p"+p1+".nc")
      a     = addfile(ifile,"r")
      tas_p = conform_dims(dimsizes(tas),rm_single_dims(a->tas_wy),(/1,2/))
      tas_m = where(tas .gt. tas_p,1,0); masked
      ifile = systemfunc("ls "+pathout+"/pr/pr_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc")
      a     = addfile(ifile,"r")
      pr    = a->pr_wy
      ifile = systemfunc("ls "+pathout+"/pr/pr_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy_p"+p2+".nc")
      a     = addfile(ifile,"r")
      pr_p  = conform_dims(dimsizes(tas),rm_single_dims(a->pr_wy),(/1,2/))
      pr_m  = where(pr .lt. pr_p,1,0); masked
      hdf   = where(tas_m .eq. 1 .and. pr_m .eq. 1,1,0)
      copy_VarCoords(tas,hdf)
      ; -- save individual model data
      system("mkdir -p "+pathout+"/"+vari)
      ofile   = pathout+"/"+vari+"/"+vari+"_mon_"+models(m)+"_"+scen+"_r1i1p1f1_g025_wy.nc"
      system("rm -rf "+ofile)
      a         = addfile(ofile,"c")
      a->hdf    = hdf
      ; --
      hdf_wy(:,m,:,:) = hdf
    end do
    system("mkdir -p "+pathout+"/"+vari)
    ofile   = pathout+"/"+vari+"/"+vari+"_mon_cmip6-ng_"+scen+"_r1i1p1f1_g025_wy.nc"
    system("rm -rf "+ofile)
    a         = addfile(ofile,"c")
    a->hdf_wy = hdf_wy
  end if


  ; -- load data and calculate epoch differences
  vari    = "tas"
  ifile   = pathout+"/"+vari+"/"+vari+"_mon_cmip6-ng_"+scen+"_r1i1p1f1_g025_wy.nc"
  a       = addfile(ifile,"r")
  tas_wy  = a->tas_wy

  vari    = "pr"
  ifile   = pathout+"/"+vari+"/"+vari+"_mon_cmip6-ng_"+scen+"_r1i1p1f1_g025_wy.nc"
  a       = addfile(ifile,"r")
  ; pr_wy   = a->pr_wy
  pr_wy_anom   = a->pr_wy_anom

  vari    = "hdf"
  ifile   = pathout+"/"+vari+"/"+vari+"_mon_cmip6-ng_"+scen+"_r1i1p1f1_g025_wy.nc"
  a       = addfile(ifile,"r")
  hdf_wy  = a->hdf_wy
  hdf_wy  = hdf_wy * 100 ; convert from fraction to percent

  ; -- mean change
  tas_change = dim_avg_n_Wrap(tas_wy(jstart:jende,:,:,:),0) - dim_avg_n_Wrap(tas_wy(istart:iende,:,:,:),0)
  pr_change = dim_avg_n_Wrap(pr_wy_anom(jstart:jende,:,:,:),0) - dim_avg_n_Wrap(pr_wy_anom(istart:iende,:,:,:),0)
  hdf_change = (dim_avg_n_Wrap(hdf_wy(jstart:jende,:,:,:),0) - dim_avg_n_Wrap(hdf_wy(istart:iende,:,:,:),0))
  ; -- multi-model mean (mmm) change
  tas_mmm_change = dim_avg_n_Wrap(dim_avg_n_Wrap(tas_wy(jstart:jende,:,:,:),0),0) - dim_avg_n_Wrap(dim_avg_n_Wrap(tas_wy(istart:iende,:,:,:),0),0)
  pr_mmm_change = dim_avg_n_Wrap(dim_avg_n_Wrap(pr_wy_anom(jstart:jende,:,:,:),0),0) - dim_avg_n_Wrap(dim_avg_n_Wrap(pr_wy_anom(istart:iende,:,:,:),0),0)
  ; -- hdf change
  hdf_abs_or_percent = 0 ; 0=absolute change, 1=percent change
  hdf_clim = dim_avg_n_Wrap(dim_avg_n_Wrap(hdf_wy(istart:iende,:,:,:),0),0)
  if hdf_abs_or_percent .eq. 0
    hdf_mmm_change = (dim_avg_n_Wrap(dim_avg_n_Wrap(hdf_wy(jstart:jende,:,:,:),0),0) - dim_avg_n_Wrap(dim_avg_n_Wrap(hdf_wy(istart:iende,:,:,:),0),0))
  else
    tmp = dim_avg_n_Wrap(dim_avg_n_Wrap(hdf_wy(istart:iende,:,:,:),0),0)
    tmp = where(tmp.eq.0,tmp@_FillValue,tmp)
    hdf_mmm_change = (dim_avg_n_Wrap(dim_avg_n_Wrap(hdf_wy(jstart:jende,:,:,:),0),0) - dim_avg_n_Wrap(dim_avg_n_Wrap(hdf_wy(istart:iende,:,:,:),0),0)) / tmp
    delete(tmp)
  end if
  copy_VarCoords(tas_wy(0,0,:,:),tas_mmm_change)
  copy_VarCoords(pr_wy_anom(0,0,:,:),pr_mmm_change)
  copy_VarCoords(hdf_wy(0,0,:,:),hdf_mmm_change)

  ; -- agreement on sign (aos) of change - as fraction of models
  ; -- tas
  mean_change_binary = where(tas_change .gt. 0, 1, 0)
  mmm_change_binary = where(tas_mmm_change .gt. 0, 1, 0)
  aos_tmp = where(mean_change_binary .eq. conform_dims(dimsizes(mean_change_binary), mmm_change_binary, (/1,2/)), 1, 0)
  aos_sum = dim_sum_n_Wrap(aos_tmp,0)
  tmp = new((/dimsizes(aos_sum(:,0)),dimsizes(aos_sum(0,:))/),float)
  tmp(:,:) = dimsizes(models)
  tas_aos = aos_sum/tmp
  delete([/tmp,aos_sum,aos_tmp/])
  copy_VarCoords(tas_wy(0,0,:,:),tas_aos(:,:))
  ; -- pr
  mean_change_binary = where(pr_change .gt. 0, 1, 0)
  mmm_change_binary = where(pr_mmm_change .gt. 0, 1, 0)
  aos_tmp = where(mean_change_binary .eq. conform_dims(dimsizes(mean_change_binary), mmm_change_binary, (/1,2/)), 1, 0)
  aos_sum = dim_sum_n_Wrap(aos_tmp,0)
  tmp = new((/dimsizes(aos_sum(:,0)),dimsizes(aos_sum(0,:))/),float)
  tmp(:,:) = dimsizes(models)
  pr_aos = aos_sum/tmp
  delete([/tmp,aos_sum,aos_tmp/])
  copy_VarCoords(pr_wy_anom(0,0,:,:),pr_aos(:,:))
  ; -- hdf
  mean_change_binary = where(hdf_change .gt. 0, 1, 0)
  mmm_change_binary = where(hdf_mmm_change .gt. 0, 1, 0)
  aos_tmp = where(mean_change_binary .eq. conform_dims(dimsizes(mean_change_binary), mmm_change_binary, (/1,2/)), 1, 0)
  aos_sum = dim_sum_n_Wrap(aos_tmp,0)
  tmp = new((/dimsizes(aos_sum(:,0)),dimsizes(aos_sum(0,:))/),float)
  tmp(:,:) = dimsizes(models)
  hdf_aos = aos_sum/tmp
  delete([/tmp,aos_sum,aos_tmp/])
  copy_VarCoords(pr_wy_anom(0,0,:,:),hdf_aos(:,:))

  ; -- magnitude (mag) of change
  tas_mag = where(tas_mmm_change .gt. 1*dim_stddev_n_Wrap(tas_change,0), 1, 0)
  copy_VarCoords(tas_wy(0,0,:,:),tas_mag(:,:))
  pr_mag = where(pr_mmm_change .gt. 1*dim_stddev_n_Wrap(pr_change,0), 1, 0)
  copy_VarCoords(pr_wy_anom(0,0,:,:),pr_mag(:,:))
  hdf_mag = where(hdf_mmm_change .gt. 1*dim_stddev_n_Wrap(hdf_change,0), 1, 0)
  copy_VarCoords(pr_wy_anom(0,0,:,:),hdf_mag(:,:))




  ;**************************************************************************
  ; PLOTTING
  print_clock("plotting...")

  varp 		    = "/home/flehner/publication/lehner22_smile_perspective/fig/cmip6_mmm_change_wy_"+scen+"_"+region+"_ref"+refstart+"-"+refende
  wks_type 	  = "eps"
  wks 		    = gsn_open_wks(wks_type,varp)
  nbplot      = 3
  plot 		    = new(nbplot,graphic)
  plot_stip1  = new(nbplot,graphic)
  plot_stip2  = new(nbplot,graphic)

  ; -- merge Adam's colormap with precip_diff_12lev
  gsn_define_colormap(wks,"amwg_blueyellowred")
  cmap = gsn_retrieve_colormap(wks)
  hsv = rgbhsv(cmap)
  hsv(10,1) = .47
  hsv(8,1) = .45    ; 1c+1d
  hsv(9,1) = .30    ; 1d
  hsv(4,:) = (/ (hsv(3,:)+hsv(5,:))/2.0 /)
  cmap2 = hsvrgb(hsv)
  gsn_define_colormap(wks,cmap2)
  newcolor = NhlNewColor(wks,0.42,0.42,0.42)
  nc0 = NhlNewColor(wks,0.95,0.95,0.95)
  nc1 = NhlNewColor(wks,0.7,0.7,0.7)
  nc3 = NhlNewColor(wks,0.5,0.5,0.5)
  nc4 = NhlNewColor(wks,0.30,0.30,0.30)
  nc2 = NhlNewColor(wks,139/256.,35/256.,35/256.)
  gsn_merge_colormaps(wks,cmap2,"precip_diff_12lev")

  res                     = True
  proj = "LambertConformal"; "Robinson"
  res@mpLambertParallel1F = p1f         ; two parallels
  res@mpLambertParallel2F = p2f
  res@mpLambertMeridianF  = mf        ; central meridian
  res@mpProjection        = proj
  res@mpFillOn = True                       ; default = True
  res@mpFillColors = (/-1,0,-1,0/)      ; map color fill ocean
  res@cnFillDrawOrder = "PreDraw"        ; color fill first
  res@mpOutlineDrawOrder = "Draw"      ; draw map outline after cnFill
  res@mpFillDrawOrder = "Draw"              ; draw map fill next
  res@mpPerimDrawOrder = "PostDraw"   ; draw map perimeter last
  res@gsnDraw              = False    ; do not draw picture
  res@gsnFrame             = False    ; do not advance frame
  res@gsnPaperOrientation  = "landscape"
  res@gsnAddCyclic = False;True
  res@mpLimitMode   = "LatLon"
  res@mpOutlineOn = True
  res@mpDataSetName = "Earth..4"
  res@mpDataBaseVersion = "MediumRes"
  res@mpGeophysicalLineColor = "gray42"
  res@mpCountyLineColor = "transparent"
  res@mpUSStateLineColor = "gray42"
  res@mpNationalLineColor = "gray22"
  res@mpGeophysicalLineThicknessF = 1.4
  res@mpOutlineBoundarySets = "National" ; "AllBoundaries" "GeophysicalAndUSStates";"USStates";"National"

  res@mpMinLatF            = minlat
  res@mpMaxLatF            = maxlat
  res@mpMinLonF            = minlon
  res@mpMaxLonF            = maxlon
  ; res@mpCenterLonF         = 255.         ; default is 0 [GM]

; Define colors for filled contours
  res@cnFillOn             = True
  res@gsnSpreadColors      = True                  ; use full colormap
  res@cnLevelSelectionMode = "ExplicitLevels"

  res@cnLinesOn           = False         ; Turn lines off
  res@cnLineLabelsOn      = False         ; Turn labels off
  res@gsnLeftString = ""
  res@gsnRightString = ""
  res@mpGridLineColor      =  -1       ; Turn off grid lat/lon lines
  res@mpGridAndLimbOn      =  True
  res@lbLabelBarOn        = True; False
  res@lbBoxEndCapStyle    = "TriangleBothEnds"; "TriangleLowEnd", "TriangleHighEnd", "TriangleBothEnds"
  res@lbTitlePosition         = "Bottom"
  res@lbLabelFontHeightF      = 0.02
  res@lbTitleFontHeightF      = 0.02
  res@pmLabelBarHeightF       = 0.08
  res@cnMissingValFillPattern = "SolidFill"
  res@cnMissingValFillColor =-1
  res@tmXBOn               = False
  res@tmXTOn               = False
  res@tmYLOn               = False
  res@tmYROn               = False

  ; -- tas resources
  rest = res
  rest@cnFillColors     = (/10,11,12,13,14,15,16/)
  rest@cnLevels         = fspan(2.5,5,6)
  rest@lbTitleString    = "Temperature change (~F34~0~F~C)"
  ; -- pr resources
  resp = res
  delete(resp@cnFillColors)
  resp@cnFillColors     = (/18,19,20,21,22,23,25,26,27,28,29,30/)
  resp@cnLevels         = fspan(-25,25,11)
  resp@lbTitleString    = "Precipitation change (%)"
  ; -- hdf resources
  resm = resp
  resm@cnFillPalette             = "BlueDarkOrange18"; "GreenMagenta16"
  delete(resm@cnFillColors)
  delete(resm@cnLevels)
  resm@cnFillColors     = (/9,10,11,12,13,14,15,16/)
  if hdf_abs_or_percent .eq. 0
    resm@cnLevels         = (/0,5,10,15,20,25/)
    resm@lbTitleString    = "fHD change (percentage points)"
  else
    resm@cnLevels         = (/0,5,10,15,20,25/); fspan(10,100,6)
    resm@lbTitleString    = "fHD change (%)"
  end if

  res2 = res
  res2@cnLinesOn           = True         ; Turn lines off
  res2@cnLineColor        = "blue";"gray10"
  res2@cnLineThicknessF      = 1.5
  res2@cnFillOn            = False
  res2@cnInfoLabelOn       = False
  res2@cnLevelSelectionMode = "ManualLevels"      ; manually set cn levels
  delete(res2@cnLevels)
  res2@cnMinLevelValF       = -10.25                 ; min level
  res2@cnMaxLevelValF       =  10.25                 ; max level
  res2@cnLevelSpacingF      =  .5                  ; contour level spacing
  res2@gsnContourNegLineDashPattern = 2

  res3=True
  res3@gsnDraw              = False        ; don't draw yet
  res3@gsnFrame             = False
  res3@cnFillOn            = True
  res3@cnLinesOn           = False
  res3@cnLineLabelsOn      = False
  res3@gsnLeftString = ""
  res3@cnLineThicknessF      = 3
  res3@cnMonoFillColor       = True
  res3@cnMonoFillPattern     = False
  res3@cnMonoFillScale       = True
  res3@cnFillColor = "black"
  res3cnFillDotSizeF = 0.003
  res3@cnLevelSelectionMode  = "ExplicitLevels"
  res3@cnFillPatterns        = (/8,-1/)
  res3@cnFillScaleF = 0.5
  res3@cnLevels              = fspan(-3,3,7)
  res3@cnInfoLabelOn         = False
  res3@lbLabelBarOn        = False

  res@tiMainFontHeightF     = 0.015
  res@tiYAxisFontHeightF    = 0.015

; ---------------------------------------------------------------------

  dum = new(nbplot,graphic)

  plot(0) 	= gsn_csm_contour_map_ce(wks,tas_mmm_change,rest)
  plot(1) 	= gsn_csm_contour_map_ce(wks,pr_mmm_change,resp)
  plot(2) 	= gsn_csm_contour_map_ce(wks,hdf_mmm_change,resm)


  ; -- aos
  ; -- overlay hatching
  resov1       = True
  resov1@mpProjection        = proj
  resov1@mpLambertParallel1F = 33.0         ; two parallels
  resov1@mpLambertParallel2F = 45.0
  resov1@mpLambertMeridianF  = -95.0        ; central meridian
  resov1@mpFillColors = (/-1,0,-1,0/)      ; map color fill ocean
  resov1@cnFillDrawOrder = "PreDraw"        ; color fill first
  resov1@mpOutlineDrawOrder = "Draw"      ; draw map outline after cnFill
  resov1@mpFillDrawOrder = "Draw"              ; draw map fill next
  resov1@mpPerimDrawOrder = "PostDraw"   ; draw map perimeter last
  resov1@gsnDraw              = False        ; don't draw yet
  resov1@gsnFrame             = False
  resov1@cnFillOn            = True
  resov1@cnLinesOn           = False
  resov1@cnLineLabelsOn      = False
  resov1@gsnLeftString = ""
  resov1@cnLineThicknessF      = 3
  resov1@cnMonoFillColor       = True
  resov1@cnMonoFillPattern     = False
  resov1@cnMonoFillScale       = True
  resov1@cnFillColor = "gray40"
  resov1cnFillDotSizeF = 0.005
  resov1@cnLevelSelectionMode  = "ExplicitLevels"
  resov1@cnFillPatterns        = (/8,-1/)
  resov1@cnFillScaleF = 1.2
  resov1@cnLevels              = (/.666/) ; .666
  resov1@cnInfoLabelOn         = False
  resov1@lbLabelBarOn        = False
  p = 0
  plot_stip1(p) =  gsn_csm_contour(wks,tas_aos,resov1)
  overlay(plot(p),plot_stip1(p))
  p = 1
  plot_stip1(p) =  gsn_csm_contour(wks,pr_aos,resov1)
  overlay(plot(p),plot_stip1(p))
  p = 2
  plot_stip1(p) =  gsn_csm_contour(wks,hdf_aos,resov1)
  overlay(plot(p),plot_stip1(p))


  ; -- mag
  ; -- overlay stippling
  resov2 = resov1
  resov2@cnFillPatterns        = (/17,-1/)
  p = 0
  plot_stip2(p) =  gsn_csm_contour(wks,tas_mag,resov2)
  overlay(plot(p),plot_stip2(p))
  p = 1
  plot_stip2(p) =  gsn_csm_contour(wks,pr_mag,resov2)
  overlay(plot(p),plot_stip2(p))
  p = 2
  plot_stip2(p) =  gsn_csm_contour(wks,hdf_mag,resov2)
  overlay(plot(p),plot_stip2(p))

  print_clock("shapefiles...")
  ; -- region mask contours (shapefiles of watersheds) -------------------------
  load "/home/flehner/scripts/ncl/shapefile_utils.ncl"
  shp_filename = "/project/ccr/flehner/shapefiles/HUC02_conus.shp"
  f = addfile(shp_filename, "r")   ; Open shapefile
  segments = f->segments
  geometry = f->geometry
  lon      = f->x
  lat      = f->y
  segsDims = dimsizes(segments)
  geomDims = dimsizes(geometry)
  ; -- Read global attributes
  geom_segIndex = f@geom_segIndex
  geom_numSegs  = f@geom_numSegs
  segs_xyzIndex = f@segs_xyzIndex
  segs_numPnts  = f@segs_numPnts
  lines2_0       = new(segsDims(0),graphic)   ; Array to hold polygons
  lines2_1       = new(segsDims(0),graphic)
  lines2_2       = new(segsDims(0),graphic)
  lines2_3       = new(segsDims(0),graphic)
  numFeatures = geomDims(0)
  segNum = 0
  plres             = True       ; resources for polylines
  plres@gsEdgesOn   = True       ; draw border around polygons
  plres@gsEdgeColor = "blue"
  plres@gsFillColor = "Transparent"
  plres@gsFillBackgroundColor = "Transparent"
  plres@gsEdgeThicknessF = 2
  plres@gsLineThicknessF = 1
  plres@tfPolyDrawOrder = "PostDraw"
  do i = 8,8;0,numFeatures-1
    startSegment = geometry(i, geom_segIndex)
    numSegments  = geometry(i, geom_numSegs)
    do seg=startSegment, startSegment+numSegments-1
      startPT = segments(seg, segs_xyzIndex)
      endPT = startPT + segments(seg, segs_numPnts) - 1
      lines2_0(segNum) = gsn_add_polygon(wks, plot(0), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines2_1(segNum) = gsn_add_polygon(wks, plot(1), lon(startPT:endPT), lat(startPT:endPT), plres)
      lines2_2(segNum) = gsn_add_polygon(wks, plot(2), lon(startPT:endPT), lat(startPT:endPT), plres)
      segNum = segNum + 1
    end do
  end do
  delete([/segments,geometry,lat,lon,numFeatures/])
  ; -- end of region mask contours (shapefiles of watersheds) ------------------

  ; -- add a box --
  ypts = (/ 36.0,  36.0,  42.0,   42.0, 36.0, 36.0/)
  xpts = (/-112.0, -106.0, -106.0, -112.0, -112.0, -106.0/)
  resp                  = True                      ; polyline mods desired
  resp@gsLineColor      = "black"                     ; color of lines
  resp@gsLineThicknessF = 1                      ; thickness of lines
  dum1 = new(5,graphic)
  dum2 = new(5,graphic)
  dum3 = new(5,graphic)
  ; -- draw each line separately. Each line must contain two points.
  do i = 0 , 4
    dum1(i) = gsn_add_polyline(wks,plot(0),xpts(i:i+1),ypts(i:i+1),resp)
    dum2(i) = gsn_add_polyline(wks,plot(1),xpts(i:i+1),ypts(i:i+1),resp)
    dum3(i) = gsn_add_polyline(wks,plot(2),xpts(i:i+1),ypts(i:i+1),resp)
  end do

  ; ---------------------------------------------------------------------

  resP1                  = True                    ; modify the panel plot
  resP1@gsnMaximize     = True
  resP1@gsnPanelRowSpec = True
  resP1@gsnPanelBottom    = 0.1
  resP1@gsnPanelTop    = 1.0

  labels =  (/"(a)","(b)","(c)"/)
  resP1@gsnPanelFigureStrings = labels
  resP1@gsnPanelFigureStringsFontHeightF = 0.015; 0.011

  resP1@txString        = "CMIP6 multi-model mean change ("+futstart+"-"+futende+")-("+refstart+"-"+refende+")"
  resP1@txFontHeightF = 0.02

  resP1@gsnPanelYWhiteSpacePercent = 1
  resP1@gsnPanelXWhiteSpacePercent = 1
  resP1@pmLabelBarOrthogonalPosF = -.04
  resP1@lbOrientation     = "horizontal"


  resP1@gsnPanelFigureStringsPerimOn = False
  resP1@amJust          = "BottomLeft"

  resP1@gsnPaperOrientation = "portrait"
  resP1@lbLabelStride = 1
  resP1@lbBoxLineColor = "gray20"
  resP1@gsnPanelLabelBar    = False;True       ; common label bar
  resP1@gsnFrame            = False      ; don't advance frame yet
  resP1@lbLabelAutoStride   = True       ; nice label bar labels

  gsn_panel(wks,plot,(/3,1/),resP1)
  delete(wks)

  system("convert "+varp+".eps "+varp+".png")


  print_clock("end script")

end
