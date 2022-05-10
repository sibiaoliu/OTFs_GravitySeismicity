#===== Script for RMBA of mid-ocean ridge-transform system=====
# 01.04.2022

#===== Input parameters' value change=====
    dataname=Blanco 
    #Clipperton, ChileRidge, Discovery, Marathon, Marion, Blanco, Oceanographer, Gofar, 
    #Garrett, Vlamingh, Atlantis, Clipperton, MarieCeleste, AtlantisII, Kane, DuTroit 
    Etas=(isov disl vp vep)
    gmt makecpt -T-23/15/1 -Z >grav_diff.cpt
    gmt makecpt -T-15/23/1 -Z >grav_diff2.cpt

# Define paths & constants
    globalDataPath=/Users/sliu/Downloads/SLiu/Gitlab-thermalStructure-OTFs/GlobalData/
    ETOPO1=${globalDataPath}/ETOPO1_Bed_g_gmt4.grd    
    faa=../inputs/${dataname}_faa.nc
    bathy_small=../inputs/${dataname}_small.nc
    bathy_ship=../inputs/${dataname}_ship.nc
    bathy_large_ship=../inputs/${dataname}_large_ship.nc
    bathy_large_Sat=../inputs/${dataname}_large_Sat.nc
    OTF=../inputs/${dataname}_OTF.txt
    ROTF=../inputs/${dataname}_ridge_otf.txt
    otffz_boxinfo=../inputs/otffz_boxinfo.txt
    mor_boxinfo=../inputs/mor_boxinfo.txt
    if [ ! -d Results ] ; then 
        mkdir Results
    fi    
    rho_crust=2700
    rho_water=1020
    rho_mantle=3300
    crust_thickness=6000

# 0. Calculate FAA & Mantle Bouger anomaly (MBA) using FFT approach in GMT
  # FAA -> Bouger gravity anomaly(BGA) -> Mantle Bouger anomaly(MBA)
  # Loading the ship-based bathymetry and satallite bathymetry (etopo1)
  # Outputs: Bathymetry & FAA for the individual ridge-transform system
  # Check the ranges of *ship.nc and large_ship.nc, if the ship.nc range 
  # is larger than large_ship.nc, then make the faa (from ship.nc) range smaller
  if [ ! -f $faa ]; then 
    echo "Warning, $faa doesn't exist, try to sample from global data"
    lon_min_small=`gmt grdinfo $bathy_large_ship | grep "x_min" | awk '{print $3+0.5}'`
    lon_max_small=`gmt grdinfo $bathy_large_ship | grep "x_max" | awk '{print $5-0.5}'`
    lat_min_small=`gmt grdinfo $bathy_large_ship | grep "y_min" | awk '{print $3+0.5}'`
    lat_max_small=`gmt grdinfo $bathy_large_ship | grep "y_max" | awk '{print $5-0.5}'`
    range_small=${lon_min_small}/${lon_max_small}/${lat_min_small}/${lat_max_small}
    gmt grdsample $bathy_ship -R$range_small -I100e -G${bathy_small}
    gmt grdsample $globalDataPath/grav.23.nc -R$bathy_small -I1m -G$faa
  fi
  # Merge bathymetry to a larger area
  if [ ! -f $bathy_large_ship ]; then 
    echo "Warning, $bathy_large_ship doesn't exist, try to merge from global bathymetry (ETOPO1)"
    dxdy=100e
    # merge ship-based bathymetry and satellite bathymetry from etopo1
    gmt grd2xyz $bathy_ship | gmt convert -bo > ship.b
    gmt nearneighbor -R$range_large -I$dxdy -S5k -Gship.nc ship.b -bi
    gmt grdsample $ETOPO1 -R$range_large -I$dxdy -Gtmp_etopo.nc
    gmt grdmath ship.nc tmp_etopo.nc AND = $bathy_large_ship 
    gmt grdsample $ETOPO1 -R$range_large -I$dxdy -G$bathy_large_Sat 
    rm ship.nc ship.b tmp_etopo.nc
  fi
  # Output the MBA
  mba=Results/${dataname}_mba.nc
  gmt gravfft $bathy_large_ship -D$(($rho_crust - $rho_water)) -W$W -Ff -fg -Gg_wc.grd -E4
  gmt gravfft ${bathy_large_ship}=+o-${crust_thickness} -D$(($rho_mantle - $rho_crust)) -fg -Gg_cm.grd -E4
  # make sure the grid size same as faa
  gmt grdsample g_wc.grd -R$faa -Gg_wc.grd -V0
  gmt grdsample g_cm.grd -R$faa -Gg_cm.grd -V0
  gmt grdmath $faa g_wc.grd SUB g_cm.grd SUB = $mba
  rm g_wc.grd g_cm.grd

# Required functions:
function gravity_xyz2grd()
{
    gravityFile=$1
    gmt grd2xyz ${faa} >${dataname}_faa.xyz
    paste ${dataname}_faa.xyz ${gravityFile}.txt  | awk '{print $1,$2,$4}' >${gravityFile}.lonlat
    gmt xyz2grd ${gravityFile}.lonlat -R$faa -G${gravityFile}.nc
    rm ${gravityFile}.lonlat ${dataname}_faa.xyz
}
function plotControlTransformFault()
{
    file_ROTF=$1
    file_OTF=$2
    if [ -z $file_ROTF ]; then 
        file_ROTF=$ROTF
    fi
    if [ -z $file_OTF ]; then 
        file_OTF=$OTF
    fi
    if [ -z $RTI_ms ]; then 
        RTI_ms=0.2
    fi

    if [ -f ${file_ROTF} ]; then 
        gmt plot ${file_ROTF} -W1p,white
        gmt plot ${file_ROTF} -W0.5p,black,-
        # gmt plot ${file_ROTF} -Sd0.1c -Gwhite -W0.5p,black
    else
        echo "please create ${dataname}_ridge_otf.txt"
    fi
    if [ ! -z $file_OTF ]; then 
        if [ -f ${file_OTF} ]; then 
            numOTF=`awk 'END{print NR}' ${file_OTF}`
            for ((i=1; i<=$numOTF; i=i+2))
            do 
                # awk -v row0=$i 'NR>=row0 && NR<=(row0+1){print }' ${file_OTF} |\
                # gmt plot -W1p,white
                awk -v row0=$i 'NR>=row0 && NR<=(row0+1){print }' ${file_OTF} |\
                gmt plot -W0.5p,black
            done
            # gmt plot ${file_OTF} -Sk@bullseye -Gblack
            gmt plot ${file_ROTF} -Sd${RTI_ms}c -Gwhite -W0.5p,black
        else
            echo "please create ${file_OTF}"
        fi
    fi
    
}
function shiftDataTo0OTF()
{
    grd_data=$1
    start_otf=`awk 'NR==1{print $1"/"$2}' ${OTF}`
    stop_otf=`awk 'END{print $1"/"$2}' ${OTF}`
    gmt grdtrack -G$grd_data -E${start_otf}/${stop_otf}+i1k+d >data_OTF.txt
    meanValue_OTF=`awk -v sum_data=0 -v num_data=0 '{ while ( getline == 1 ) { sum_data+=$4; num_data+=1; }printf "%.0f", sum_data/num_data} ' data_OTF.txt`
    gmt grdmath ${grd_data} ${meanValue_OTF} SUB = ${grd_data}
    str_shiftData=", shift ${meanValue_OTF} mGal"
}
function makecpt_grd()
{
    grdfile=$1
    masterCPT_grav=../../../basecpt_grav.cpt
    data_min=`gmt grdinfo $grdfile | grep "v_min" | awk '{printf "%.1f", $3}'`
    data_max=`gmt grdinfo $grdfile | grep "v_max" | awk '{printf "%.1f", $5}'`
    cpt_min=`echo ${data_min} | awk '{printf "%.1f", $1}'`
    cpt_max=`echo ${data_max} | awk '{printf "%.1f", $1}'`
    if (( $(echo "$data_min > 0" |bc -l) )); then
        gmt makecpt -C${masterCPT_grav} -T${data_min}/${data_max}
    else
        gmt makecpt -C${masterCPT_grav}+h -T${cpt_min}/${cpt_max}
    fi
}
function makecpt_grd_basecpt()
{
    grdfile=$1
    basecpt=../../../basecpt_thermal.cpt
    data_min=`gmt grdinfo $grdfile | grep "v_max" | awk '{printf "%.1f", $3}'`
    data_max=`gmt grdinfo $grdfile | grep "v_max" | awk '{printf "%.1f", $5}'`
    cpt_min=`echo ${data_min} | awk '{printf "%.1f", $1*1.5}'`
    cpt_max=`echo ${data_max} | awk '{printf "%.1f", $1*1.5}'`
    if (( $(echo "$data_min > 0" |bc -l) )); then
        gmt makecpt -C${basecpt} -T${data_min}/${data_max}
    else
        gmt makecpt -C${basecpt}+h -T${cpt_min}/${cpt_max}
    fi
    
}
function getAverageBox_TF_FZ()
{
    boxname=$1
    ind_OTF=$2
    loc_frac_box=$3 #-1.2 #-1 means left FZ; +1 means right FZ; 0 means OTF
    length_frac_box=$4 # 0.8 # percentage of TF length
    w_box=$5
    halfWidth_box=`echo $w_box | awk '{print $1/2*1000}'` #m
    additional_roation=$6
    lon_start_otf=`awk -v indOTF=$ind_OTF 'NR==(indOTF*2-1){print $1}' $OTF`
    lon_end_otf=`awk -v indOTF=$ind_OTF 'NR==(indOTF*2){print $1}' $OTF`
    lat_start_otf=`awk -v indOTF=$ind_OTF 'NR==(indOTF*2-1){print $2}' $OTF`
    lat_end_otf=`awk -v indOTF=$ind_OTF 'NR==(indOTF*2){print $2}' $OTF`
    # start point of OTF assumes to (0,0)
    # OTF/FZ center point: lon/lat & x/y 
    lonC_otf=`echo $lon_start_otf $lon_end_otf | awk '{print ($1+$2)/2.0}'`
    latC_otf=`echo $lat_start_otf $lat_end_otf | awk '{print ($1+$2)/2.0}'`
    xCyC_otf=`echo $lonC_otf $latC_otf | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 | awk '{print $1, $2}'`
    XmaxYmax_otf=`echo $lon_end_otf $lat_end_otf | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 | awk '{print $1, $2}'`
    length_otf=`echo $XmaxYmax_otf | awk '{print sqrt($1*$1+$2*$2)}'`
    length_box=`echo $length_otf $length_frac_box | awk '{print $1*$2}'`
    costheta=`echo $XmaxYmax_otf | awk '{printf "%.8f", $1/sqrt($1*$1+$2*$2)}'`
    sintheta=`echo $XmaxYmax_otf | awk '{printf "%.8f",$2/sqrt($1*$1+$2*$2)}'`
    angle_rot=`echo $sintheta $costheta | awk '{printf "%.1f", atan2($1, $2)/3.141592653*180}'`
    # set final center xy: local coordinate on OTF/FZ box is [-0.5, 0.5]
    xCyC_box=`echo $xCyC_otf $costheta $sintheta $length_otf $loc_frac_box $length_box | awk '{print $1+$6*$5*$3, $2+$5*$6*$4}'`
    lonClatC_box=`echo $xCyC_box | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}'`

    # xy coordinate of OTF/FZ box start/end point 
    xy_start_box=`echo $xCyC_box $costheta $sintheta $length_box | awk '{print $1-$5/2.0*$3, $2-$5/2.0*$4}'`
    xy_end_box=`echo $xCyC_box $costheta $sintheta $length_box | awk '{print $1+$5/2.0*$3, $2+$5/2.0*$4}'`
    
    # xy coordinate of averate box nodes
    lonlat_start_box=`echo $xy_start_box | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}'`
    echo $xy_start_box $costheta $sintheta $halfWidth_box | awk '{print $1-$5*$4, $2+$5*$3}' >../inputs/${boxname}.xy
    echo $xy_start_box $costheta $sintheta $halfWidth_box | awk '{print $1+$5*$4, $2-$5*$3}' >>../inputs/${boxname}.xy
    echo $xy_end_box $costheta $sintheta $halfWidth_box | awk '{print $1+$5*$4, $2-$5*$3}' >>../inputs/${boxname}.xy
    echo $xy_end_box $costheta $sintheta $halfWidth_box | awk '{print $1-$5*$4, $2+$5*$3}' >>../inputs/${boxname}.xy
    # additional rotation with rotate center one of RTIs: if loc_frac_box<0, the first RTI; if loc_frac_box>0 the second RTI
    row_RTI=`echo $ind_OTF | awk '{print $1*2}'`
    if [ `echo "$loc_frac_box < 0" | bc` -eq 1 ]; then 
        row_RTI=`echo $ind_OTF | awk '{print $1*2-1}'`
    fi 
    x0y0=(`awk -v row=$row_RTI 'NR==row{print $1, $2}' $OTF | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 | awk '{print $1, $2}'`)
    costheta=`echo $additional_roation | awk '{printf "%.6f", cos(($1)/180*3.141592653)}'`
    sintheta=`echo $additional_roation | awk '{printf "%.6f", sin(($1)/180*3.141592653)}'`
    awk -v co=$costheta -v si=$sintheta -v x0=${x0y0[0]} -v y0=${x0y0[1]} '{print co*($1-x0)+si*($2-y0) +x0, -si*($1-x0)+co*($2-y0) + y0, $3}' ../inputs/${boxname}.xy >tmp.xy
    mv tmp.xy ../inputs/${boxname}.xy
    # convert box coordinate (x,y) to (lon,lat)
    cat ../inputs/${boxname}.xy | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}' >../inputs/${boxname}.lonlat
    # lonClatC_box
    lonClatC_box=`awk '{{meanlon+=$1; meanlat+=$2}; if(NR==4){print meanlon/4, meanlat/4}}' ../inputs/${boxname}.lonlat`
    angle_rot=`echo "$angle_rot - $additional_roation" | bc `
    #echo $boxname $lonClatC_box $angle_rot
    echo $lonClatC_box $angle_rot | awk '{print $1, $2, $3}' >../inputs/BoxC_${boxname}.lonlat
}
function getAverageBox_Ridge()
{
    boxname=$1
    ind_Ridge=$2
    loc_frac_box=$3 #-1.2 #-1 means left FZ; +1 means right FZ; 0 means OTF
    length_frac_box=$4 # 0.8 # percentage of MOR length
    w_box=$5
    halfWidth_box=`echo $w_box | awk '{print $1/2*1000}'` #m
    additional_roation=$6
    file_ridge_index=../inputs/index_ridge.txt
    file_ridge_otf=$ROTF

    index_start_ridge=`awk -v indRidge=$ind_Ridge 'NR==indRidge{print $1}' $file_ridge_index`
    index_end_ridge=`awk -v indRidge=$ind_Ridge 'NR==indRidge{print $2}' $file_ridge_index`

    lon_start_ridge=`awk -v row=$index_start_ridge 'NR==row{print $1}' $file_ridge_otf`
    lon_end_ridge=`awk -v row=$index_end_ridge 'NR==row{print $1}' $file_ridge_otf`
    lat_start_ridge=`awk -v row=$index_start_ridge 'NR==row{print $2}' $file_ridge_otf`
    lat_end_ridge=`awk -v row=$index_end_ridge 'NR==row{print $2}' $file_ridge_otf`

    # start point of Ridge with coordinate of (0,0)
    # center point of Ridge
    lonC_ridge=`echo $lon_start_ridge $lon_end_ridge | awk '{print ($1+$2)/2.0}'`
    latC_ridge=`echo $lat_start_ridge $lat_end_ridge | awk '{print ($1+$2)/2.0}'`
    xCyC_ridge=`echo $lonC_ridge $latC_ridge | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_ridge +lat_0=$lat_start_ridge +x_0=0 +y_0=0 | awk '{print $1, $2}'`
    XmaxYmax_ridge=`echo $lon_end_ridge $lat_end_ridge | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_ridge +lat_0=$lat_start_ridge +x_0=0 +y_0=0 | awk '{print $1, $2}'`
    length_ridge=`echo $XmaxYmax_ridge | awk '{print sqrt($1*$1+$2*$2)}'`
    length_box=`echo $length_ridge $length_frac_box | awk '{print $1*$2}'`
    costheta=`echo $XmaxYmax_ridge | awk '{printf "%.8f", $1/sqrt($1*$1+$2*$2)}'`
    sintheta=`echo $XmaxYmax_ridge | awk '{printf "%.8f",$2/sqrt($1*$1+$2*$2)}'`
    angle_rot=`echo $sintheta $costheta | awk '{printf "%.1f", atan2($1, $2)/3.141592653*180}'`
    # set final center xy: local coordinate on Ridge is [-0.5, 0.5]
    xCyC_box=`echo $xCyC_ridge $costheta $sintheta $length_ridge $loc_frac_box $length_box | awk '{print $1+$6*$5*$3, $2+$5*$6*$4}'`
    lonClatC_box=`echo $xCyC_box | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_ridge +lat_0=$lat_start_ridge +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}'`
    # xy coordinate of Ridge start point 
    xy_start_box=`echo $xCyC_box $costheta $sintheta $length_box | awk '{print $1-$5/2.0*$3, $2-$5/2.0*$4}'`
    # xy coordinate of Ridge end point 
    xy_end_box=`echo $xCyC_box $costheta $sintheta $length_box | awk '{print $1+$5/2.0*$3, $2+$5/2.0*$4}'`
    
    # xy coordinate of averate box nodes
    lonlat_start_box=`echo $xy_start_box | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_ridge +lat_0=$lat_start_ridge +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}'`
    echo $xy_start_box $costheta $sintheta $halfWidth_box | awk '{print $1-$5*$4, $2+$5*$3}' >../inputs/${boxname}.xy
    echo $xy_start_box $costheta $sintheta $halfWidth_box | awk '{print $1+$5*$4, $2-$5*$3}' >>../inputs/${boxname}.xy
    echo $xy_end_box $costheta $sintheta $halfWidth_box | awk '{print $1+$5*$4, $2-$5*$3}' >>../inputs/${boxname}.xy
    echo $xy_end_box $costheta $sintheta $halfWidth_box | awk '{print $1-$5*$4, $2+$5*$3}' >>../inputs/${boxname}.xy

    cat ../inputs/${boxname}.xy | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_ridge +lat_0=$lat_start_ridge +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}' >../inputs/${boxname}.lonlat
    # lonClatC_box
    lonClatC_box=`awk '{{meanlon+=$1; meanlat+=$2}; if(NR==4){print meanlon/4, meanlat/4}}' ../inputs/${boxname}.lonlat`
    angle_rot=`echo "$angle_rot - $additional_roation" | bc `
    echo $lonClatC_box $angle_rot | awk '{print $1, $2, $3}' >../inputs/BoxC_${boxname}.lonlat
}
# 2. Generate the meanRMBA boxes for OTF, FZ, MOR, IC, and OC
  # Call functions: (1) getAverageBox_TF_FZ; (2)getAverageBox_Ridge;
  # (3) getAverageBox_ICOC
  # 2.1 OTF-FZ boxes
  nbox_otffz=`awk 'END{print NR}' $otffz_boxinfo`
  for (( j=2; j<=$nbox_otffz; j++ ));
    do 
    boxname=averageBox_`awk -v row=$j 'NR==row{print $1}' $otffz_boxinfo`
    ind_OTF=`awk -v row=$j 'NR==row{print $2}' $otffz_boxinfo`
    loc_box=`awk -v row=$j 'NR==row{print $3}' $otffz_boxinfo`
    l_box=`awk -v row=$j 'NR==row{print $4}' $otffz_boxinfo`
    w_box=`awk -v row=$j 'NR==row{print $5}' $otffz_boxinfo`
    rot_angle=`awk -v row=$j 'NR==row{print $6}' $otffz_boxinfo`
    boxcolor=`awk -v row=$j 'NR==row{print $7}' $otffz_boxinfo`
    getAverageBox_TF_FZ $boxname $ind_OTF $loc_box $l_box $w_box $rot_angle
  done
  # 2.2 MOR boxes
  nbox_mor=`awk 'END{print NR}' $mor_boxinfo`
  for (( k=2; k<=$nbox_mor; k++ ));
    do 
    boxname=averageBox_`awk -v row=$k 'NR==row{print $1}' $mor_boxinfo`
    ind_mor=`awk -v row=$k 'NR==row{print $2}' $mor_boxinfo`
    loc_box=`awk -v row=$k 'NR==row{print $3}' $mor_boxinfo`
    l_box=`awk -v row=$k 'NR==row{print $4}' $mor_boxinfo`
    w_box=`awk -v row=$k 'NR==row{print $5}' $mor_boxinfo`
    rot_angle=`awk -v row=$k 'NR==row{print $6}' $mor_boxinfo`
    boxcolor=`awk -v row=$k 'NR==row{print $7}' $mor_boxinfo`
    getAverageBox_Ridge $boxname $ind_mor $loc_box $l_box $w_box $rot_angle
  done  
  # 2.3 IC & OC boxes

# 3. Calculate RMBA field with rheologic-dependent thermal model
  # Call the functions: (1) gravity_xyz2grd(); (2) shiftDataTo0OTF;
  # (3) plotControlTransformFault; (4) makecpt_grd;
  # (5) makecpt_grd_basecpt  
  # Figure includes 1st row (4): gravity of rheologic-dependent thermal model
  # 2nd row (4): RMBA for rheoglical differing models
  # 3rd row (4): RMBA difference betwen them: a. disl-isov
  # b. vp-isov, c. vp-disl, d. vep-vp
    model1=Isoviscous
    model2=Viscous
    model3=Visco-plastic
    model4=Visco-elasto-plastic
    #afg_grav=a20f5  #Using in RMBA  
    gmt begin Results/${dataname}_grav_models pdf
        #gmt set PS_MEDIA A3 #if you use other canvas sizes for the PS format
        gmt subplot begin 3x4 -A+JTL+o0.5c -Fs15c/9c -M0.5c/1.5c -R$faa -JM15c -Ba1df30m -BWSne -T"Gravity anomaly with rheologic-dependent thermal model: "$dataname
            gmt subplot set 0,0
                modelname=$model1
                etaname=${Etas[0]}
                #Thermal gravity anomaly
                gravity_xyz2grd grav_${etaname}
                grav_therm=grav_${etaname}.nc
                grav_therm_minusMeanValue=grav_${etaname}_minusMean.nc
                # mean value of the thermal gravity
                meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
                gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
                makecpt_grd_basecpt $grav_therm_minusMeanValue
                gmt grdimage $grav_therm_minusMeanValue -BwsEN #-I${grav_therm_minusMeanValue}.grad
                gmt grdcontour $grav_therm_minusMeanValue -C10
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal"
                plotControlTransformFault
            gmt subplot set 0,1
                modelname=$model2
                etaname=${Etas[1]}
                #Thermal gravity anomaly
                gravity_xyz2grd grav_${etaname}
                grav_therm=grav_${etaname}.nc
                grav_therm_minusMeanValue=grav_${etaname}_minusMean.nc
                # mean value of the thermal gravity
                meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
                gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
                makecpt_grd_basecpt $grav_therm_minusMeanValue
                gmt grdimage $grav_therm_minusMeanValue -BwseN #-I${grav_therm_minusMeanValue}.grad
                gmt grdcontour $grav_therm_minusMeanValue -C10
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal"
                plotControlTransformFault
            gmt subplot set 0,2
                modelname=$model3
                etaname=${Etas[2]}
                #Thermal gravity anomaly
                gravity_xyz2grd grav_${etaname}
                grav_therm=grav_${etaname}.nc
                grav_therm_minusMeanValue=grav_${etaname}_minusMean.nc
                # mean value of the thermal gravity
                meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
                gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
                makecpt_grd_basecpt $grav_therm_minusMeanValue
                gmt grdimage $grav_therm_minusMeanValue -BwsEN #-I${grav_therm_minusMeanValue}.grad
                gmt grdcontour $grav_therm_minusMeanValue -C10
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal"
                plotControlTransformFault 
            gmt subplot set 0,3
                modelname=$model4
                etaname=${Etas[3]}
                #Thermal gravity anomaly
                gravity_xyz2grd grav_${etaname}
                grav_therm=grav_${etaname}.nc
                grav_therm_minusMeanValue=grav_${etaname}_minusMean.nc
                # mean value of the thermal gravity
                meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
                gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
                makecpt_grd_basecpt $grav_therm_minusMeanValue
                gmt grdimage $grav_therm_minusMeanValue -BwsEN #-I${grav_therm_minusMeanValue}.grad
                gmt grdcontour $grav_therm_minusMeanValue -C10
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal"
                plotControlTransformFault
            # 2nd row --- RMBA maps
            gmt subplot set 1,0
                # First, get RMBA
                modelname=$model1
                etaname=${Etas[0]}
                rmba=Results/${dataname}_rmba_${etaname}.nc
                gmt grdmath ${mba} grav_${etaname}.nc SUB = ${rmba}
                # If required, make a shif of RMBA along OTF close to zero 
                shiftDataTo0OTF $rmba
                gmt grdgradient ${rmba} -A30 -Nt0.6 -Qc -G${rmba}.grad
                # Second, plot RMBA
                gmt basemap -BwseN -Ba 
                makecpt_grd $rmba
                gmt grdimage $rmba -I${rmba}.grad
                gmt grdcontour $rmba -C10
                gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                #if [ -f $grid_mask ]; then 
                #    gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                #fi 
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}                
                # Third, plot the mean RMBA in OTF-FZ boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >Results/averageRMBA_TFZ_${etaname}.txt
                nbox_otffz=`awk 'END{print NR}' $otffz_boxinfo`
                for (( ii=2; ii<=$nbox_otffz; ii++ ));
                    do 
                    boxname=averageBox_`awk -v row=$ii 'NR==row{print $1}' $otffz_boxinfo`
                    boxcolor=`awk -v row=$ii 'NR==row{print $7}' $otffz_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $rmba -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/averageRMBA_TFZ_${etaname}.txt
                    # plot text of mean rmba
                    lonClatC_box_otffz=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_otffz=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`
                    echo $lonClatC_box_otffz $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otffz
                done
                # Fourth, plot the mean RMBA in MOR boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >Results/averageRMBA_MOR_${etaname}.txt
                nbox_mor=`awk 'END{print NR}' $mor_boxinfo`
                for (( jj=2; jj<=$nbox_mor; jj++ ));
                    do 
                    boxname=averageBox_`awk -v row=$jj 'NR==row{print $1}' $mor_boxinfo`
                    boxcolor=`awk -v row=$jj 'NR==row{print $7}' $mor_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $rmba -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/averageRMBA_MOR_${etaname}.txt
                    # plot text of mean rmba
                    lonClatC_box_mor=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_mor=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`                    
                    echo $lonClatC_box_mor $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor
                done
                plotControlTransformFault
                unset str_shiftData
            gmt subplot set 1,1
                # First, get RMBA
                modelname=$model2
                etaname=${Etas[1]}
                rmba=Results/${dataname}_rmba_${etaname}.nc
                gmt grdmath ${mba} grav_${etaname}.nc SUB = ${rmba}
                # If required, make a shif of RMBA along OTF close to zero 
                shiftDataTo0OTF $rmba
                gmt grdgradient ${rmba} -A30 -Nt0.6 -Qc -G${rmba}.grad
                # Second, plot RMBA
                gmt basemap -BwseN -Ba 
                makecpt_grd $rmba
                gmt grdimage $rmba -I${rmba}.grad
                gmt grdcontour $rmba -C10
                gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                #if [ -f $grid_mask ]; then 
                #    gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                #fi 
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}                
                # Third, plot the mean RMBA in OTF-FZ boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >Results/averageRMBA_TFZ_${etaname}.txt
                nbox_otffz=`awk 'END{print NR}' $otffz_boxinfo`
                for (( ii=2; ii<=$nbox_otffz; ii++ ));
                    do 
                    boxname=averageBox_`awk -v row=$ii 'NR==row{print $1}' $otffz_boxinfo`
                    boxcolor=`awk -v row=$ii 'NR==row{print $7}' $otffz_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $rmba -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/averageRMBA_TFZ_${etaname}.txt
                    # plot text of mean rmba
                    lonClatC_box_otffz=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_otffz=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`
                    echo $lonClatC_box_otffz $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otffz
                done
                # Fourth, plot the mean RMBA in MOR boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >Results/averageRMBA_MOR_${etaname}.txt
                nbox_mor=`awk 'END{print NR}' $mor_boxinfo`
                for (( jj=2; jj<=$nbox_mor; jj++ ));
                    do 
                    boxname=averageBox_`awk -v row=$jj 'NR==row{print $1}' $mor_boxinfo`
                    boxcolor=`awk -v row=$jj 'NR==row{print $7}' $mor_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $rmba -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/averageRMBA_MOR_${etaname}.txt
                    # plot text of mean rmba
                    lonClatC_box_mor=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_mor=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`                    
                    echo $lonClatC_box_mor $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor
                done
                plotControlTransformFault
                unset str_shiftData                
            gmt subplot set 1,2
                # First, get RMBA
                modelname=$model3
                etaname=${Etas[2]}
                rmba=Results/${dataname}_rmba_${etaname}.nc
                gmt grdmath ${mba} grav_${etaname}.nc SUB = ${rmba}
                # If required, make a shif of RMBA along OTF close to zero 
                shiftDataTo0OTF $rmba
                gmt grdgradient ${rmba} -A30 -Nt0.6 -Qc -G${rmba}.grad
                # Second, plot RMBA
                gmt basemap -BwseN -Ba 
                makecpt_grd $rmba
                gmt grdimage $rmba -I${rmba}.grad
                gmt grdcontour $rmba -C10
                gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                #if [ -f $grid_mask ]; then 
                #    gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                #fi 
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}                
                # Third, plot the mean RMBA in OTF-FZ boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >Results/averageRMBA_TFZ_${etaname}.txt
                nbox_otffz=`awk 'END{print NR}' $otffz_boxinfo`
                for (( ii=2; ii<=$nbox_otffz; ii++ ));
                    do 
                    boxname=averageBox_`awk -v row=$ii 'NR==row{print $1}' $otffz_boxinfo`
                    boxcolor=`awk -v row=$ii 'NR==row{print $7}' $otffz_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $rmba -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/averageRMBA_TFZ_${etaname}.txt
                    # plot text of mean rmba
                    lonClatC_box_otffz=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_otffz=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`
                    echo $lonClatC_box_otffz $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otffz
                done
                # Fourth, plot the mean RMBA in MOR boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >Results/averageRMBA_MOR_${etaname}.txt
                nbox_mor=`awk 'END{print NR}' $mor_boxinfo`
                for (( jj=2; jj<=$nbox_mor; jj++ ));
                    do 
                    boxname=averageBox_`awk -v row=$jj 'NR==row{print $1}' $mor_boxinfo`
                    boxcolor=`awk -v row=$jj 'NR==row{print $7}' $mor_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $rmba -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/averageRMBA_MOR_${etaname}.txt
                    # plot text of mean rmba
                    lonClatC_box_mor=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_mor=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`                    
                    echo $lonClatC_box_mor $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor
                done
                plotControlTransformFault
                rm ${rmba}.grad
                unset str_shiftData
            gmt subplot set 1,3
                # First, get RMBA
                modelname=$model4
                etaname=${Etas[3]}
                rmba=Results/${dataname}_rmba_${etaname}.nc
                gmt grdmath ${mba} grav_${etaname}.nc SUB = ${rmba}
                # If required, make a shif of RMBA along OTF close to zero 
                shiftDataTo0OTF $rmba
                gmt grdgradient ${rmba} -A30 -Nt0.6 -Qc -G${rmba}.grad
                # Second, plot RMBA
                gmt basemap -BwseN -Ba 
                makecpt_grd $rmba
                gmt grdimage $rmba -I${rmba}.grad
                gmt grdcontour $rmba -C10
                gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                #if [ -f $grid_mask ]; then 
                #    gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                #fi 
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}                
                # Third, plot the mean RMBA in OTF-FZ boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >Results/averageRMBA_TFZ_${etaname}.txt
                nbox_otffz=`awk 'END{print NR}' $otffz_boxinfo`
                for (( ii=2; ii<=$nbox_otffz; ii++ ));
                    do 
                    boxname=averageBox_`awk -v row=$ii 'NR==row{print $1}' $otffz_boxinfo`
                    boxcolor=`awk -v row=$ii 'NR==row{print $7}' $otffz_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $rmba -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/averageRMBA_TFZ_${etaname}.txt
                    # plot text of mean rmba
                    lonClatC_box_otffz=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_otffz=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`
                    echo $lonClatC_box_otffz $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otffz
                done
                # Fourth, plot the mean RMBA in MOR boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >Results/averageRMBA_MOR_${etaname}.txt
                nbox_mor=`awk 'END{print NR}' $mor_boxinfo`
                for (( jj=2; jj<=$nbox_mor; jj++ ));
                    do 
                    boxname=averageBox_`awk -v row=$jj 'NR==row{print $1}' $mor_boxinfo`
                    boxcolor=`awk -v row=$jj 'NR==row{print $7}' $mor_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $rmba -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/averageRMBA_MOR_${etaname}.txt
                    # plot text of mean rmba
                    lonClatC_box_mor=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_mor=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`                    
                    echo $lonClatC_box_mor $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor
                done
                plotControlTransformFault
                unset str_shiftData                
            # 3rd row --- RMBA difference between different rheologies   
            gmt subplot set 2,0
                m1=$model1
                m2=$model2
                grav_therm1=grav_${Etas[0]}.nc
                grav_therm2=grav_${Etas[1]}.nc
                # RMBA_diff_disl-isov = grav_diff_isov-disl!
                grav_diff=rmba_diff_${Etas[1]}-${Etas[0]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                #shiftDataTo0OTF $grav_diff #not necessary
                gmt grdimage $grav_diff -BWseN -Cgrav_diff2.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                #gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity difference: $m2 - $m1,$str_shiftData" -By+l"mGal" -Cgrav_diff2.cpt -G${data_min}/${data_max}
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference (no shift): $m2 - $m1" -By+l"mGal" -Cgrav_diff2.cpt 
                plotControlTransformFault
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
                # plot the mean RMBA in OTF-FZ boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2 in the RMBA difference map" >Results/meanRMBA_TFZ_${Etas[1]}-${Etas[0]}.txt
                nbox_otffz=`awk 'END{print NR}' $otffz_boxinfo`
                for (( ii=2; ii<=$nbox_otffz; ii++ ));
                    do 
                    boxname=averageBox_`awk -v row=$ii 'NR==row{print $1}' $otffz_boxinfo`
                    boxcolor=`awk -v row=$ii 'NR==row{print $7}' $otffz_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $grav_diff -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/meanRMBA_TFZ_${Etas[1]}-${Etas[0]}.txt
                    # plot text of mean rmba
                    lonClatC_box_otffz=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_otffz=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`
                    echo $lonClatC_box_otffz $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otffz
                done
                # Fourth, plot the mean RMBA in MOR boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2 in the RMBA difference map" >Results/meanRMBA_MOR_${Etas[1]}-${Etas[0]}.txt
                nbox_mor=`awk 'END{print NR}' $mor_boxinfo`
                for (( jj=2; jj<=$nbox_mor; jj++ ));
                    do 
                    boxname=averageBox_`awk -v row=$jj 'NR==row{print $1}' $mor_boxinfo`
                    boxcolor=`awk -v row=$jj 'NR==row{print $7}' $mor_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $grav_diff -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/meanRMBA_MOR_${Etas[1]}-${Etas[0]}.txt
                    # plot text of mean rmba
                    lonClatC_box_mor=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_mor=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`                    
                    echo $lonClatC_box_mor $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor
                done                
            gmt subplot set 2,1
                m1=$model1
                m2=$model3
                grav_therm1=grav_${Etas[0]}.nc
                grav_therm2=grav_${Etas[2]}.nc
                grav_diff=rmba_diff_${Etas[2]}-${Etas[0]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BwseN -Cgrav_diff2.cpt #-I${grav_diff}.grad  
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference (no shift): $m2 - $m1" -By+l"mGal"  -Cgrav_diff2.cpt
                plotControlTransformFault
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
                # plot the mean RMBA in OTF-FZ boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2 in the RMBA difference map" >Results/meanRMBA_TFZ_${Etas[2]}-${Etas[0]}.txt
                nbox_otffz=`awk 'END{print NR}' $otffz_boxinfo`
                for (( ii=2; ii<=$nbox_otffz; ii++ ));
                    do 
                    boxname=averageBox_`awk -v row=$ii 'NR==row{print $1}' $otffz_boxinfo`
                    boxcolor=`awk -v row=$ii 'NR==row{print $7}' $otffz_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $grav_diff -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/meanRMBA_TFZ_${Etas[2]}-${Etas[0]}.txt
                    # plot text of mean rmba
                    lonClatC_box_otffz=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_otffz=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`
                    echo $lonClatC_box_otffz $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otffz
                done
                # Fourth, plot the mean RMBA in MOR boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2 in the RMBA difference map" >Results/meanRMBA_MOR_${Etas[2]}-${Etas[0]}.txt
                nbox_mor=`awk 'END{print NR}' $mor_boxinfo`
                for (( jj=2; jj<=$nbox_mor; jj++ ));
                    do 
                    boxname=averageBox_`awk -v row=$jj 'NR==row{print $1}' $mor_boxinfo`
                    boxcolor=`awk -v row=$jj 'NR==row{print $7}' $mor_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $grav_diff -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/meanRMBA_MOR_${Etas[2]}-${Etas[0]}.txt
                    # plot text of mean rmba
                    lonClatC_box_mor=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_mor=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`                    
                    echo $lonClatC_box_mor $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor
                done                
            gmt subplot set 2,2
                m1=$model2
                m2=$model3
                grav_therm1=grav_${Etas[1]}.nc
                grav_therm2=grav_${Etas[2]}.nc
                grav_diff=rmba_diff_${Etas[2]}-${Etas[1]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff2.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference (no shift): $m2 - $m1" -By+l"mGal" -Cgrav_diff2.cpt
                plotControlTransformFault
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
                # plot the mean RMBA in OTF-FZ boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2 in the RMBA difference map" >Results/meanRMBA_TFZ_${Etas[2]}-${Etas[1]}.txt
                nbox_otffz=`awk 'END{print NR}' $otffz_boxinfo`
                for (( ii=2; ii<=$nbox_otffz; ii++ ));
                    do 
                    boxname=averageBox_`awk -v row=$ii 'NR==row{print $1}' $otffz_boxinfo`
                    boxcolor=`awk -v row=$ii 'NR==row{print $7}' $otffz_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $grav_diff -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/meanRMBA_TFZ_${Etas[2]}-${Etas[1]}.txt
                    # plot text of mean rmba
                    lonClatC_box_otffz=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_otffz=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`
                    echo $lonClatC_box_otffz $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otffz
                done
                # Fourth, plot the mean RMBA in MOR boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2 in the RMBA difference map" >Results/meanRMBA_MOR_${Etas[2]}-${Etas[1]}.txt
                nbox_mor=`awk 'END{print NR}' $mor_boxinfo`
                for (( jj=2; jj<=$nbox_mor; jj++ ));
                    do 
                    boxname=averageBox_`awk -v row=$jj 'NR==row{print $1}' $mor_boxinfo`
                    boxcolor=`awk -v row=$jj 'NR==row{print $7}' $mor_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $grav_diff -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/meanRMBA_MOR_${Etas[2]}-${Etas[1]}.txt
                    # plot text of mean rmba
                    lonClatC_box_mor=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_mor=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`                    
                    echo $lonClatC_box_mor $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor
                done
            gmt subplot set 2,3
                m1=$model3
                m2=$model4
                grav_therm1=grav_${Etas[2]}.nc
                grav_therm2=grav_${Etas[3]}.nc
                grav_diff=rmba_diff_${Etas[3]}-${Etas[2]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BwseN -Cgrav_diff2.cpt #-I${grav_diff}.grad  
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference (no shift): $m2 - $m1" -By+l"mGal" -Cgrav_diff2.cpt
                plotControlTransformFault
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
                # plot the mean RMBA in OTF-FZ boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2 in the RMBA difference map" >Results/meanRMBA_TFZ_${Etas[3]}-${Etas[2]}.txt
                nbox_otffz=`awk 'END{print NR}' $otffz_boxinfo`
                for (( ii=2; ii<=$nbox_otffz; ii++ ));
                    do 
                    boxname=averageBox_`awk -v row=$ii 'NR==row{print $1}' $otffz_boxinfo`
                    boxcolor=`awk -v row=$ii 'NR==row{print $7}' $otffz_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $grav_diff -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/meanRMBA_TFZ_${Etas[3]}-${Etas[2]}.txt
                    # plot text of mean rmba
                    lonClatC_box_otffz=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_otffz=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`
                    echo $lonClatC_box_otffz $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otffz
                done
                # Fourth, plot the mean RMBA in MOR boxes and save the data
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2 in the RMBA difference map" >Results/meanRMBA_MOR_${Etas[3]}-${Etas[2]}.txt
                nbox_mor=`awk 'END{print NR}' $mor_boxinfo`
                for (( jj=2; jj<=$nbox_mor; jj++ ));
                    do 
                    boxname=averageBox_`awk -v row=$jj 'NR==row{print $1}' $mor_boxinfo`
                    boxcolor=`awk -v row=$jj 'NR==row{print $7}' $mor_boxinfo`
                    gmt plot ../inputs/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
                    # cut data and get the mean RMBA in the box
                    gmt grdcut $grav_diff -F../inputs/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    # write to file
                    echo $meanRMBA >>Results/meanRMBA_MOR_${Etas[3]}-${Etas[2]}.txt
                    # plot text of mean rmba
                    lonClatC_box_mor=`awk '{print $1, $2}' ../inputs/BoxC_${boxname}.lonlat`
                    angle_rot_mor=`awk '{print $3}' ../inputs/BoxC_${boxname}.lonlat`                    
                    echo $lonClatC_box_mor $meanRMBA | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor
                done
        gmt subplot end
    gmt end show

rm gmt.history *.grad #grav_*.nc
