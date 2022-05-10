#===== Script for RMBA of mid-ocean ridge-transform system=====
# 01.04.2022

#===== Change part1: Input parameters' value=====
    dataname=Oceanographer
    #Clipperton, ChileRidge, Discovery, Marathon, Marion, Blanco, Oceanographer, Gofar, 
    #Garrett, Vlamingh, Atlantis, Clipperton, MarieCeleste, AtlantisII, Kane, DuTroit 
    Etas=(hsc isov disl vp vep)
    gmt makecpt -C../../../romaO.cpt+h -I -T-15/10/1 -Z >grav_diff.cpt
    gmt makecpt -C../../../romaO.cpt+h -I -T-12/15/1 -Z >grav_diff2.cpt
    gmt makecpt -C../../../roma.cpt+h -I -T-3/3/1 -Z >grav_moho.cpt
    gmt makecpt -C../../../vikO.cpt+h -T-40/40 >grav_therm.cpt
    gmt makecpt -C../../../basecpt_grav.cpt+h -T-30/25 >grav_rmba.cpt
   
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
    mor_index=../inputs/index_ridge.txt   
    otffz_boxinfo=../inputs/otffz_boxinfo.txt
    mor_boxinfo=../inputs/mor_boxinfo.txt
    icoc_boxinfo=../inputs/box_IC_OC.txt
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
function shiftDatabyMean()
{
    grd_data=$1
    meanValue_OTF=`gmt grdmath ${grd_data} MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
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
    file_ridge_index=$mor_index
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
function getNode_Box()
{
    # Role: extract the lon/lat of the vertex
    x0=$1
    y0=$2
    costheta=$3
    sintheta=$4
    len_along=$5  #length along the vector
    len_oth=$6    #length along the orthogonal direction
    # Outputs: X_out & Y_out
    xa=`echo $x0 $costheta $len_along | awk '{print $1+$2*$3}'`
    ya=`echo $y0 $sintheta $len_along | awk '{print $1+$2*$3}'`
    X_out=`echo $x0 $sintheta $len_oth $xa | awk '{print $1-$2*$3 + $4}'`
    Y_out=`echo $y0 $costheta $len_oth $ya | awk '{print $1+$2*$3 + $4}'`
}
function getAverageBox_ICOC()
{
    # Get the four vertex Lon/Lat of IC/OC boxes
    file_xy=$1
    file_box_def=$2
    nBoxes=`awk 'END{print (NR-1)/4}' $file_box_def`
    file_box_lonlat="${file_box_def%.*}".lonlat
    rm $file_box_lonlat
    for (( i=1; i<=$nBoxes; i++ ));
    do 
        ind_start=`awk -v indBox=$i 'NR==(indBox-1)*4+3{print $1}' $file_box_def`
        ind_end=`awk -v indBox=$i 'NR==(indBox-1)*4+3{print $2}' $file_box_def`
        width=`awk -v indBox=$i 'NR==(indBox-1)*4+4{print $1}' $file_box_def`
        height=`awk -v indBox=$i 'NR==(indBox-1)*4+4{print $2}' $file_box_def`
        len2ridge=`awk -v indBox=$i 'NR==(indBox-1)*4+5{print $1}' $file_box_def`
        len2TF=`awk -v indBox=$i 'NR==(indBox-1)*4+5{print $2}' $file_box_def`
        lonlat_start=(`awk -v row=$ind_start 'NR==row{print $1, $2}' $file_xy`)
        lonlat_end=`awk -v row=$ind_end 'NR==row{print $1, $2}' $file_xy`
        lon0=${lonlat_start[0]}
        lat0=${lonlat_start[1]}
        x0=0
        y0=0
        x1y1=(`echo $lonlat_end | cs2cs +proj=latlon +to +proj=tmerc +datum=WGS84 +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0`)
        x1=${x1y1[0]}
        y1=${x1y1[1]}
        length_TF=`echo $x0 $y0 $x1 $y1 | awk '{print sqrt(($3-$1)*($3-$1) + ($4-$2)*($4-$2))}'`
        costheta=`echo $x0 $y0 $x1 $y1 $length_TF | awk '{print ($3-$1)/$5}'`
        sintheta=`echo $x0 $y0 $x1 $y1 $length_TF | awk '{print ($4-$2)/$5}'`
        angle_rot=`echo $sintheta $costheta | awk '{printf "%.1f", atan2($1, $2)/3.141592653*180}'`
        if [[ "`basename $file_box_def`" == "zone_OTF_FZ_"*".txt" ]]; then
            height=`echo $height $length_TF | awk '{print $1+$2}'`
        fi
        # four vertexs
        # 1: closest to the starting point
        getNode_Box $x0 $y0 $costheta $sintheta $len2TF $len2ridge
        echo $X_out $Y_out | cs2cs +proj=latlon +to +proj=tmerc +datum=WGS84 +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 -I -f "%.8f" | awk '{print $1, $2}' >tmp.box
        # 2
        getNode_Box $x0 $y0 $costheta $sintheta $len2TF `echo $len2ridge $width | awk '{print $1+$1/sqrt($1*$1)*$2}'`
        echo $X_out $Y_out | cs2cs +proj=latlon +to +proj=tmerc +datum=WGS84 +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 -I -f "%.8f" | awk '{print $1, $2}' >>tmp.box
        # 3
        getNode_Box $x0 $y0 $costheta $sintheta `echo $len2TF $height | awk '{print $1+$1/sqrt($1*$1)*$2}'` `echo $len2ridge $width | awk '{print $1+$1/sqrt($1*$1)*$2}'`
        echo $X_out $Y_out | cs2cs +proj=latlon +to +proj=tmerc +datum=WGS84 +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 -I -f "%.8f" | awk '{print $1, $2}' >>tmp.box
        # 4
        getNode_Box $x0 $y0 $costheta $sintheta `echo $len2TF $height | awk '{print $1+$1/sqrt($1*$1)*$2}'` $len2ridge
        echo $X_out $Y_out | cs2cs +proj=latlon +to +proj=tmerc +datum=WGS84 +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 -I -f "%.8f" | awk '{print $1, $2}' >>tmp.box
        # centre point of the box
        lonclatc=`awk '{lonc+=$1; latc+=$2}; END {print lonc/NR, latc/NR}' tmp.box`
        # write to file
        awk -v indBox=$i -v lonclatc="$lonclatc" -v angle_rot=${angle_rot} 'NR==(indBox-1)*4+2{print $1, $2, $3, lonclatc, angle_rot}' $file_box_def >>$file_box_lonlat
        cat tmp.box >>$file_box_lonlat
        rm tmp.box
    done
}

if false; then
# 2. Generate all boxes
  # 2.1 OTF, FZ boxes
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
    getAverageBox_ICOC $ROTF ${icoc_boxinfo}    
  # Quick check the box
    gmt begin Results/checkbox pdf
    gmt basemap -R$faa -JM15c -Ba -BWseN
    gmt grdimage $mba
    for box_names in `ls ../inputs/averageBox_*.lonlat`
        do 
        gmt plot ${box_names} -L -W0.25p -Gwhite@50
    done
    gmt plot ../inputs/box_IC_OC.lonlat -L -W0.25p -Gwhite@50
    plotControlTransformFault
    gmt end show
#exit      
fi
# 3. Calculate RMBA field with rheology-dependence mantle thermal model
  # And generate the meanRMBA boxes for OTF, FZ, MOR, IC, and OC
  # Call functions: (1) gravity_xyz2grd(); (2) shiftDataTo0OTF;
  # (3) getAverageBox_TF_FZ; (4)getAverageBox_Ridge; (5) getAverageBox_ICOC    
    echo "OTFname rheology RMBAshift OTF FZ1 FZ2 MOR1 MOR2 IC1 OC1 IC2 OC2" >Results/averageRMBA_${dataname}.txt
    echo "OTFname rheology mean_thermal OTF FZ1 FZ2 MOR1 MOR2 IC1 OC1 IC2 OC2" >Results/averageThermal_${dataname}.txt 
    echo "OTFname rheology OTF FZ1 FZ2 MOR1 MOR2 IC1 OC1 IC2 OC2 OTF_min OTF_max MOR1_min MOR1_max MOR2_min MOR2_max" >Results/averageMoho_${dataname}.txt
    # we consider the hsc model as the reference model
    gravity_xyz2grd grav_hsc
    grav_therm_hsc=grav_hsc.nc
    # reference mean value from the hsc_thermal gravity
    meangrav_thermal=`gmt grdmath $grav_therm_hsc MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
    gmt grdmath $grav_therm_hsc $meangrav_thermal SUB = grav_therm_minusMeanValue_hsc.nc
    gmt grdmath ${mba} grav_therm_minusMeanValue_hsc.nc SUB = rmba_hsc.nc
    shiftDatabyMean rmba_hsc.nc

    for i in {0..4}; do
        #3.1 Calculate gravitational effect of rheology-dependence thermal model and then calculate the RMBA
        etaname=${Etas[i]}
        #Thermal gravity anomaly
        gravity_xyz2grd grav_${etaname}
        grav_therm=grav_${etaname}.nc
        grav_therm_minusMeanValue=grav_${etaname}_minusMean.nc
        # mean value of the thermal gravity
        #meangrav_thermal=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
        gmt grdmath $grav_therm $meangrav_thermal SUB = $grav_therm_minusMeanValue
        gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
        echo $meangrav_thermal >meangrav_${etaname}.txt
        # Calculate the full RMBA (no shift) for moho inversion only
        fullrmba=Results/${dataname}_fullrmba_${etaname}.nc
        rmba=Results/${dataname}_rmba_${etaname}.nc
        gmt grdmath ${mba} $grav_therm_minusMeanValue SUB = ${fullrmba}
        # We shift RMBA with the same shift value in the hsc model (mean value)
        gmt grdmath ${fullrmba} ${meanValue_OTF} SUB = ${rmba}
        gmt grdgradient ${rmba} -A30 -Nt0.6 -Qc -G${rmba}.grad
        echo $meanValue_OTF >meanRMBA_${etaname}.txt
        # Moho
        moho=Results/${dataname}_moho_${etaname}.nc
        # 3.2  Extract average thermal gravity & RMBA & Moho values for all boxes
          # OTF1, get the mean values in the OTF boxes
            gmt grdcut $grav_therm_minusMeanValue -F../inputs/averageBox_OTF1.lonlat -Gtmp_cut_box.nc
            meanThermal_otf1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $rmba -F../inputs/averageBox_OTF1.lonlat -Gtmp_cut_box.nc
            meanRMBA_otf1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $moho -F../inputs/averageBox_OTF1.lonlat -Gtmp_cut_box.nc
            meanMoho_otf1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc`
            minMoho_otf1=`gmt grdinfo tmp_cut_box.nc | grep "v_max" | awk '{printf "%.1f", $3}'`
            maxMoho_otf1=`gmt grdinfo tmp_cut_box.nc | grep "v_max" | awk '{printf "%.1f", $5}' && rm tmp_cut_box.nc`
          # FZ1
            gmt grdcut $grav_therm_minusMeanValue -F../inputs/averageBox_FZ1.lonlat -Gtmp_cut_box.nc
            meanThermal_fz1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $rmba -F../inputs/averageBox_FZ1.lonlat -Gtmp_cut_box.nc
            meanRMBA_fz1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $moho -F../inputs/averageBox_FZ1.lonlat -Gtmp_cut_box.nc
            meanMoho_fz1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`  
          # FZ2
            gmt grdcut $grav_therm_minusMeanValue -F../inputs/averageBox_FZ2.lonlat -Gtmp_cut_box.nc
            meanThermal_fz2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`            
            gmt grdcut $rmba -F../inputs/averageBox_FZ2.lonlat -Gtmp_cut_box.nc
            meanRMBA_fz2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $moho -F../inputs/averageBox_FZ2.lonlat -Gtmp_cut_box.nc
            meanMoho_fz2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`  
          # MOR1
            gmt grdcut $grav_therm_minusMeanValue -F../inputs/averageBox_Ridge1.lonlat -Gtmp_cut_box.nc
            meanThermal_mor1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $rmba -F../inputs/averageBox_Ridge1.lonlat -Gtmp_cut_box.nc
            meanRMBA_mor1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $moho -F../inputs/averageBox_Ridge1.lonlat -Gtmp_cut_box.nc
            meanMoho_mor1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc`
            minMoho_mor1=`gmt grdinfo tmp_cut_box.nc | grep "v_max" | awk '{printf "%.1f", $3}'`
            maxMoho_mor1=`gmt grdinfo tmp_cut_box.nc | grep "v_max" | awk '{printf "%.1f", $5}' && rm tmp_cut_box.nc`
          # MOR2
            gmt grdcut $grav_therm_minusMeanValue -F../inputs/averageBox_Ridge2.lonlat -Gtmp_cut_box.nc
            meanThermal_mor2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $rmba -F../inputs/averageBox_Ridge2.lonlat -Gtmp_cut_box.nc
            meanRMBA_mor2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $moho -F../inputs/averageBox_Ridge2.lonlat -Gtmp_cut_box.nc
            meanMoho_mor2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc`
            minMoho_mor2=`gmt grdinfo tmp_cut_box.nc | grep "v_max" | awk '{printf "%.1f", $3}'`
            maxMoho_mor2=`gmt grdinfo tmp_cut_box.nc | grep "v_max" | awk '{printf "%.1f", $5}' && rm tmp_cut_box.nc`
          # IC1
            nbox_icoc=../inputs/box_IC_OC.lonlat
            awk -v indBox=1 '{if(NR>(1+(indBox-1)*5) && NR<=(indBox*5)){print}}' $nbox_icoc >tmp_ic1.box
            gmt grdcut $grav_therm_minusMeanValue -Ftmp_ic1.box -Gtmp_cut_box.nc
            meanThermal_ic1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $rmba -Ftmp_ic1.box -Gtmp_cut_box.nc
            meanRMBA_ic1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $moho -Ftmp_ic1.box -Gtmp_cut_box.nc
            meanMoho_ic1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
          # OC1
            awk -v indBox=2 '{if(NR>(1+(indBox-1)*5) && NR<=(indBox*5)){print}}' $nbox_icoc >tmp_oc1.box
            gmt grdcut $grav_therm_minusMeanValue -Ftmp_oc1.box -Gtmp_cut_box.nc
            meanThermal_oc1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $rmba -Ftmp_oc1.box -Gtmp_cut_box.nc
            meanRMBA_oc1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $moho -Ftmp_oc1.box -Gtmp_cut_box.nc
            meanMoho_oc1=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
          # IC2
            nbox_icoc=../inputs/box_IC_OC.lonlat
            awk -v indBox=3 '{if(NR>(1+(indBox-1)*5) && NR<=(indBox*5)){print}}' $nbox_icoc >tmp_ic2.box
            gmt grdcut $grav_therm_minusMeanValue -Ftmp_ic2.box -Gtmp_cut_box.nc
            meanThermal_ic2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $rmba -Ftmp_ic2.box -Gtmp_cut_box.nc
            meanRMBA_ic2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $moho -Ftmp_ic2.box -Gtmp_cut_box.nc
            meanMoho_ic2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
          # OC2
            awk -v indBox=4 '{if(NR>(1+(indBox-1)*5) && NR<=(indBox*5)){print}}' $nbox_icoc >tmp_oc2.box
            gmt grdcut $grav_therm_minusMeanValue -Ftmp_oc2.box -Gtmp_cut_box.nc
            meanThermal_oc2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $rmba -Ftmp_oc2.box -Gtmp_cut_box.nc
            meanRMBA_oc2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
            gmt grdcut $moho -Ftmp_oc2.box -Gtmp_cut_box.nc
            meanMoho_oc2=`gmt grdmath tmp_cut_box.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box.nc`
        # write to file
        echo ${dataname} ${etaname} ${meanValue_OTF} ${meanRMBA_otf1} ${meanRMBA_fz1} ${meanRMBA_fz2} ${meanRMBA_mor2} ${meanRMBA_mor1} ${meanRMBA_ic2} ${meanRMBA_oc2} ${meanRMBA_ic1} ${meanRMBA_oc1} >>Results/averageRMBA_${dataname}.txt
        echo ${dataname} ${etaname} ${meangrav_thermal} ${meanThermal_otf1} ${meanThermal_fz1} ${meanThermal_fz2} ${meanThermal_mor2} ${meanThermal_mor1} ${meanThermal_ic2} ${meanThermal_oc2} ${meanThermal_ic1} ${meanThermal_oc1} >>Results/averageThermal_${dataname}.txt        
        echo ${dataname} ${etaname} ${meanMoho_otf1} ${meanMoho_fz1} ${meanMoho_fz2} ${meanMoho_mor2} ${meanMoho_mor1} ${meanMoho_ic2} ${meanMoho_oc2} ${meanMoho_ic1} ${meanMoho_oc1} ${minMoho_otf1} ${maxMoho_otf1} ${minMoho_mor2} ${maxMoho_mor2} ${minMoho_mor1} ${maxMoho_mor1} >>Results/averageMoho_${dataname}.txt
        # Check ranges of thermal gravity, RMBA, and Moho for the figure
        zmin_therm=`gmt grdinfo $grav_therm_minusMeanValue | grep "v_max" | awk '{printf "%.1f", $3}'`
        zmax_therm=`gmt grdinfo $grav_therm_minusMeanValue | grep "v_max" | awk '{printf "%.1f", $5}'`
        echo ${etaname} "Thermal - min: $zmin_therm mGal, max: $zmax_therm mGal"
        zmin_rmba=`gmt grdinfo $rmba | grep "v_max" | awk '{printf "%.1f", $3}'`
        zmax_rmba=`gmt grdinfo $rmba | grep "v_max" | awk '{printf "%.1f", $5}'`
        echo ${etaname} "RMBA - min: $zmin_rmba mGal, max: $zmax_rmba mGal"
        zmin_moho=`gmt grdinfo $moho | grep "v_max" | awk '{printf "%.1f", $3}'`
        zmax_moho=`gmt grdinfo $moho | grep "v_max" | awk '{printf "%.1f", $5}'`
        echo ${etaname} "Moho - min: $zmin_moho km, max: $zmax_moho km"
    done

sed 's/[ ][ ]*/,/g' Results/averageRMBA_${dataname}.txt > Results/averageRMBA_${dataname}.csv
sed 's/[ ][ ]*/,/g' Results/averageThermal_${dataname}.txt > Results/averageThermal_${dataname}.csv
sed 's/[ ][ ]*/,/g' Results/averageMoho_${dataname}.txt > Results/averageMoho_${dataname}.csv     
#exit
#fi

# 4. Plot the figure
  # Call the functions: (1) plotControlTransformFault; (2) makecpt_grd;
  # (3) makecpt_grd_basecpt  
  # Figure includes 1st row (5): gravity of rheologic-dependent thermal model
  # 2nd row (5): RMBA for rheoglical differing models
  # 3rd row (5): RMBA difference betwen them: a. disl-isov
  # b. vp-isov, c. vep-isov, d. vp-disl, e. vep-disl
  # 4th row (5): RMBA difference betweem: a. isov-hsc, b. disl-hsc, 
  # c. vp-hsc, vep-hsc, e. vep-vp    
    model1=Half_space_cooling    
    model2=Isoviscous
    model3=Viscous
    model4=Visco-plastic
    model5=Visco-elasto-plastic
      gmt begin Results/${dataname}_gravity_pattern pdf
        gmt subplot begin 5x5 -A+JTL+o0.5c -Fs15c/9c -M0.5c/1.5c -R$faa -JM15c -Ba1df30m -BWSne -T"Gravity anomaly estimated from rheology-dependence thermal models: "$dataname
            gmt subplot set 0,0
                modelname=$model1
                etaname=${Etas[0]}
                #Thermal gravity anomaly -hsc
                #makecpt_grd_basecpt grav_${etaname}_minusMean.nc
                gmt grdimage grav_${etaname}_minusMean.nc -BwsEN -Cgrav_therm.cpt #-I${grav_therm_minusMeanValue}.grad
                gmt grdcontour grav_${etaname}_minusMean.nc -C10
                meanGrav=`awk '{printf "%.1f", $1}' meangrav_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal" -Cgrav_therm.cpt
                plotControlTransformFault
            gmt subplot set 0,1
                modelname=$model2
                etaname=${Etas[1]}
                #Thermal gravity anomaly -isov
                #makecpt_grd_basecpt grav_${etaname}_minusMean.nc
                gmt grdimage grav_${etaname}_minusMean.nc -BwsEN -Cgrav_therm.cpt #-I${grav_therm_minusMeanValue}.grad
                gmt grdcontour grav_${etaname}_minusMean.nc -C10
                meanGrav=`awk '{printf "%.1f", $1}' meangrav_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal" -Cgrav_therm.cpt
                plotControlTransformFault
            gmt subplot set 0,2
                modelname=$model3
                etaname=${Etas[2]}
                #Thermal gravity anomaly -disl
                #makecpt_grd_basecpt grav_${etaname}_minusMean.nc
                gmt grdimage grav_${etaname}_minusMean.nc -BwsEN -Cgrav_therm.cpt #-I${grav_therm_minusMeanValue}.grad
                gmt grdcontour grav_${etaname}_minusMean.nc -C10
                meanGrav=`awk '{printf "%.1f", $1}' meangrav_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal" -Cgrav_therm.cpt
                plotControlTransformFault                
            gmt subplot set 0,3
                modelname=$model4
                etaname=${Etas[3]}
                #Thermal gravity anomaly -vp
                #makecpt_grd_basecpt grav_${etaname}_minusMean.nc
                gmt grdimage grav_${etaname}_minusMean.nc -BwsEN -Cgrav_therm.cpt #-I${grav_therm_minusMeanValue}.grad
                gmt grdcontour grav_${etaname}_minusMean.nc -C10
                meanGrav=`awk '{printf "%.1f", $1}' meangrav_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal" -Cgrav_therm.cpt
                plotControlTransformFault
            gmt subplot set 0,4
                modelname=$model5
                etaname=${Etas[4]}
                #Thermal gravity anomaly -vep
                #makecpt_grd_basecpt grav_${etaname}_minusMean.nc
                gmt grdimage grav_${etaname}_minusMean.nc -BwsEN -Cgrav_therm.cpt #-I${grav_therm_minusMeanValue}.grad
                gmt grdcontour grav_${etaname}_minusMean.nc -C10
                meanGrav=`awk '{printf "%.1f", $1}' meangrav_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal" -Cgrav_therm.cpt
                plotControlTransformFault
              # plot the all boxes
                gmt plot ../inputs/averageBox_OTF1.lonlat -W0.5p,red -L -Gwhite@50
                gmt plot ../inputs/averageBox_FZ1.lonlat -W0.5p,blue -L -Gwhite@50
                gmt plot ../inputs/averageBox_FZ2.lonlat -W0.5p,blue -L -Gwhite@50
                gmt plot ../inputs/averageBox_Ridge1.lonlat -W0.5p,purple -L -Gwhite@50
                gmt plot ../inputs/averageBox_Ridge2.lonlat -W0.5p,purple -L -Gwhite@50
                gmt plot tmp_ic1.box -W0.5p,white -L -Gwhite@50
                gmt plot tmp_oc1.box -W0.5p,black -L -Gwhite@50
                gmt plot tmp_ic2.box -W0.5p,white -L -Gwhite@50
                gmt plot tmp_oc2.box -W0.5p,black -L -Gwhite@50                
                # plot text of mean thermal
                lonClatC_box_otf1=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_OTF1.lonlat`
                angle_rot_otf1=`awk '{print $3}' ../inputs/BoxC_averageBox_OTF1.lonlat`
                echo $lonClatC_box_otf1 $meanThermal_otf1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otf1              
                #fz
                lonClatC_box_FZ1=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_FZ1.lonlat`
                angle_rot_FZ1=`awk '{print $3}' ../inputs/BoxC_averageBox_FZ1.lonlat`
                echo $lonClatC_box_FZ1 $meanThermal_fz1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_FZ1
                lonClatC_box_FZ2=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_FZ2.lonlat`
                angle_rot_FZ2=`awk '{print $3}' ../inputs/BoxC_averageBox_FZ2.lonlat`
                echo $lonClatC_box_FZ2 $meanThermal_fz2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_FZ2
                # MOR
                lonClatC_box_mor1=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_Ridge1.lonlat`
                angle_rot_mor1=`awk '{print $3}' ../inputs/BoxC_averageBox_Ridge1.lonlat`
                echo $lonClatC_box_mor1 $meanThermal_mor1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor1
                lonClatC_box_mor2=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_Ridge2.lonlat`
                angle_rot_mor2=`awk '{print $3}' ../inputs/BoxC_averageBox_Ridge2.lonlat`
                echo $lonClatC_box_mor2 $meanThermal_mor2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor2
                #IC/OC
                lonClatC_ic1=`awk -v indBox=1 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_ic1=`awk -v indBox=1 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_ic1 $meanThermal_ic1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_ic1               
                lonClatC_oc1=`awk -v indBox=2 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_oc1=`awk -v indBox=2 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_oc1 $meanThermal_oc1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_oc1
                lonClatC_ic2=`awk -v indBox=3 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_ic2=`awk -v indBox=3 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_ic2 $meanThermal_ic2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_ic2
                lonClatC_oc2=`awk -v indBox=4 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_oc2=`awk -v indBox=4 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_oc2 $meanThermal_oc2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_oc2               
            # 2nd row --- RMBA maps
            gmt subplot set 1,0
                modelname=$model1
                etaname=${Etas[0]}
                rmba=Results/${dataname}_rmba_${etaname}.nc
                gmt basemap -BwseN -Ba 
                #makecpt_grd $rmba
                gmt grdimage $rmba -Cgrav_rmba.cpt #-I${rmba}.grad
                #gmt grdcontour $rmba -C10
                #gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                shiftRMBA=`awk '{printf "%.1f", $1}' meanRMBA_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $shiftRMBA" -By+l"mGal" -Cgrav_rmba.cpt #-G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 1,1
                modelname=$model2
                etaname=${Etas[1]}
                rmba=Results/${dataname}_rmba_${etaname}.nc
                gmt basemap -BwseN -Ba 
                #makecpt_grd $rmba
                gmt grdimage $rmba -Cgrav_rmba.cpt #-I${rmba}.grad
                #gmt grdcontour $rmba -C10
                #gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                shiftRMBA=`awk '{printf "%.1f", $1}' meanRMBA_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $shiftRMBA" -By+l"mGal" -Cgrav_rmba.cpt #-G${data_min}/${data_max}
                plotControlTransformFault                
            gmt subplot set 1,2
                modelname=$model3
                etaname=${Etas[2]}
                rmba=Results/${dataname}_rmba_${etaname}.nc
                gmt basemap -BwseN -Ba 
                #makecpt_grd $rmba
                gmt grdimage $rmba -Cgrav_rmba.cpt #-I${rmba}.grad
                #gmt grdcontour $rmba -C10
                #gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                shiftRMBA=`awk '{printf "%.1f", $1}' meanRMBA_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $shiftRMBA" -By+l"mGal" -Cgrav_rmba.cpt #-G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 1,3
                modelname=$model4
                etaname=${Etas[3]}
                rmba=Results/${dataname}_rmba_${etaname}.nc
                gmt basemap -BwseN -Ba 
                #makecpt_grd $rmba
                gmt grdimage $rmba -Cgrav_rmba.cpt #-I${rmba}.grad
                #gmt grdcontour $rmba -C10
                #gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                shiftRMBA=`awk '{printf "%.1f", $1}' meanRMBA_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $shiftRMBA" -By+l"mGal" -Cgrav_rmba.cpt #-G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 1,4
                modelname=$model5
                etaname=${Etas[4]}
                rmba=Results/${dataname}_rmba_${etaname}.nc
                gmt basemap -BwseN -Ba 
                #makecpt_grd $rmba
                gmt grdimage $rmba -Cgrav_rmba.cpt #-I${rmba}.grad
                #gmt grdcontour $rmba -C10
                #gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                shiftRMBA=`awk '{printf "%.1f", $1}' meanRMBA_${etaname}.txt`
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $shiftRMBA" -By+l"mGal" -Cgrav_rmba.cpt #-G${data_min}/${data_max}
              # plot the all boxes
                gmt plot ../inputs/averageBox_OTF1.lonlat -W0.5p,red -L -Gwhite@50
                gmt plot ../inputs/averageBox_FZ1.lonlat -W0.5p,blue -L -Gwhite@50
                gmt plot ../inputs/averageBox_FZ2.lonlat -W0.5p,blue -L -Gwhite@50
                gmt plot ../inputs/averageBox_Ridge1.lonlat -W0.5p,purple -L -Gwhite@50
                gmt plot ../inputs/averageBox_Ridge2.lonlat -W0.5p,purple -L -Gwhite@50
                gmt plot tmp_ic1.box -W0.5p,white -L -Gwhite@50
                gmt plot tmp_oc1.box -W0.5p,black -L -Gwhite@50
                gmt plot tmp_ic2.box -W0.5p,white -L -Gwhite@50
                gmt plot tmp_oc2.box -W0.5p,black -L -Gwhite@50                
                # plot text of mean rmba
                lonClatC_box_otf1=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_OTF1.lonlat`
                angle_rot_otf1=`awk '{print $3}' ../inputs/BoxC_averageBox_OTF1.lonlat`
                echo $lonClatC_box_otf1 $meanRMBA_otf1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otf1              
                #fz
                lonClatC_box_FZ1=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_FZ1.lonlat`
                angle_rot_FZ1=`awk '{print $3}' ../inputs/BoxC_averageBox_FZ1.lonlat`
                echo $lonClatC_box_FZ1 $meanRMBA_fz1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_FZ1
                lonClatC_box_FZ2=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_FZ2.lonlat`
                angle_rot_FZ2=`awk '{print $3}' ../inputs/BoxC_averageBox_FZ2.lonlat`
                echo $lonClatC_box_FZ2 $meanRMBA_fz2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_FZ2
                # MOR
                lonClatC_box_mor1=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_Ridge1.lonlat`
                angle_rot_mor1=`awk '{print $3}' ../inputs/BoxC_averageBox_Ridge1.lonlat`
                echo $lonClatC_box_mor1 $meanRMBA_mor1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor1
                lonClatC_box_mor2=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_Ridge2.lonlat`
                angle_rot_mor2=`awk '{print $3}' ../inputs/BoxC_averageBox_Ridge2.lonlat`
                echo $lonClatC_box_mor2 $meanRMBA_mor2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor2
                #IC/OC
                lonClatC_ic1=`awk -v indBox=1 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_ic1=`awk -v indBox=1 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_ic1 $meanRMBA_ic1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_ic1               
                lonClatC_oc1=`awk -v indBox=2 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_oc1=`awk -v indBox=2 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_oc1 $meanRMBA_oc1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_oc1
                lonClatC_ic2=`awk -v indBox=3 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_ic2=`awk -v indBox=3 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_ic2 $meanRMBA_ic2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_ic2
                lonClatC_oc2=`awk -v indBox=4 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_oc2=`awk -v indBox=4 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_oc2 $meanRMBA_oc2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_oc2
                plotControlTransformFault
            # 3rd row --- Moho maps
            gmt subplot set 2,0
                modelname=$model1
                etaname=${Etas[0]}
                moho=Results/${dataname}_moho_${etaname}.nc
                gmt basemap -BwseN -Ba 
                gmt grdimage $moho -Cgrav_moho.cpt #-I${moho}.grad
                #gmt grdcontour $moho -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Moho depth ($modelname) " -By+l"km" -Cgrav_moho.cpt #-G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 2,1
                modelname=$model2
                etaname=${Etas[1]}
                moho=Results/${dataname}_moho_${etaname}.nc
                gmt basemap -BwseN -Ba 
                gmt grdimage $moho -Cgrav_moho.cpt #-I${moho}.grad
                #gmt grdcontour $moho -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Moho depth ($modelname) " -By+l"km" -Cgrav_moho.cpt #-G${data_min}/${data_max}
                plotControlTransformFault            
            gmt subplot set 2,2
                modelname=$model3
                etaname=${Etas[2]}
                moho=Results/${dataname}_moho_${etaname}.nc
                gmt basemap -BwseN -Ba 
                gmt grdimage $moho -Cgrav_moho.cpt #-I${moho}.grad
                #gmt grdcontour $moho -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Moho depth ($modelname) " -By+l"km" -Cgrav_moho.cpt #-G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 2,3
                modelname=$model4
                etaname=${Etas[3]}
                moho=Results/${dataname}_moho_${etaname}.nc
                gmt basemap -BwseN -Ba 
                gmt grdimage $moho -Cgrav_moho.cpt #-I${moho}.grad
                #gmt grdcontour $moho -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Moho depth ($modelname) " -By+l"km" -Cgrav_moho.cpt #-G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 2,4
                modelname=$model5
                etaname=${Etas[4]}
                moho=Results/${dataname}_moho_${etaname}.nc
                gmt basemap -BwseN -Ba 
                gmt grdimage $moho -Cgrav_moho.cpt #-I${moho}.grad
                #gmt grdcontour $moho -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Moho depth ($modelname) " -By+l"km" -Cgrav_moho.cpt #-G${data_min}/${data_max}
                plotControlTransformFault
              # plot the all boxes
                gmt plot ../inputs/averageBox_OTF1.lonlat -W0.5p,red -L -Gwhite@50
                gmt plot ../inputs/averageBox_FZ1.lonlat -W0.5p,blue -L -Gwhite@50
                gmt plot ../inputs/averageBox_FZ2.lonlat -W0.5p,blue -L -Gwhite@50
                gmt plot ../inputs/averageBox_Ridge1.lonlat -W0.5p,purple -L -Gwhite@50
                gmt plot ../inputs/averageBox_Ridge2.lonlat -W0.5p,purple -L -Gwhite@50
                gmt plot tmp_ic1.box -W0.5p,white -L -Gwhite@50
                gmt plot tmp_oc1.box -W0.5p,black -L -Gwhite@50
                gmt plot tmp_ic2.box -W0.5p,white -L -Gwhite@50
                gmt plot tmp_oc2.box -W0.5p,black -L -Gwhite@50                
                # plot text of mean moho
                lonClatC_box_otf1=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_OTF1.lonlat`
                angle_rot_otf1=`awk '{print $3}' ../inputs/BoxC_averageBox_OTF1.lonlat`
                echo $lonClatC_box_otf1 $meanMoho_otf1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_otf1              
                #fz
                lonClatC_box_FZ1=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_FZ1.lonlat`
                angle_rot_FZ1=`awk '{print $3}' ../inputs/BoxC_averageBox_FZ1.lonlat`
                echo $lonClatC_box_FZ1 $meanMoho_fz1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_FZ1
                lonClatC_box_FZ2=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_FZ2.lonlat`
                angle_rot_FZ2=`awk '{print $3}' ../inputs/BoxC_averageBox_FZ2.lonlat`
                echo $lonClatC_box_FZ2 $meanMoho_fz2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_FZ2
                # MOR
                lonClatC_box_mor1=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_Ridge1.lonlat`
                angle_rot_mor1=`awk '{print $3}' ../inputs/BoxC_averageBox_Ridge1.lonlat`
                echo $lonClatC_box_mor1 $meanMoho_mor1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor1
                lonClatC_box_mor2=`awk '{print $1, $2}' ../inputs/BoxC_averageBox_Ridge2.lonlat`
                angle_rot_mor2=`awk '{print $3}' ../inputs/BoxC_averageBox_Ridge2.lonlat`
                echo $lonClatC_box_mor2 $meanMoho_mor2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_mor2
                #IC/OC
                lonClatC_ic1=`awk -v indBox=1 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_ic1=`awk -v indBox=1 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_ic1 $meanMoho_ic1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_ic1               
                lonClatC_oc1=`awk -v indBox=2 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_oc1=`awk -v indBox=2 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_oc1 $meanMoho_oc1 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_oc1
                lonClatC_ic2=`awk -v indBox=3 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_ic2=`awk -v indBox=3 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_ic2 $meanMoho_ic2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_ic2
                lonClatC_oc2=`awk -v indBox=4 'NR==(1+(indBox-1)*5){print $4, $5}' $nbox_icoc`
                angle_rot_oc2=`awk -v indBox=4 'NR==(1+(indBox-1)*5){if($6>90 || $6<-90){print $6-180}else{print $6}}' $nbox_icoc`
                echo $lonClatC_oc2 $meanMoho_oc2 | gmt text -F+f12p,Helvetica-Bold,black+a$angle_rot_oc2
                plotControlTransformFault
            # 4th row --- RMBA difference between different rheologies   
            gmt subplot set 3,0
                m1=$model1
                m2=$model2
                grav_therm1=grav_${Etas[0]}.nc
                grav_therm2=grav_${Etas[1]}.nc
                # RMBA_diff_isov-hsc = grav_diff_hsc-isov!
                grav_diff=rmba_diff_${Etas[1]}-${Etas[0]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
                plotControlTransformFault                
            gmt subplot set 3,1
                m1=$model1
                m2=$model3
                grav_therm1=grav_${Etas[0]}.nc
                grav_therm2=grav_${Etas[2]}.nc
                grav_diff=rmba_diff_${Etas[2]}-${Etas[0]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
                plotControlTransformFault
            gmt subplot set 3,2
                m1=$model1
                m2=$model4
                grav_therm1=grav_${Etas[0]}.nc
                grav_therm2=grav_${Etas[3]}.nc
                grav_diff=rmba_diff_${Etas[3]}-${Etas[0]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
                plotControlTransformFault
            gmt subplot set 3,3
                m1=$model1
                m2=$model5
                grav_therm1=grav_${Etas[0]}.nc
                grav_therm2=grav_${Etas[4]}.nc
                grav_diff=rmba_diff_${Etas[4]}-${Etas[0]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
                plotControlTransformFault
            gmt subplot set 3,4
                m1=$model4
                m2=$model5
                grav_therm1=grav_${Etas[3]}.nc
                grav_therm2=grav_${Etas[4]}.nc
                grav_diff=rmba_diff_${Etas[4]}-${Etas[3]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
                plotControlTransformFault                                
            # 5th row --- RMBA difference between different rheologies   
            gmt subplot set 4,0
                m1=$model2
                m2=$model3
                grav_therm1=grav_${Etas[1]}.nc
                grav_therm2=grav_${Etas[2]}.nc
                # RMBA_diff_disl-isov = grav_diff_isov-disl!
                grav_diff=rmba_diff_${Etas[2]}-${Etas[1]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff2.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff2.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20               
                plotControlTransformFault                
            gmt subplot set 4,1
                m1=$model2
                m2=$model4
                grav_therm1=grav_${Etas[1]}.nc
                grav_therm2=grav_${Etas[3]}.nc
                # RMBA_diff_disl-isov = grav_diff_isov-disl!
                grav_diff=rmba_diff_${Etas[3]}-${Etas[1]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff2.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff2.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20               
                plotControlTransformFault
            gmt subplot set 4,2
                m1=$model2
                m2=$model5
                grav_therm1=grav_${Etas[1]}.nc
                grav_therm2=grav_${Etas[4]}.nc
                # RMBA_diff_disl-isov = grav_diff_isov-disl!
                grav_diff=rmba_diff_${Etas[4]}-${Etas[1]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff2.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff2.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20               
                plotControlTransformFault
            gmt subplot set 4,3
                m1=$model3
                m2=$model4
                grav_therm1=grav_${Etas[2]}.nc
                grav_therm2=grav_${Etas[3]}.nc
                grav_diff=rmba_diff_${Etas[3]}-${Etas[2]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff2.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff2.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20               
                plotControlTransformFault
            gmt subplot set 4,4
                m1=$model3
                m2=$model5
                grav_therm1=grav_${Etas[2]}.nc
                grav_therm2=grav_${Etas[4]}.nc
                grav_diff=rmba_diff_${Etas[4]}-${Etas[2]}.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -Cgrav_diff2.cpt #-I${grav_diff}.grad 
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: $m2 - $m1" -By+l"mGal" -Cgrav_diff2.cpt 
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_max" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20               
                plotControlTransformFault                                                        
        gmt subplot end
    gmt end show

rm gmt.* *.grad tmp* mean* *.cpt grav_*.nc rmba*.nc Results/*.grad Results/average*.txt
