. baseInfo.sh
. gmt_shell_functions.sh
figwidth=16
# figheight=`gmt_get_map_height -R$bathy  -JM${figwidth}c`
# In case the picture is too large, for example AtlantisII
# figwidth=`echo $figwidth $figheight | awk '{if($2/$1 > 1){ print $1/($2/$1) } else { print $1 } }'`
figfmt="pdf"
thermalModelType="" # _multiOTF, or empty
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
# UTM CS
function plotLengthAngle_ROTF()
{
    file_ROTF=$1
    num_RTI=`gmt gmtinfo $file_ROTF | awk '{print $4}'`
    # 1. calculate spreading rate vector
    file_ridge_index=${dataname}/input/ridge_index.txt
    # if [ -f $file_ridge_index ]; then 
    #     num_ridge=`gmt gmtinfo $file_ridge_index | awk '{print $4}'`
    #     for ((i=1; i<=num_ridge; i++))
    #     do 
    #         i_start=`awk -v row=$i 'NR==(row){print $1}' $file_ridge_index`
    #         i_end=`awk -v row=$i 'NR==(row){print $2}' $file_ridge_index`
    #         echo "Ridge $i: "$i_start $i_end
    #         point1_orig=`awk -v row=$i_start 'NR==(row){print $1, $2}' $file_ROTF` 
    #         point2_orig=`awk -v row=$i_end 'NR==(row){print $1, $2}' $file_ROTF` 
    #         # make sure y of point1 less than y of point 2
    #         point1=`echo $point1_orig $point2_orig | awk '{if($2>$4){print $3, $4}else{print $1, $2}}'`
    #         point2=`echo $point1_orig $point2_orig | awk '{if($2>$4){print $1, $2}else{print $3, $4}}'`
    #         # right arrow
    #         nv=`echo $point1 $point2 | awk '{len=sqrt(($1-$3)*($1-$3) + ($2-$4)*($2-$4)); print ($3-$1)/len, ($4-$2)/len}'`
    #         costheta_xaxis=`echo $nv | awk '{print $1}'`
    #         theta_xaxis=`echo $costheta_xaxis | awk '{printf "%f", atan2(sqrt(1-$1*$1), $1)/3.141592653*180}'`
    #         alpha=`echo $theta_xaxis | awk '{print $1-90}'`
    #         offset=5000
    #         x0y0=`echo $point1 $point2 $alpha $offset | awk '{print ($1+$3)/2.0+cos($5/180.0*3.141592653)*$6, ($2+$4)/2.0+sin($5/180.0*3.141592653)*$6}'`
    #         echo "$x0y0 $alpha 2c" | gmt plot -Sv0.4c/0.5c/0.4c -W1p,red -Gyellow --MAP_VECTOR_SHAPE=0.5 
    #         offset_text=10000
    #         x0y0_text=`echo $point1 $point2 $alpha $offset_text | awk '{print ($1+$3)/2.0+cos($5/180.0*3.141592653)*$6, ($2+$4)/2.0+sin($5/180.0*3.141592653)*$6}'`
    #         text_U=`echo $alpha | awk '{theta=$1/180.0*3.141592653; printf "(%.3f, %.3f) U",cos(theta),sin(theta)}'`
    #         echo "$x0y0_text 6p,Times-Roman $alpha LM $text_U" | gmt text -Dj0p -F+f+a+j
    #         # left arrow
    #         alpha2=`echo $alpha | awk '{print $1+180}'`
    #         x0y0=`echo $point1 $point2 $alpha2 $offset | awk '{print ($1+$3)/2.0+cos($5/180.0*3.141592653)*$6, ($2+$4)/2.0+sin($5/180.0*3.141592653)*$6}'`
    #         echo "$x0y0 $alpha2 2c" | gmt plot -Sv0.4c/0.5c/0.4c -W1p,red -Gyellow --MAP_VECTOR_SHAPE=0.5 
    #         x0y0_text=`echo $point1 $point2 $alpha2 $offset_text | awk '{print ($1+$3)/2.0+cos($5/180.0*3.141592653)*$6, ($2+$4)/2.0+sin($5/180.0*3.141592653)*$6}'`
    #         text_U=`echo $alpha2 | awk '{theta=$1/180.0*3.141592653; printf "(%.3f, %.3f) U",cos(theta),sin(theta)}'`
    #         echo "$x0y0_text 6p,Times-Roman $alpha RM $text_U" | gmt text -Dj0p -F+f+a+j
    #     done
    # fi
    # 2. calculate and plot length
    for ((i=1; i<num_RTI; i++))
    do 
        point1=`awk -v row=$i 'NR==row{print $1, $2}' $file_ROTF` 
        point2=`awk -v row=$i 'NR==(row+1){print $1, $2}' $file_ROTF` 
        length=`echo $point1 $point2 | awk '{printf "%.0f", sqrt(($1-$3)*($1-$3) + ($2-$4)*($2-$4))/1000}'`
        x0y0=`echo $point1 $point2 | awk '{print ($1+$3)/2.0, ($2+$4)/2.0}'`
        alpha=`echo $point1 $point2 | awk '{alpha=atan2($2-$4, $1-$3)/3.141592653*180;if(sqrt((alpha-180)*(alpha-180))<0.1)alpha=0; if(alpha<-90){alpha+=180};printf "%.4f", alpha}'`
        echo "$x0y0 3p,Times-Roman $alpha CM ${length} km " | gmt text -Gcyan@10 -Dj0p -F+f+a+j
    done
    # 3. calculate angle begween ridge and OTF
    for ((i=2; i<num_RTI; i++))
    do 
        point1=`awk -v row=$i 'NR==(row-1){print $1, $2}' $file_ROTF` 
        point0=`awk -v row=$i 'NR==(row){print $1, $2}' $file_ROTF` 
        point2=`awk -v row=$i 'NR==(row+1){print $1, $2}' $file_ROTF` 
        nv1=`echo $point0 $point1 | awk '{len=sqrt(($1-$3)*($1-$3) + ($2-$4)*($2-$4)); print ($3-$1)/len, ($4-$2)/len}'`
        nv2=`echo $point0 $point2 | awk '{len=sqrt(($1-$3)*($1-$3) + ($2-$4)*($2-$4)); print ($3-$1)/len, ($4-$2)/len}'`
        costheta=`echo $nv1 $nv2 | awk '{print ($1*$3)+($2*$4)}'`
        theta=`echo $costheta | awk '{printf "%.0f", atan2(sqrt(1-$1*$1), $1)/3.141592653*180}'`
        align_h=`echo $nv1 $nv2 | awk '{x0=($1+$3)/2.0; if(x0>0){print "R"}else{print "L"}}'`
        align_v=`echo $nv1 $nv2 | awk '{x0=($2+$4)/2.0; if(x0>0){print "T"}else{print "B"}}'`
        echo "$point0 3p,Times-Roman 0 ${align_h}${align_v} ${theta}@. " | gmt text -Ggreen@10 -Dj0p -F+f+a+j
    done
}
function transformCoordinates()
{
    controlPointsFile_lonlat=$1
    translate_x=$2
    translate_y=$3
    z0=$4
    # 0. convert this transform fault line parallel to x axis in ASPECT modeling
    # suggest: set one endpoint of the transform fault as (0,0) in ASPECT modeling coordinate system !!
    lon0=`awk 'NR==1{print $1}' ${OTF}` 
    lat0=`awk 'NR==1{print $2}' ${OTF}` 
    x0=0
    y0=0
    cat $controlPointsFile_lonlat | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 >${dataname}_OTF_UTM.txt
    # awk 'NR==2{print "length of control line: "sqrt($1*$1+$2*$2)}' input/${dataname}_OTF_UTM.txt
    cat $ROTF | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 >tmp_ROTF_UTM.txt
    # calculate cos(theta) and sin(theta)
    costheta=`awk 'NR==2{print $1/sqrt($1*$1+$2*$2)}' ${dataname}_OTF_UTM.txt`
    sintheta=`awk 'NR==2{print $2/sqrt($1*$1+$2*$2)}' ${dataname}_OTF_UTM.txt`
    angle_rot=`echo $sintheta $costheta | awk '{printf "%.1f", atan2($1, $2)/3.141592653*180}'`
    echo "Rotation angle: "${angle_rot}" degree"
    # 1. get x, y coordinate of all gravity calculation points
    gmt grd2xyz $faa | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 >${dataname}_UTM.xyz
    # get range of data points
    xlim=`gmt gmtinfo ${dataname}_UTM.xyz | awk '{print $5}' | cut -d "<" -f 2 | cut -d ">" -f 1`
    ylim=`gmt gmtinfo ${dataname}_UTM.xyz | awk '{print $6}' | cut -d "<" -f 2 | cut -d ">" -f 1`
    range_UTM=$xlim/$ylim
    # 2. rotate x,y
    awk -v co=$costheta -v si=$sintheta '{print co*$1+si*$2, -si*$1+co*$2, $3}' ${dataname}_UTM.xyz >${dataname}_UTM_rot.xyz
    awk -v co=$costheta -v si=$sintheta '{print co*$1+si*$2, -si*$1+co*$2, $3}' ${dataname}_OTF_UTM.txt >${dataname}_OTF_UTM_rot.txt
    awk -v co=$costheta -v si=$sintheta '{print co*$1+si*$2, -si*$1+co*$2, $3}' tmp_ROTF_UTM.txt >tmp_ROTF_UTM_rot.txt
    length_transform=`awk 'NR==2{printf "%.0f", sqrt($1*$1+$2*$2)/1000}' ${dataname}_OTF_UTM_rot.txt`
    echo "length of control line after transforming: "$length_transform" km"
    # translate 
    awk -v dx=$translate_x -v dy=$translate_y '{print $1+dx, $2+dy, $3}' ${dataname}_UTM_rot.xyz >${dataname}_UTM_rot_translate.xyz
    awk -v dx=$translate_x -v dy=$translate_y '{print $1+dx, $2+dy, $3}' ${dataname}_OTF_UTM_rot.txt >${dataname}_OTF_UTM_rot_translate.txt
    awk -v dx=$translate_x -v dy=$translate_y '{print $1+dx, $2+dy, $3}' tmp_ROTF_UTM_rot.txt >tmp_ROTF_UTM_rot_translate.txt
    awk -v dx=$translate_x -v dy=$translate_y '{print $1+dx, $2+dy, $3}' ${dataname}_OTF_UTM.txt >${dataname}_OTF_UTM_translate.txt
    # output sites.xyz
    awk -v z0=$z0 '{print $1, $2, z0}' ${dataname}_UTM_rot_translate.xyz >${dataname}/input/${dataname}_sites_${z0}.xyz

    # plot gravity calculation points before and after transformation
    figwidth_orig=$figwidth
    figheight_orig=`gmt_get_map_height -R$bathy  -JM${figwidth_orig}c`
    # in case figure too big
    figwidth_orig=`echo $figwidth $figheight_orig | awk '{if($1>$2){print $1} else{print $1/($2/$1)}}'`
    figheight_orig=`echo $figwidth $figheight_orig | awk '{if($1>$2){print $2} else{print $1}}'`

    gmt grd2cpt $faa -Crainbow > tmp.cpt
    gmt begin ${dataname}/figures/${dataname}_rotation pdf
        xyrange=`gmt gmtinfo ${dataname}_UTM.xyz | awk '{print $5, $6}'`
        xmin=`echo $xyrange | awk -F'[</>]' '{print $2}'`
        xmax=`echo $xyrange | awk -F'[</>]' '{print $3}'`
        ymin=`echo $xyrange | awk -F'[</>]' '{print $5}'`
        ymax=`echo $xyrange | awk -F'[</>]' '{print $6}'`
        xyrange=$xmin/$xmax/$ymin/$ymax
        gmt basemap -JX${figwidth_orig}c/${figheight_orig}c -R$xyrange -Bafg -B+t"Gravity calculation points in UTM CS" --FONT_TITLE=14p,Helvetica,blue
        gmt plot ${dataname}_UTM.xyz -Sc0.05c -Gblack  -Ctmp.cpt  -t20
        plotControlTransformFault tmp_ROTF_UTM.txt
        awk 'NR==1{print $1, $2, $3}' ${dataname}_OTF_UTM.txt | gmt plot -Sk@bullseye -Gred
        awk 'END{print $1, $2, $3}' ${dataname}_OTF_UTM.txt | gmt plot -Sk@bullseye -Gblack
        gmt plot ${dataname}_OTF_UTM.txt -W1p,blue
        
        
        # transformed result
        xy_extend_axes=10000
        xyrange=`gmt gmtinfo ${dataname}_UTM_rot_translate.xyz | awk '{print $5, $6}'`
        xmin=`echo $xyrange $xy_extend_axes | awk -F'[</>]' '{print $2-$7}'`
        xmax=`echo $xyrange $xy_extend_axes | awk -F'[</>]' '{print $3+$7}'`
        ymin=`echo $xyrange $xy_extend_axes | awk -F'[</>]' '{print $5-$7}'`
        ymax=`echo $xyrange $xy_extend_axes | awk -F'[</>]' '{print $6+$7}'`
        xyrange=$xmin/$xmax/$ymin/$ymax

        figheight_trans=`echo $xmin $xmax $ymin $ymax ${figwidth} | awk '{print ($4-$3)/($2-$1)*$5}'`
        figwidth_trans=`echo $figwidth $figheight_trans | awk '{if($1>$2){print $1} else{print $1/($2/$1)}}'`
        figheight_trans=`echo $figwidth $figheight_trans | awk '{if($1>$2){print $2} else{print $1}}'`

        move_x=`echo ${figwidth_orig} | awk '{print $1+1}'`
        move_y=`echo $figheight_orig $figheight_trans | awk '{print ($1-$2)/2}'`
        gmt basemap -JX${figwidth_trans}c/${figheight_trans}c -R$xyrange -X${move_x}c -Y${move_y}c -Bafg -BwSEn+t"Gravity calculation points in UTM CS" --FONT_TITLE=14p,Helvetica,red
        gmt plot ${dataname}_UTM_rot_translate.xyz -Sc0.05c -Gblack  -Ctmp.cpt -t80
        plotControlTransformFault tmp_ROTF_UTM_rot_translate.txt
        awk 'END{print $1, $2, $3}' ${dataname}_OTF_UTM_translate.txt | gmt plot  -Sk@bullseye -Gblack
        gmt plot ${dataname}_OTF_UTM_translate.txt -W1p,blue,-
        awk 'NR==1{print $1, $2, $3}' ${dataname}_OTF_UTM_rot_translate.txt | gmt plot -Sk@bullseye -Gred
        awk 'END{print $1, $2, $3}' ${dataname}_OTF_UTM_rot_translate.txt | gmt plot -Sk@bullseye -Gblack
        gmt plot ${dataname}_OTF_UTM_rot_translate.txt -W1p,red
        echo "Transform fault length: "$length_transform" km,  Rotation angle: "${angle_rot}" degree" | gmt text -F+cTL+f12 -Dj0.1c/0.1c -Gwhite@20 
        echo 0 0 1c $angle_rot 0 | gmt plot -Sm0.4c -W1p,- -Gblack
        # plot length of OTF and ridge, angle between ridge and OTF
        plotLengthAngle_ROTF tmp_ROTF_UTM_rot_translate.txt
        # save ROTF coordinates to input folder
        awk '{print $1, $2}' tmp_ROTF_UTM_rot_translate.txt > ${dataname}/input/${dataname}_RTI_ASPECT.txt
    gmt end show
    rm ${dataname}_OTF_*.txt *.xyz tmp*
}
function makecpt_grd()
{
    grdfile=$1
    data_min=`gmt grdinfo $grdfile | grep "v_min" | awk '{printf "%.1f", $3}'`
    data_max=`gmt grdinfo $grdfile | grep "v_min" | awk '{printf "%.1f", $5}'`
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
    basecpt=$2
    data_min=`gmt grdinfo $grdfile | grep "v_min" | awk '{printf "%.1f", $3}'`
    data_max=`gmt grdinfo $grdfile | grep "v_min" | awk '{printf "%.1f", $5}'`
    cpt_min=`echo ${data_min} | awk '{printf "%.1f", $1*1.5}'`
    cpt_max=`echo ${data_max} | awk '{printf "%.1f", $1*1.5}'`
    if (( $(echo "$data_min > 0" |bc -l) )); then
        gmt makecpt -C${basecpt} -T${data_min}/${data_max}
    else
        gmt makecpt -C${basecpt}+h -T${cpt_min}/${cpt_max}
    fi
}
function gravity_xyz2grd()
{
    gravityFile=$1
    gmt grd2xyz ${faa} >${dataname}.xyz
    paste ${dataname}.xyz ${gravityFile}.txt  | awk '{print $1,$2,$4}' >${gravityFile}.lonlat
    gmt xyz2grd ${gravityFile}.lonlat -R$faa -G${gravityFile}.nc
    rm ${gravityFile}.lonlat ${dataname}.xyz
}
# shift gravity to make gravity along OTF close to zero
function shiftDataTo0OTF()
{
    grd_data=$1
    start_otf=`awk 'NR==1{print $1"/"$2}' ${OTF}`
    stop_otf=`awk 'END{print $1"/"$2}' ${OTF}`
    gmt grdtrack -G$grd_data -E$start_otf/$stop_otf+i1k+d >data_OTF.txt
    meanValue_OTF=`awk -v sum_data=0 -v num_data=0 '{ while ( getline == 1 ) { sum_data+=$4; num_data+=1; }printf "%.0f", sum_data/num_data} ' data_OTF.txt`
    gmt grdmath $grd_data $meanValue_OTF SUB = $grd_data
    str_shiftData=", shift $meanValue_OTF mGal"
}
function FAA2MBA_FFT()
{
    g_wc=${dataname}/Results_FFT/${dataname}_wc_g.nc
    bga=${dataname}/Results_FFT/${dataname}_bga.nc
    g_cm=${dataname}/Results_FFT/${dataname}_cm_g.nc
    mba=${dataname}/Results_FFT/${dataname}_mba.nc
    rmba=${dataname}/Results_FFT/${dataname}_rmba.nc
    moho=${dataname}/Results_FFT/${dataname}_moho.nc
    # 1. Water/crust interface
    gmt gravfft $bathy_large -D$(($rho_crust - $rho_water)) -W$W -Ff -fg -G$g_wc -E4
    if [ ! -f $faa ]; then 
        echo "$faa doesn't exist, try to sample from global data, please check and run again"
        gmt grdsample $globalDataPath/grav/grav.23.nc -R$bathy -I1m -G$faa
        exit
    fi
    # 3. calculate gravity of crust/mantle interface
    gmt gravfft ${bathy_large}=+o-${crust_thickness} -D$(($rho_mantle - $rho_crust)) -fg -G$g_cm -E4
    # =======!!!! Let's calculate MBA in one step ==========
    # 4. mantle bouguer anomaly (MBA)
    # make sure the grid size same as faa
    gmt grdsample $g_wc -R$faa -G$g_wc -V0
    gmt grdsample $g_cm -R$faa -G$g_cm -V0
    gmt grdmath $faa $g_wc SUB $g_cm SUB = $mba
    
    # plot
    cpt_bathy=bathy.cpt
    gmt grd2cpt $bathy -Cbasecpt_bathy.cpt -Z >$cpt_bathy
    gmt begin ${dataname}/figures/FFT/${dataname}_${dataSource}_FFT_W${W}${thermalModelType} ${figfmt}
        gmt subplot begin 2x4 -A+JTL+o0.5c -Fs15c/9c -M0.5c/1.5c -R$bathy -JM15c -Ba -BWSne -T"${dataname}: FFT approach, -W$W"
            gmt subplot set 0,0
                grid_mask=tmp_mask.nc
                alpha_mask=50
                gmt grdmath $bathy_ship $bathy_ship NAN = NAN.nc 
                gmt grdmath $bathy_ship NAN.nc XOR -99999 ADD = $grid_mask
                gmt grdimage $bathy_large -BWseN -Cbathy.cpt
                if [ ! -f ${bathy}.grad ]; then 
                    gmt grdgradient $bathy -A30 -Nt1 -Qc -G${bathy}.grad
                fi
                gmt grdimage $bathy -BWseN -C$cpt_bathy -I${bathy}.grad -t50
                # gmt grdcontour $bathy -C500
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Bathymetry" -By+l"m" -C$cpt_bathy
                plotControlTransformFault
            gmt subplot set 0,1
                gmt grdgradient $faa -A30 -Nt0.6 -Qc -G${faa}.grad
                makecpt_grd $faa
                gmt grdimage $faa -BwseN -I${faa}.grad
                gmt grdcontour $faa -C10
                gmt grdcontour $faa -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                if [ -f $grid_mask ]; then 
                    gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                fi 
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"FAA" -By+l"mGal" -G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 0,2
                gmt grdgradient $g_wc -A30 -Nt0.6 -Qc -G${g_wc}.grad
                makecpt_grd $g_wc 
                gmt grdimage $g_wc -BwseN -I${g_wc}.grad
                gmt grdcontour $g_wc -C10
                gmt grdcontour $g_wc -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity of water/crust interface" -By+l"mGal" -G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 0,3
                gmt grdgradient $g_cm -A30 -Nt0.6 -Qc -G${g_cm}.grad
                makecpt_grd $g_cm 
                gmt grdimage $g_cm -BwsEN -I${g_cm}.grad
                gmt grdcontour $g_cm -C10
                gmt grdcontour $g_cm -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                if [ -f $grid_mask ]; then 
                    gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                fi 
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity of crust/mantle interface" -By+l"mGal" # -G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 1,0
                shiftDataTo0OTF $mba
                gmt grdgradient $mba -A30 -Nt0.6 -Qc -G${mba}.grad
                makecpt_grd $mba
                gmt grdimage $mba -BWseN -I${mba}.grad
                gmt grdcontour $mba -C10
                gmt grdcontour $mba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                if [ -f $grid_mask ]; then 
                    gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                fi 
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"MBA $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 1,1
                gmt basemap -Ba -BwseN
                modelname=Visco-plastic${thermalModelType}
                gravfile=${dataname}/grav_Thermal/grav_${modelname}
                if [ -f ${gravfile}.txt ]; then 
                    gravity_xyz2grd $gravfile
                    grav_therm=${dataname}/grav_Thermal/grav_${modelname}.nc
                    grav_therm_minusMeanValue=grav_therm.nc
                    # mean value of the thermal gravity
                    meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
                    gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                    gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
                    makecpt_grd $grav_therm_minusMeanValue
                    gmt grdimage $grav_therm_minusMeanValue -I${grav_therm_minusMeanValue}.grad
                    gmt grdcontour $grav_therm_minusMeanValue -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                    if [ -f $grid_mask ]; then 
                        gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                    fi 
                    gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity of thermal model: $modelname (minus mean value $meanGrav)" -By+l"mGal" -G${data_min}/${data_max}
                else
                    echo ${gravfile}.txt "does not exist"
                fi
                plotControlTransformFault
            gmt subplot set 1,2
                gmt basemap -BwseN -Ba 
                if [ -f ${gravfile}.txt ]; then 
                    # make the grid size same as mba
                    gmt grdmath $mba $grav_therm_minusMeanValue SUB = $rmba
                    shiftDataTo0OTF $rmba
                    gmt grdgradient $rmba -A30 -Nt0.6 -Qc -G${rmba}.grad
                    # get min max of gravity
                    makecpt_grd $rmba
                    gmt grdimage $rmba -I${rmba}.grad
                    gmt grdcontour $rmba -C10
                    gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                    if [ -f $grid_mask ]; then 
                        gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                    fi 
                    gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}
                fi
                plotControlTransformFault
            gmt subplot set 1,3
                gmt basemap -Ba  -BwsEN
                if [ -f ${moho}.a ]; then 
                    grid=tmp_relative_moho.nc
                    gmt grdmath $moho $crust_thickness ADD -1E-3 MUL = $grid
                    # gmt grdgradient $grid -A30 -Nt0.6 -Qc -G${grid}.grad
                    # cut values close to boundaries
                    range_cut=`gmt_get_gridregion $grid | awk -F '[/]' -v dlondlat=0.05 '{print $1+dlondlat"/"$2-dlondlat"/"$3+dlondlat"/"$4-dlondlat}'`
                    gmt grdcut $grid -R$range_cut -Gtmp_cut_bd.nc
                    makecpt_grd_basecpt tmp_cut_bd.nc $basecpt_moho
                    gmt grdimage tmp_cut_bd.nc -BwseN #-I${grid}.grad
                    gmt grdcontour tmp_cut_bd.nc -C0.5
                    if [ -f $grid_mask ]; then 
                        gmt grdimage $grid_mask -C$cpt_bathy -Q -t$alpha_mask
                    fi 
                    gmt colorbar -DJCB+o0/0.5c -Bxa1f0.5+l"Variation of crust thickness" -By+l"km"  -G${data_min}/${data_max}
                elif [ -f ${dataname}/slice/slice_stressxx_top.xyz ]; then 
                    gmt makecpt -Cbasecpt_stress.cpt+h -T-100/100 -Z
                    gmt contour ${dataname}/slice/slice_stressxx_top.xyz -C -I 
                    # gmt contour ${csv_profile}.xyz -C$clevels -Wthin  -A${clevels}+f7p,Helvetica,black
                    gmt colorbar -DJCB+o0/0.5c -Bxaf+l"stress xx" -By+l"MPa"
                fi
                plotControlTransformFault
        gmt subplot end
    gmt end show
    rm *.nc data_OTF.txt *.grad *.history ${dataname}/input/*.grad ${dataname}/Results_FFT/*.grad *.grd
}
function FAA2MBA_Spatial()
{
    # convert forward gravity xyz to grd
    gravity_xyz2grd ${dataname}/Results_Spatial/${dataname}_grav_waterLayer_$dataSource
    gravity_xyz2grd ${dataname}/Results_Spatial/${dataname}_grav_crustLayer_$dataSource
    
    g_waterLayer=${dataname}/Results_Spatial/${dataname}_grav_waterLayer_$dataSource.nc  #density is rho_m - rho_w
    g_crustLayer=${dataname}/Results_Spatial/${dataname}_grav_crustLayer_$dataSource.nc  #density if rho_m -rho_c
    mba=${dataname}/Results_Spatial/${dataname}_mba.nc
    rmba=${dataname}/Results_Spatial/${dataname}_rmba.nc
    # 4. mantle bouguer anomaly (MBA)
    gmt grdmath $faa $g_waterLayer SUB $g_crustLayer SUB $gconst ADD = $mba
    # plot
    gmt begin ${dataname}/figures/Spatial/${dataname}_${dataSource}_Spatial${thermalModelType} ${figfmt}
        gmt subplot begin 2x4 -A+JTL+o0.5c -Fs15c/9c -M0.5c/1.5c -R$bathy -JM15c -Ba1df30m -BWSne -T"${dataname}: Spatial domain approach"
            gmt subplot set 0,0
                gmt grdgradient $bathy_large -A30 -Nt0.6 -Qc -G${bathy_large}.grad
                gmt grdimage $bathy_large -BWseN -Csealand -I${bathy_large}.grad
                gmt grdcontour $bathy_large -C500
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Bathymetry" -By+l"m"
            gmt subplot set 0,1
                gmt grdgradient $faa -A30 -Nt0.6 -Qc -G${faa}.grad
                makecpt_grd $faa
                gmt grdimage $faa -BwseN -I${faa}.grad
                gmt grdcontour $faa -C10
                gmt grdcontour $faa -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"FAA" -By+l"mGal" -G${data_min}/${data_max}  -G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 0,2
                gmt grdgradient $g_waterLayer -A30 -Nt0.6 -Qc -G${g_waterLayer}.grad
                gmt grdimage $g_waterLayer -BwsEN -I${g_waterLayer}.grad
                gmt grdcontour $g_waterLayer -C10
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity of water layer with density of rho_w - rho_m" -By+l"mGal"
                plotControlTransformFault
            gmt subplot set 0,3
                gmt grdgradient $g_crustLayer -A30 -Nt0.6 -Qc -G${g_crustLayer}.grad
                gmt grdimage $g_crustLayer -BwsEN -I${g_crustLayer}.grad
                gmt grdcontour $g_crustLayer -C10
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity of crust layer with density of rho_c - rho_m" -By+l"mGal"
                plotControlTransformFault
            gmt subplot set 1,0
                shiftDataTo0OTF $mba
                gmt grdgradient $mba -A30 -Nt0.6 -Qc -G${mba}.grad
                makecpt_grd $mba
                gmt grdimage $mba -BWseN -I${mba}.grad
                gmt grdcontour $mba -C10
                gmt grdcontour $mba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"MBA $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 1,1
                gmt basemap -Ba -BwseN
                modelname=Viscous-plastic${thermalModelType}
                gravfile=${dataname}/grav_Thermal/grav_${modelname}
                if [ -f ${gravfile}.txt ]; then 
                    gravity_xyz2grd $gravfile
                    grav_therm=${dataname}/grav_Thermal/grav_${modelname}.nc
                    grav_therm_minusMeanValue=grav_therm.nc
                    # mean value of the thermal gravity
                    meanGrav=`grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
                    gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                    gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
                    makecpt_grd $grav_therm_minusMeanValue
                    gmt grdimage $grav_therm_minusMeanValue -I${grav_therm_minusMeanValue}.grad
                    gmt grdcontour $grav_therm_minusMeanValue -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                    gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity of thermal model: $modelname (minus mean value $meanGrav)" -By+l"mGal" -G${data_min}/${data_max}
                fi
                plotControlTransformFault
            gmt subplot set 1,2
                gmt basemap -BwseN -Ba 
                if [ -f ${gravfile}.txt ]; then 
                    # make the grid size same as mba
                    gmt grdmath $mba $grav_therm_minusMeanValue SUB = $rmba
                    shiftDataTo0OTF $rmba
                    gmt grdgradient $rmba -A30 -Nt0.6 -Qc -G${rmba}.grad
                    # get min max of gravity
                    makecpt_grd $rmba
                    gmt grdimage $rmba -I${rmba}.grad
                    gmt grdcontour $rmba -C10
                    gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                    gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}
                fi
                plotControlTransformFault
            gmt subplot set 1,3
                gmt basemap -Ba  -BWseN
                modelname=Viscous${thermalModelType}
                gravfile=${dataname}/grav_Thermal/grav_${modelname}
                if [ -f ${gravfile}.txt ]; then 
                    gravity_xyz2grd ${dataname}/grav_Thermal/grav_$modelname
                    grav_therm=${dataname}/grav_Thermal/grav_${modelname}.nc
                    grav_therm_minusMeanValue=grav_therm.nc
                    # mean value of the thermal gravity
                    meanGrav=`grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{print $3}' && rm tmp.nc`
                    gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                    gmt grdmath $mba $grav_therm_minusMeanValue SUB = $rmba
                    shiftDataTo0OTF $rmba
                    gmt grdgradient $rmba -A30 -Nt0.6 -Qc -G${rmba}.grad
                    # get min max of gravity
                    makecpt_grd $rmba
                    gmt grdimage $rmba -I${rmba}.grad 
                    gmt grdcontour $rmba -C10
                    gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                    gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA ($modelname) $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}
                fi
                plotControlTransformFault
        gmt subplot end
    gmt end show
    rm *.nc data_OTF.txt *.grad *.history ${dataname}/input/*.grad ${dataname}/Results_Spatial/*.grad 
}
function diff_Gravity_models()
{
    model1=Viscous
    model2=Viscous-plastic
    model3=Viscous-elastic-plastic
    gmt begin ../figures/${dataname}_grav_models ${figfmt}
        gmt subplot begin 2x3 -A+JTL+o0.5c -Fs15c/9c -M0.5c/1.5c -R$bathy -JM15c -Ba1df30m -BWSne -T"Gravity of thermal model: "$dataname
            gmt subplot set 0,0
                modelname=$model1
                gravity_xyz2grd grav_Thermal/grav_$modelname
                grav_therm=grav_Thermal/grav_${modelname}.nc
                grav_therm_minusMeanValue=grav_therm.nc
                # mean value of the thermal gravity
                meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
                gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
                gmt grdimage $grav_therm_minusMeanValue -BWseN -I${grav_therm_minusMeanValue}.grad
                gmt grdcontour $grav_therm_minusMeanValue -C10
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal"
                plotControlTransformFault
            gmt subplot set 0,1
                modelname=$model2
                gravity_xyz2grd grav_Thermal/grav_$modelname
                grav_therm=grav_Thermal/grav_${modelname}.nc
                grav_therm_minusMeanValue=grav_therm_${modelname}.nc
                # mean value of the thermal gravity
                meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
                gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
                gmt grdimage $grav_therm_minusMeanValue -BwseN -I${grav_therm_minusMeanValue}.grad
                gmt grdcontour $grav_therm_minusMeanValue -C10
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal"
                plotControlTransformFault
            gmt subplot set 0,2
                modelname=$model3
                gravity_xyz2grd grav_Thermal/grav_$modelname
                grav_therm=grav_Thermal/grav_${modelname}.nc
                grav_therm_minusMeanValue=grav_therm_${modelname}.nc
                # mean value of the thermal gravity
                meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
                gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
                gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
                gmt grdimage $grav_therm_minusMeanValue -BwsEN -I${grav_therm_minusMeanValue}.grad
                gmt grdcontour $grav_therm_minusMeanValue -C10
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$modelname (minus mean value $meanGrav mGal)" -By+l"mGal"
                plotControlTransformFault
            gmt subplot set 1,0
                m1=$model1
                m2=$model2
                grav_therm1=grav_Thermal/grav_${m1}.nc
                grav_therm2=grav_Thermal/grav_${m2}.nc
                grav_diff=grav_diff.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BWseN -I${grav_diff}.grad -Cgrav_diff.cpt
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity difference: $m1 - $m2" -By+l"mGal"  -Cgrav_diff.cpt
                plotControlTransformFault
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
            gmt subplot set 1,1
                m1=$model1
                m2=$model3
                grav_therm1=grav_Thermal/grav_${m1}.nc
                grav_therm2=grav_Thermal/grav_${m2}.nc
                grav_diff=grav_diff.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BwseN -I${grav_diff}.grad  -Cgrav_diff.cpt
                # gmt grdcontour $grav_diff -C1
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity difference: $m1 - $m2" -By+l"mGal"  -Cgrav_diff.cpt
                plotControlTransformFault
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
            gmt subplot set 1,2
                m1=$model2
                m2=$model3
                grav_therm1=grav_Thermal/grav_${m1}.nc
                grav_therm2=grav_Thermal/grav_${m2}.nc
                grav_diff=grav_diff.nc
                gmt grdmath $grav_therm1 $grav_therm2 SUB = ${grav_diff}
                gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
                gmt grdimage $grav_diff -BwsEN -I${grav_diff}.grad  -Cgrav_diff.cpt
                # gmt grdcontour $grav_diff -C0, -A0,+f10p,red -Wathickest,black -Wcthinnest,blue,-
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity difference: $m1 - $m2" -By+l"mGal"  -Cgrav_diff.cpt
                plotControlTransformFault
                zmin=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $3}'`
                zmax=`gmt grdinfo $grav_diff | grep "v_min" | awk '{printf "%.1f", $5}'`
                echo "min: $zmin mGal, max: $zmax mGal" | gmt text -F+cTL+f16,black -Dj1c/1c -Gwhite@20
        gmt subplot end
    gmt end show
}
function mergeShipSatelliteBathy()
{
    dxdy=$1
    # merge ship-based bathymetry and satellite bathymetry
    gmt grd2xyz $bathy_ship | gmt convert -bo > ship.b
    gmt nearneighbor -R$range_large -I$dxdy -S5k -Gship.nc ship.b -bi
    gmt grdsample $ETOPO1 -R$range_large -I$dxdy -Gtmp_etopo.nc
    gmt grdmath ship.nc tmp_etopo.nc AND = ${dataname}/input/${dataname}_large_ship.nc 
    gmt grdsample $ETOPO1 -R$range_large -I$dxdy -G${dataname}/input/${dataname}_large_Sat.nc
    # gmt grdcut  $ETOPO1 -R$range_large  -G${dataname}/input/${dataname}_large_Sat.nc
    rm ship.nc ship.b tmp*.nc
}
function addColors2SegmentLines()
{
    segmentFile=$1
    colors=$2
    index_segment=0
    awk -v index_segment=0 -v colors="${colors[*]}" '
        {
            if ( ($1 == ">") )
            {
                index_segment=index_segment+1;
                split(colors,cArray," ");
                print "> -W2p,"cArray[index_segment]
            }
            else{ print $0 }
        }
    ' $segmentFile >tmp
    mv tmp $segmentFile
}
function plotProfiles()
{
    data=$1 
    label=$2
    method=$3
    title=$4
    figheight=`gmt_get_map_height -R$data  -JM${figwidth}c`
    # plot
    gmt begin ${dataname}/figures/${method}/${dataname}_${dataSource}_Profiles_${label}_${method} ${figfmt}
        gmt subplot begin 1x4 -A+JTL+o0.5c -Fs${figwidth}c/${figheight}c -M0.5c/1.5c -R$data -JM${figwidth}c -Ba1df30m -BWSne -T"${title}"
            gmt subplot set 0,0
                gmt grdgradient $data -A30 -Nt0.6 -Qc -G${mba}.grad
                makecpt_grd $data
                gmt grdimage $data -BWseN -I${mba}.grad -t50
                gmt grdcontour $data -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"$label" -By+l"mGal" -G${data_min}/${data_max}
                plotControlTransformFault
                # # track profiles perpendicular to the OTF
                awk '{print $1, $2}' ${OTF} | gmt grdtrack  -G$data -C60k/0.2k/20k+v > prof_OTF.txt
                colors=(green@100 red blue cyan purple pink cyan)
                addColors2SegmentLines prof_OTF.txt $colors
                gmt plot -W1p,+ve0.5c+h0.5 prof_OTF.txt
                awk '{print $1, $2}' ${dataname}/input/${dataname}_ridgeN.txt | gmt grdtrack  -G$data -C60k/0.2k/30k+v > prof_ridgeN.txt
                addColors2SegmentLines prof_ridgeN.txt $colors
                gmt plot -W1p,+ve0.5c+h0.5 prof_ridgeN.txt
                awk '{print $1, $2}' ${dataname}/input/${dataname}_ridgeS.txt | gmt grdtrack  -G$data -C100k/0.2k/30k+v > prof_ridgeS.txt
                addColors2SegmentLines prof_ridgeS.txt $colors
                gmt plot -W1p,+ve0.5c+h0.5 prof_ridgeS.txt
            gmt subplot set 0,1
                gmt basemap -R-30/30/${data_min}/${data_max} -JX$figwidth/$figheight -Bafg -By+l"MBA (mGal)" -Bx+l"Distance from OTF (km)" -BWSen+t"OTF"
                echo "0 -999" >tmp && echo "0 999" >>tmp && cat tmp | gmt plot -W3p,black,-
                gmt plot prof_OTF.txt -i2,4 
            gmt subplot set 0,2
                gmt basemap -R-30/30/${data_min}/${data_max} -JX$figwidth/$figheight -Bafg -Bx+l"Distance from ridge (km)" -BwSen+t"Northern ridge"
                echo "0 -999" >tmp && echo "0 999" >>tmp && cat tmp | gmt plot -W3p,black,-
                gmt plot prof_ridgeN.txt -i2,4 -W1p,blue 
            gmt subplot set 0,3
                gmt basemap -R-50/50/${data_min}/${data_max} -JX$figwidth/$figheight -Bafg -By+l"MBA (mGal)" -Bx+l"Distance from ridge (km)" -BwSEn+t"Southern ridge"
                echo "0 -999" >tmp && echo "0 999" >>tmp && cat tmp | gmt plot -W3p,black,-
                gmt plot prof_ridgeS.txt -i2,4 -W1p,blue 
        gmt subplot end
    gmt end show
    rm prof* *.history tmp ${dataname}/Results_${method}/*.grad 
}
function PrepareData()
{
    # 1.1 cut faa
    gmt grdsample $globalDataPath/grav/grav.23.nc -R$bathy -I1m -G$faa
    # 1.2 merge bathymetry to a larger area
    mergeShipSatelliteBathy 100e
}
function RMBA2Moho()
{
    error=$1
    faa=${dataname}/input/${dataname}_faa.nc
    rmba=${dataname}/Results_FFT/${dataname}_rmba.nc
    wc_g=${dataname}/Results_FFT/${dataname}_wc_g.nc
    cm_g=${dataname}/Results_FFT/${dataname}_cm_g.nc
    therm_g=${dataname}/grav_Thermal/grav_Visco-plastic.nc
    if [ ! -f $therm_g ]; then 
        echo "Thermal gravity does not exit: $therm_g"
        exit
    fi
    grav=tmp_bga.nc 
    drho=$(($rho_mantle - $rho_crust))
    # ----------------------
    gmt grdmath $rmba $cm_g ADD = $grav 
    lon0=`gmt grdinfo ${grav} | grep "x_min" | awk '{print $3}'`
    lat0=`gmt grdinfo ${grav} | grep "y_min" | awk '{print $3}'`
    gmt grd2xyz ${grav} | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=0 +y_0=0 >${grav}.xyz
    xlim=`gmt gmtinfo ${grav}.xyz | awk '{print $5}' | cut -d "<" -f 2 | cut -d ">" -f 1`
    ylim=`gmt gmtinfo ${grav}.xyz | awk '{print $6}' | cut -d "<" -f 2 | cut -d ">" -f 1`
    # gdps only support Cartesian grid of surfer 6 format
    gmt grdedit ${grav} -R$xlim/$ylim -G${grav}.grd=sf

    # moho inversion using Parker's method
    gdps -D $drho -H $crust_thickness -f ${grav}.grd -O ${grav}_moho.grd -r $error
    # convert to lonlat CS
    gmt grdedit ${grav}_moho.grd -R$grav -G${dataname}/Results_FFT/${dataname}_moho.nc -fg
    # calculate data fitting
    gmt gravfft ${dataname}/Results_FFT/${dataname}_moho.nc -D$drho -Nf -fg -E4 -G${dataname}/Results_FFT/${dataname}_predict_bga.nc -W0 # -V
    # calculate misfit
    gmt grdmath ${dataname}/Results_FFT/${dataname}_predict_bga.nc $grav SUB = ${dataname}/Results_FFT/${dataname}_err.nc
    rm tmp*
}
function RMBA2Moho_v2()
{
    start_iter=$1
    numIter=$2
    path_tmp=$dataname/Results_FFT/tmp_Moho
    if [ ! -d $path_tmp ]; then 
        mkdir $path_tmp
    fi
    if [ -z $start_iter ]; then
        start_iter=0
    fi
    if [ -z $numIter ]; then 
        numIter=600
    fi
    rmba=$dataname/Results_FFT/${dataname}_rmba.nc
    cm_g=$dataname/Results_FFT/${dataname}_cm_g.nc
    # ------------------------------------
    grav=tmp_bga.nc
    gmt grdmath $rmba $cm_g ADD = $grav  #use bga to inversion
    # cp $rmba $grav   #use rmba directly
    # ------------------------------------
    geo="-fg"
    w=10 # iternation step
    drho=$(($rho_mantle - $rho_crust))
    z0=-$crust_thickness
    lon0=`gmt grdinfo ${grav} | grep "x_min" | awk '{print $3}'`
    lat0=`gmt grdinfo ${grav} | grep "y_min" | awk '{print $3}'`
    gmt grd2xyz ${grav} | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=0 +y_0=0 >${grav}.xyz
    xlim=`gmt gmtinfo ${grav}.xyz | awk '{print $5}' | cut -d "<" -f 2 | cut -d ">" -f 1`
    ylim=`gmt gmtinfo ${grav}.xyz | awk '{print $6}' | cut -d "<" -f 2 | cut -d ">" -f 1`
    range_cart=$xlim/$ylim
    # ===============begin =================================
    h_old=$path_tmp/Moho_${start_iter}.nc
    if [ ! -f $h_old ]; then 
        gmt grdmath $grav 0 MUL $z0 ADD = $h_old
        start_iter=0
    fi
    Norm2_inputdata=`gmt grd2xyz ${grav} | awk '{sum += ($3*$3); N+=1} END {print sqrt(sum)}'`
    fit_smoth=$path_tmp/fit_smoth.txt
    rm $fit_smoth && touch $fit_smoth
    for tau in $(seq ${start_iter} 1 $((${start_iter}+${numIter}))) ;
    do 
        KU=${grav}_predict.nc
        R=tmp_R.nc 
        KR=tmp_KR.nc 
        wKR=tmp_wKR.nc
        # forward KU
        gmt gravfft $h_old -D$drho -Ff -Nf+a  -W0 -G$KU
        gmt grdmath $grav $KU SUB = $R
        # forward KR
        gmt gravfft $R -D$drho -Ff -Nf+a  -W0 -G$KR
        gmt grdmath $KR $w MUL = $wKR
        h_new=$path_tmp/Moho_${tau}.nc
        gmt grdmath $h_old $wKR ADD = $h_new
        h_old=$h_new
        
        # fitting error
        FIT=`gmt grd2xyz ${R} | awk -v norm2_data=$Norm2_inputdata '{sum += ($3*$3); N+=1} END {print sqrt(sum)/norm2_data}'`
        gmt grdedit ${h_new} -R$range_cart -Gh_new_cart.grd=sf
        gmt grdmath h_new_cart.grd DDX = tmp_ddx.nc 
        gmt grdmath h_new_cart.grd DDY = tmp_ddy.nc 
        grad=tmp_grad.nc
        gmt grdmath tmp_ddx.nc tmp_ddy.nc HYPOT = $grad
        SMOTH=`gmt grd2xyz ${grad} | awk '{sum += ($3*$3); N+=1} END {print sqrt(sum)}'`
        # calculate err between the true solution
        if [ ! -z $solution ]; then 
            gmt grdmath $solution $h_new SUB = tmp_err.nc 
            ERR=`gmt grd2xyz tmp_err.nc | awk '{sum += ($3*$3); N+=1} END {print sqrt(sum/N)}'`
        else 
            ERR=0
        fi
        echo $tau $FIT $SMOTH $ERR
        echo $tau $FIT $SMOTH $ERR >>$fit_smoth

    done
    rm tmp* h_new_cart.*
}
function plotMohoInversion()
{
    moho=$1
    if [ ! -f $moho ]; then 
        echo "The moho file does not exist: $moho"
        exit
    fi
    cp $moho ${dataname}/Results_FFT/${dataname}_moho.nc
    rmba=${dataname}/Results_FFT/${dataname}_rmba.nc
    g_cm=${dataname}/Results_FFT/${dataname}_cm_g.nc
    bga_predict=${dataname}/Results_FFT/${dataname}_predict_bga.nc
    err=${dataname}/Results_FFT/${dataname}_err.nc
    gmt gravfft $moho -D$(($rho_mantle - $rho_crust)) -Ff -Nf+a -fg  -W0 -G$bga_predict
    figheight=`gmt_get_map_height -R$rmba  -JM${figwidth}c`
    gmt begin ${dataname}/figures/FFT/${dataname}_${dataSource}_FFT_W${W}${thermalModelType}_Moho ${figfmt}
        gmt subplot begin 2x3 -A+JTL+o0.5c -Fs${figwidth}c/${figheight}c -M0.5c/1.5c -R$bathy -JM${figwidth}c -Ba30mf30m -BWSne -T"${dataname}: FFT approach, -W$W"
            gmt subplot set 0,0
                gmt basemap -BwseN -Ba 
                if [ -f $rmba ]; then 
                    # make the grid size same as mba
                    shiftDataTo0OTF $rmba
                    gmt grdgradient $rmba -A30 -Nt0.6 -Qc -G${rmba}.grad
                    # get min max of gravity
                    makecpt_grd $rmba
                    gmt grdimage $rmba -I${rmba}.grad
                    gmt grdcontour $rmba -C10
                    gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                    gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA" -By+l"mGal" -G${data_min}/${data_max}
                fi
                plotControlTransformFault
            gmt subplot set 0,1
                gmt grdgradient $g_cm -A30 -Nt0.6 -Qc -G${g_cm}.grad
                makecpt_grd $g_cm 
                gmt grdimage $g_cm -BwseN -I${g_cm}.grad
                gmt grdcontour $g_cm -C10
                gmt grdcontour $g_cm -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Gravity of crust/mantle interface" -By+l"mGal"  -G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 0,2
                grid=tmp_bga.nc
                gmt grdmath $rmba $g_cm ADD = $grid
                # calculate err
                gmt grdmath $bga_predict $grid SUB = $err
                gmt grdgradient $grid -A30 -Nt0.6 -Qc -G${grid}.grad
                makecpt_grd $grid 
                gmt grdimage $grid -BwsEN -I${grid}.grad
                gmt grdcontour $grid -C10
                gmt grdcontour $grid -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"BGA (remove thermal effect)" -By+l"mGal" # -G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 1,0
                grid=tmp_relative_moho.nc
                gmt grdmath $moho $crust_thickness ADD 1E-3 MUL = $grid
                # gmt grdgradient $grid -A30 -Nt0.6 -Qc -G${grid}.grad
                # cut values close to boundaries
                range_cut=`gmt_get_gridregion $grid | awk -F '[/]' -v dlondlat=0.1 '{print $1+dlondlat"/"$2-dlondlat"/"$3+dlondlat"/"$4-dlondlat}'`
                gmt grdcut $grid -R$range_cut -Gtmp_cut_bd.nc
                makecpt_grd_basecpt tmp_cut_bd.nc $basecpt_moho
                gmt grdimage $grid -BwseN #-I${grid}.grad
                gmt grdcontour $grid -C0.5
                gmt colorbar -DJCB+o0/0.5c -Bxa1f0.5+l"Moho variation (relative to a constant crust thickness)" -By+l"km"  -G${data_min}/${data_max}
                plotControlTransformFault
            gmt subplot set 1,2
                meanErr=`gmt grdmath $err MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{print $3}' && rm tmp.nc`
                # shift a mean value
                echo "shift mean value of $meanErr"
                gmt grd2xyz $err | awk -v mean=$meanErr '{print $1, $2, $3-mean}' >tmp.xyz
                xigma=`awk '{sum+=$3*$3; n++;} END {printf "%.1f", (sqrt(sum/n))}' tmp.xyz`
                gmt histogram tmp.xyz -Bxa5f1+l"Data fitting error (mGal)" -Bya5f1+l"Frequency"+u"%"  \
                    -BwSEn+glightblue -R-15/18/0/15 -JX${figwidth}c/${figheight}c -Gorange -W1p -T0.5 -i2 -Z1 
                echo "@%12%\163@%% = $xigma mGal" | gmt text -F+cTL+f14 -Dj0.1c/0.1c -Gwhite@20
                
            gmt subplot set 1,1
                grid=$bga_predict
                # shift a mean value
                gmt grdmath $bga_predict $g_cm SUB $meanErr SUB = tmp_rmba_predict.nc 
                grid=tmp_rmba_predict.nc

                gmt grdgradient $grid -A30 -Nt0.6 -Qc -G${grid}.grad
                # makecpt_grd $grid 
                makecpt_grd $rmba
                gmt grdimage $grid -BwseN -I${grid}.grad
                gmt grdcontour $grid -C10
                gmt grdcontour $grid -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                gmt colorbar -DJCB+o0/0.5c -Bxaf+l"Predicted RMBA" -By+l"mGal"  -G${data_min}/${data_max}
                plotControlTransformFault
                rm ${grid}.grad
            
        gmt subplot end
    gmt end show
    rm *.nc data_OTF.txt *.grad *.history ${dataname}/input/*.grad ${dataname}/Results_FFT/*.grad tmp*
}
# convert ASPECT results to lonlat
function xy2lonlat_vtu_ASPECT()
{
    # the input file must be csv format exported from vtu/vtk
    controlPointsFile_lonlat=$1
    translate_x=$2
    translate_y=$3
    csv_profile=$4
    # 0. convert this transform fault line parallel to x axis in ASPECT modeling
    # suggest: set one endpoint of the transform fault as (0,0) in ASPECT modeling coordinate system !!
    lon0=`awk 'NR==1{print $1}' ${OTF}` 
    lat0=`awk 'NR==1{print $2}' ${OTF}` 
    x0=0
    y0=0
    cat $controlPointsFile_lonlat | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 >${dataname}_OTF_UTM.txt
    # awk 'NR==2{print "length of control line: "sqrt($1*$1+$2*$2)}' input/${dataname}_OTF_UTM.txt
    cat $ROTF | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 >tmp_ROTF_UTM.txt
    # calculate cos(theta) and sin(theta)
    costheta=`awk 'NR==2{print $1/sqrt($1*$1+$2*$2)}' ${dataname}_OTF_UTM.txt`
    sintheta=`awk 'NR==2{print $2/sqrt($1*$1+$2*$2)}' ${dataname}_OTF_UTM.txt`
    angle_rot=`echo $sintheta $costheta | awk '{printf "%.1f", atan2($1, $2)/3.141592653*180}'`
    echo "Rotation angle: "${angle_rot}" degree"
    
    # # translate 
    # awk -v dx=$translate_x -v dy=$translate_y '{print $1+dx, $2+dy, $3}' ${dataname}_UTM_rot.xyz >${dataname}_UTM_rot_translate.xyz
    # awk -v dx=$translate_x -v dy=$translate_y '{print $1+dx, $2+dy, $3}' ${dataname}_OTF_UTM_rot.txt >${dataname}_OTF_UTM_rot_translate.txt
    # awk -v dx=$translate_x -v dy=$translate_y '{print $1+dx, $2+dy, $3}' tmp_ROTF_UTM_rot.txt >tmp_ROTF_UTM_rot_translate.txt
    # awk -v dx=$translate_x -v dy=$translate_y '{print $1+dx, $2+dy, $3}' ${dataname}_OTF_UTM.txt >${dataname}_OTF_UTM_translate.txt
    
    # 1. rotate x,y
    awk -v co=$costheta -v si=$sintheta -F ',' 'NR>1{{x=$2; y=$3;data=$1}; print co*x-si*y, si*x+co*y, data}' ${csv_profile}.csv >${csv_profile}_UTM.xyz
    # 2. xy to lonlat
    cat ${csv_profile}_UTM.xyz | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 -I -f "%.8E" | awk '{print $1, $2, $3/1E6}' >${csv_profile}.xyz
    gmt begin ${dataname}/Results_FFT/slice pdf
        T=$(gmt info -T5+c2 ${csv_profile}.xyz)
        T=-T-100/100
        echo $T
        gmt makecpt -Cbasecpt_stress.cpt+h "$T" -Z
        gmt basemap -JM18c -R$faa -Ba
        gmt contour ${csv_profile}.xyz -C -I 
        # gmt contour ${csv_profile}.xyz -C$clevels -Wthin  -A${clevels}+f7p,Helvetica,black
        plotControlTransformFault
        gmt colorbar -DJCB+o0/1c -Bxaf+l"stress xx" -By+l"MPa"
    gmt end show
    rm ${dataname}_OTF_*.txt *.xyz tmp* rm
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
    # start point of TF with coordinate of (0,0)
    # center point of OTF
    lonC_otf=`echo $lon_start_otf $lon_end_otf | awk '{print ($1+$2)/2.0}'`
    latC_otf=`echo $lat_start_otf $lat_end_otf | awk '{print ($1+$2)/2.0}'`
    xCyC_otf=`echo $lonC_otf $latC_otf | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 | awk '{print $1, $2}'`
    XmaxYmax_otf=`echo $lon_end_otf $lat_end_otf | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 | awk '{print $1, $2}'`
    length_otf=`echo $XmaxYmax_otf | awk '{print sqrt($1*$1+$2*$2)}'`
    length_box=`echo $length_otf $length_frac_box | awk '{print $1*$2}'`
    # echo "length_box: "$length_box
    costheta=`echo $XmaxYmax_otf | awk '{printf "%.8f", $1/sqrt($1*$1+$2*$2)}'`
    sintheta=`echo $XmaxYmax_otf | awk '{printf "%.8f",$2/sqrt($1*$1+$2*$2)}'`
    angle_rot=`echo $sintheta $costheta | awk '{printf "%.1f", atan2($1, $2)/3.141592653*180}'`
    # set final center xy: local coordinate on OTF is [-0.5, 0.5]
    xCyC_box=`echo $xCyC_otf $costheta $sintheta $length_otf $loc_frac_box $length_box | awk '{print $1+$6*$5*$3, $2+$5*$6*$4}'`
    lonClatC_box=`echo $xCyC_box | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}'`
    # xy coordinate of OTF start point 
    xy_start_box=`echo $xCyC_box $costheta $sintheta $length_box | awk '{print $1-$5/2.0*$3, $2-$5/2.0*$4}'`
    # xy coordinate of OTF end point 
    xy_end_box=`echo $xCyC_box $costheta $sintheta $length_box | awk '{print $1+$5/2.0*$3, $2+$5/2.0*$4}'`
    
    # xy coordinate of averate box nodes
    lonlat_start_box=`echo $xy_start_box | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}'`
    echo $xy_start_box $costheta $sintheta $halfWidth_box | awk '{print $1-$5*$4, $2+$5*$3}' >${dataname}/input/${boxname}.xy
    echo $xy_start_box $costheta $sintheta $halfWidth_box | awk '{print $1+$5*$4, $2-$5*$3}' >>${dataname}/input/${boxname}.xy
    echo $xy_end_box $costheta $sintheta $halfWidth_box | awk '{print $1+$5*$4, $2-$5*$3}' >>${dataname}/input/${boxname}.xy
    echo $xy_end_box $costheta $sintheta $halfWidth_box | awk '{print $1-$5*$4, $2+$5*$3}' >>${dataname}/input/${boxname}.xy
    # additional rotation with rotate center one of RTIs: if loc_frac_box<0, the first RTI; if loc_frac_box>0 the second RTI
    row_RTI=`echo $ind_OTF | awk '{print $1*2}'`
    if [ `echo "$loc_frac_box < 0" | bc` -eq 1 ]; then 
        row_RTI=`echo $ind_OTF | awk '{print $1*2-1}'`
    fi 
    x0y0=(`awk -v row=$row_RTI 'NR==row{print $1, $2}' $OTF | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 | awk '{print $1, $2}'`)
    costheta=`echo $additional_roation | awk '{printf "%.6f", cos(($1)/180*3.141592653)}'`
    sintheta=`echo $additional_roation | awk '{printf "%.6f", sin(($1)/180*3.141592653)}'`
    awk -v co=$costheta -v si=$sintheta -v x0=${x0y0[0]} -v y0=${x0y0[1]} '{print co*($1-x0)+si*($2-y0) +x0, -si*($1-x0)+co*($2-y0) + y0, $3}' ${dataname}/input/${boxname}.xy >tmp.xy
    mv tmp.xy ${dataname}/input/${boxname}.xy
    # convert box coordinate (x,y) to (lon,lat)
    cat ${dataname}/input/${boxname}.xy | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}' >${dataname}/input/${boxname}.lonlat
    # lonClatC_box
    lonClatC_box=`awk '{{meanlon+=$1; meanlat+=$2}; if(NR==4){print meanlon/4, meanlat/4}}' ${dataname}/input/${boxname}.lonlat`
    angle_rot=`echo "$angle_rot - $additional_roation" | bc `
}
function getAverageBox_Ridge()
{
    boxname=$1
    ind_Ridge=$2
    loc_frac_box=$3 #-1.2 #-1 means left FZ; +1 means right FZ; 0 means OTF
    length_frac_box=$4 # 0.8 # percentage of TF length
    w_box=$5
    # echo $boxname $ind_Ridge $loc_frac_box $length_frac_box $w_box
    halfWidth_box=`echo $w_box | awk '{print $1/2*1000}'` #m
    additional_roation=$6
    file_ridge_index=${dataname}/input/ridge_index.txt
    file_ridge_otf=${dataname}/input/${dataname}_ridge_otf.txt
    index_start_ridge=`awk -v indRidge=$ind_Ridge 'NR==indRidge{print $1}' $file_ridge_index`
    index_end_ridge=`awk -v indRidge=$ind_Ridge 'NR==indRidge{print $2}' $file_ridge_index`
    # echo $index_start_ridge $index_end_ridge
    lon_start_ridge=`awk -v row=$index_start_ridge 'NR==row{print $1}' $file_ridge_otf`
    lon_end_ridge=`awk -v row=$index_end_ridge 'NR==row{print $1}' $file_ridge_otf`
    lat_start_ridge=`awk -v row=$index_start_ridge 'NR==row{print $2}' $file_ridge_otf`
    lat_end_ridge=`awk -v row=$index_end_ridge 'NR==row{print $2}' $file_ridge_otf`
    # echo $lon_start_ridge $lon_end_ridge
    # start point of Ridge with coordinate of (0,0)
    # center point of Ridge
    lonC_ridge=`echo $lon_start_ridge $lon_end_ridge | awk '{print ($1+$2)/2.0}'`
    latC_ridge=`echo $lat_start_ridge $lat_end_ridge | awk '{print ($1+$2)/2.0}'`
    xCyC_ridge=`echo $lonC_ridge $latC_ridge | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_ridge +lat_0=$lat_start_ridge +x_0=0 +y_0=0 | awk '{print $1, $2}'`
    XmaxYmax_ridge=`echo $lon_end_ridge $lat_end_ridge | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_ridge +lat_0=$lat_start_ridge +x_0=0 +y_0=0 | awk '{print $1, $2}'`
    length_ridge=`echo $XmaxYmax_ridge | awk '{print sqrt($1*$1+$2*$2)}'`
    length_box=`echo $length_ridge $length_frac_box | awk '{print $1*$2}'`
    # echo "length_box: "$length_box
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
    echo $xy_start_box $costheta $sintheta $halfWidth_box | awk '{print $1-$5*$4, $2+$5*$3}' >${dataname}/input/${boxname}.xy
    echo $xy_start_box $costheta $sintheta $halfWidth_box | awk '{print $1+$5*$4, $2-$5*$3}' >>${dataname}/input/${boxname}.xy
    echo $xy_end_box $costheta $sintheta $halfWidth_box | awk '{print $1+$5*$4, $2-$5*$3}' >>${dataname}/input/${boxname}.xy
    echo $xy_end_box $costheta $sintheta $halfWidth_box | awk '{print $1-$5*$4, $2+$5*$3}' >>${dataname}/input/${boxname}.xy
    # additional rotation with rotate center one of RTIs: if loc_frac_box<0, the first RTI; if loc_frac_box>0 the second RTI
    # row_RTI=`echo $ind_OTF | awk '{print $1*2}'`
    # if [ `echo "$loc_frac_box < 0" | bc` -eq 1 ]; then 
    #     row_RTI=`echo $ind_OTF | awk '{print $1*2-1}'`
    # fi 
    # x0y0=(`awk -v row=$row_RTI 'NR==row{print $1, $2}' $OTF | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_otf +lat_0=$lat_start_otf +x_0=0 +y_0=0 | awk '{print $1, $2}'`)
    # costheta=`echo $additional_roation | awk '{printf "%.6f", cos(($1)/180*3.141592653)}'`
    # sintheta=`echo $additional_roation | awk '{printf "%.6f", sin(($1)/180*3.141592653)}'`
    # awk -v co=$costheta -v si=$sintheta -v x0=${x0y0[0]} -v y0=${x0y0[1]} '{print co*($1-x0)+si*($2-y0) +x0, -si*($1-x0)+co*($2-y0) + y0, $3}' ${dataname}/input/${boxname}.xy >tmp.xy
    # mv tmp.xy ${dataname}/input/${boxname}.xy
    # convert box coordinate (x,y) to (lon,lat)
    cat ${dataname}/input/${boxname}.xy | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_start_ridge +lat_0=$lat_start_ridge +x_0=0 +y_0=0 -I -f "%.8f" | awk '{print $1, $2}' >${dataname}/input/${boxname}.lonlat
    # lonClatC_box
    lonClatC_box=`awk '{{meanlon+=$1; meanlat+=$2}; if(NR==4){print meanlon/4, meanlat/4}}' ${dataname}/input/${boxname}.lonlat`
    angle_rot=`echo "$angle_rot - $additional_roation" | bc `
}
function calAverageRMBA()
{
    rmba=${dataname}/Results_FFT/${dataname}_rmba.nc
    echo $figwidth
    figheight=`gmt_get_map_height -R$rmba  -JM${figwidth}c`
    gmt begin ${dataname}/figures/FFT/${dataname}_dRMBA ${figfmt}
        gmt subplot begin 1x2 -A+JTL+o0.5c -Fs${figwidth}c/${figheight}c -M0.5c/1.5c -R$rmba -JM${figwidth}c -Ba30mf30m -BWSne #-T"${dataname}: FFT approach, -W$W"
            gmt subplot set 0,0
                grid_mask=tmp_mask.nc
                alpha_mask=50
                gmt grdmath $bathy_ship $bathy_ship NAN = NAN.nc 
                gmt grdmath $bathy_ship NAN.nc XOR 99999 ADD = $grid_mask
                gmt basemap -JM${figwidth}c -R$rmba -BWSen -Bafg
                makecpt_grd $rmba
                # gmt grdgradient $rmba -A30 -Nt0.6 -Qc -G${rmba}.grad
                # gmt grdimage $rmba -I${rmba}.grad
                # gmt grdcontour $rmba -C10
                # gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                # if [ -f $grid_mask ]; then 
                #     gmt grdimage $grid_mask -C$cpt_bathy -Q # -t$alpha_mask
                # fi 
                gmt grdgradient $bathy_ship -A30 -Nt0.6 -Qc -G${bathy_ship}.grad
                gmt grd2cpt $bathy_ship -Cbasecpt_bathy.cpt 
                gmt grdimage $bathy_ship -I${bathy_ship}.grad
                
                # gmt colorbar -DJCB+o0/1c -Bxaf+l"RMBA ($modelname) $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}
                plotControlTransformFault
                gmt basemap -LjBR+o0.5c/0.5c+w50k+f+u
                # # ============ extract RMBA ==================
                # echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >${dataname}/Results_FFT/averageRMBA.txt
                # # boxnames=(averageBox_OTF averageBox_lFZ averageBox_rFZ averageBox_lOTF averageBox_rOTF)
                # file_boxinfo=${dataname}/input/averageRMA_boxinfo.txt
                # nbox=`awk 'END{print NR}' $file_boxinfo`
                # for (( i=2; i<=$nbox; i++ ));
                # do 
                #     boxname=averageBox_`awk -v row=$i 'NR==row{print $1}' $file_boxinfo`
                #     ind_OTF=`awk -v row=$i 'NR==row{print $2}' $file_boxinfo`
                #     loc_box=`awk -v row=$i 'NR==row{print $3}' $file_boxinfo`
                #     l_box=`awk -v row=$i 'NR==row{print $4}' $file_boxinfo`
                #     w_box=`awk -v row=$i 'NR==row{print $5}' $file_boxinfo`
                #     rot_angle=`awk -v row=$i 'NR==row{print $6}' $file_boxinfo`
                #     boxcolor=`awk -v row=$i 'NR==row{print $7}' $file_boxinfo`
                #     getAverageBox_TF_FZ $boxname $ind_OTF $loc_box $l_box $w_box $rot_angle
                #     gmt plot ${dataname}/input/${boxname}.lonlat -W2p,$boxcolor -L -Gwhite@50
                #     # cut data
                #     gmt grdcut $rmba -F${dataname}/input/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                #     meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                #     echo $meanRMBA
                #     # plot text of mean rmba
                #     echo $lonClatC_box $meanRMBA | gmt text -F+f12p,Helvetica-Bold,$boxcolor=0.1p,black+a$angle_rot
                #     # write to file
                #     meanRMBAs[$i]=$meanRMBA
                #     echo $meanRMBA >>${dataname}/Results_FFT/averageRMBA.txt
                # done
                # gmt basemap -Bafg
                # # calculate delta RMBA and error bar
                # data_dRMBA=`awk 'NR>1{if(NR==2){rmba_FZ1=$1}else{rmba_OTF+=$1; nOTF+=1;}} END{{rmba_OTF-=$1; nOTF-=1; rmba_FZ2=$1; meanRMBA_OTF=rmba_OTF/nOTF; dRMBA=meanRMBA_OTF-(rmba_FZ1+rmba_FZ2)/2.0; err1=meanRMBA_OTF-rmba_FZ1-dRMBA;err1=sqrt(err1*err1)};print 0, dRMBA, err1}' ${dataname}/Results_FFT/averageRMBA.txt`
                # # gmt inset begin -D`awk 'NR==1{print $2}' $file_boxinfo` -Ba -F+gwhite+p1p,black # -F+gwhite+p1p+c0.1c+s
                # #     gmt basemap -JX3c/3c -R-1/1/`echo $data_dRMBA | awk '{print $2-$3*1.5"/"$2+$3*1.5}'` -Ba -BwsEn+gwhite 
                # #     echo $data_dRMBA | gmt plot -Sc0.2c -Gred -W0.5p,black -Ey+p0.5p
                # #     echo $data_dRMBA | awk '{print "@%12%\104@%%RMBA@-TF-FZ@- = "$2" @%12%\261@%% "$3" mGal"}' | gmt text -F+cTL+f5p,Helvetica-Bold -Dj0.1c/0.1
                # # gmt inset end
            gmt subplot set 0,1
                grid_mask=tmp_mask.nc
                alpha_mask=50
                gmt grdmath $bathy_ship $bathy_ship NAN = NAN.nc 
                gmt grdmath $bathy_ship NAN.nc XOR 99999 ADD = $grid_mask
                gmt basemap -JM${figwidth}c -R$rmba -BwSEn -Ba
                makecpt_grd $rmba
                gmt grdgradient $rmba -A30 -Nt0.6 -Qc -G${rmba}.grad
                gmt grdimage $rmba -I${rmba}.grad
                gmt grdcontour $rmba -C10
                gmt grdcontour $rmba -C0, -Wathick,white,- -Wcthin,blue -A0,+f10p,Helvetica,white
                if [ -f $grid_mask ]; then 
                    gmt grdimage $grid_mask -C$cpt_bathy -Q # -t$alpha_mask
                fi 
                # gmt grdgradient $bathy_ship -A30 -Nt0.6 -Qc -G${bathy_ship}.grad
                # gmt grd2cpt $bathy_ship -Cbasecpt_bathy.cpt 
                # gmt grdimage $bathy_ship -I${bathy_ship}.grad
                
                # gmt colorbar -DJCB+o0/1c -Bxaf+l"RMBA ($modelname) $str_shiftData" -By+l"mGal" -G${data_min}/${data_max}
                plotControlTransformFault
                gmt basemap -LjBR+o0.5c/0.5c+w50k+f+u
                # ============ extract RMBA ==================
                echo "#FZ1 OTF1, OTF2, ..., OTFN, FZ2" >${dataname}/Results_FFT/averageRMBA.txt
                # boxnames=(averageBox_OTF averageBox_lFZ averageBox_rFZ averageBox_lOTF averageBox_rOTF)
                file_boxinfo=${dataname}/input/averageRMA_boxinfo.txt
                nbox=`awk 'END{print NR}' $file_boxinfo`
                for (( i=2; i<=$nbox; i++ ));
                do 
                    boxname=averageBox_`awk -v row=$i 'NR==row{print $1}' $file_boxinfo`
                    ind_OTF=`awk -v row=$i 'NR==row{print $2}' $file_boxinfo`
                    loc_box=`awk -v row=$i 'NR==row{print $3}' $file_boxinfo`
                    l_box=`awk -v row=$i 'NR==row{print $4}' $file_boxinfo`
                    w_box=`awk -v row=$i 'NR==row{print $5}' $file_boxinfo`
                    rot_angle=`awk -v row=$i 'NR==row{print $6}' $file_boxinfo`
                    boxcolor=`awk -v row=$i 'NR==row{print $7}' $file_boxinfo`
                    getAverageBox_TF_FZ $boxname $ind_OTF $loc_box $l_box $w_box $rot_angle
                    gmt plot ${dataname}/input/${boxname}.lonlat -W2p,$boxcolor -L -Gwhite@50
                    # cut data
                    gmt grdcut $rmba -F${dataname}/input/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    echo $meanRMBA
                    # plot text of mean rmba
                    echo $lonClatC_box $meanRMBA | gmt text -F+f12p,Helvetica-Bold,$boxcolor=0.1p,black+a$angle_rot
                    # write to file
                    meanRMBAs[$i]=$meanRMBA
                    echo $meanRMBA >>${dataname}/Results_FFT/averageRMBA.txt
                done
                # calculate delta RMBA and error bar
                data_dRMBA=`awk 'NR>1{if(NR==2){rmba_FZ1=$1}else{rmba_OTF+=$1; nOTF+=1;}} END{{rmba_OTF-=$1; nOTF-=1; rmba_FZ2=$1; meanRMBA_OTF=rmba_OTF/nOTF; dRMBA=meanRMBA_OTF-(rmba_FZ1+rmba_FZ2)/2.0; err1=meanRMBA_OTF-rmba_FZ1-dRMBA;err1=sqrt(err1*err1)};print 0, dRMBA, err1}' ${dataname}/Results_FFT/averageRMBA.txt`
                # gmt inset begin -D`awk 'NR==1{print $2}' $file_boxinfo` -Ba -F+gwhite+p1p,black # -F+gwhite+p1p+c0.1c+s
                #     gmt basemap -JX3c/3c -R-1/1/`echo $data_dRMBA | awk '{print $2-$3*1.5"/"$2+$3*1.5}'` -Ba -BwsEn+gwhite 
                #     echo $data_dRMBA | gmt plot -Sc0.2c -Gred -W0.5p,black -Ey+p0.5p
                #     echo $data_dRMBA | awk '{print "@%12%\104@%%RMBA@-TF-FZ@- = "$2" @%12%\261@%% "$3" mGal"}' | gmt text -F+cTL+f5p,Helvetica-Bold -Dj0.1c/0.1
                # gmt inset end
                echo $data_dRMBA | awk '{print "@%12%\104@%%RMBA@-TF-FZ@- = "$2" @%12%\261@%% "$3" mGal"}' | gmt text -F+cTL+f14p,Helvetica-Bold -Dj0.5c/0.5
                # --------- extract RMBA along ridges ---------
                echo "#Ridge1, Ridge2" >${dataname}/Results_FFT/averageRMBA_Ridge.txt
                file_boxinfo_ridge=${dataname}/input/averageRMA_boxinfo_Ridge.txt
                nbox=`awk 'END{print NR}' $file_boxinfo_ridge`
                echo $nbox
                for (( i=2; i<=$nbox; i++ ));
                do 
                    boxname=averageBox_`awk -v row=$i 'NR==row{print $1}' $file_boxinfo_ridge`
                    ind_OTF=`awk -v row=$i 'NR==row{print $2}' $file_boxinfo_ridge`
                    loc_box=`awk -v row=$i 'NR==row{print $3}' $file_boxinfo_ridge`
                    l_box=`awk -v row=$i 'NR==row{print $4}' $file_boxinfo_ridge`
                    w_box=`awk -v row=$i 'NR==row{print $5}' $file_boxinfo_ridge`
                    rot_angle=`awk -v row=$i 'NR==row{print $6}' $file_boxinfo_ridge`
                    boxcolor=`awk -v row=$i 'NR==row{print $7}' $file_boxinfo_ridge`
                    getAverageBox_Ridge $boxname $ind_OTF $loc_box $l_box $w_box $rot_angle
                    gmt plot ${dataname}/input/${boxname}.lonlat -W2p,$boxcolor -L -Gwhite@50
                    # cut data
                    gmt grdcut $rmba -F${dataname}/input/${boxname}.lonlat -Gtmp_cut_box_otf.nc
                    meanRMBA=`gmt grdmath tmp_cut_box_otf.nc MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{printf "%.1f", $3}' && rm tmp.nc tmp_cut_box_otf.nc`
                    echo "Mean RMBA: " $meanRMBA
                    # plot text of mean rmba
                    echo $lonClatC_box $meanRMBA | gmt text -F+f12p,Helvetica-Bold,$boxcolor=0.1p,black+a$angle_rot
                    # write to file
                    meanRMBAs[$i]=$meanRMBA
                    echo $meanRMBA >>${dataname}/Results_FFT/averageRMBA_Ridge.txt
                done
        gmt subplot end
    gmt end show
    rm NAN.nc tmp_mask.nc tmp*
}
function makeBox()
{
    tmp_xmin=$1
    tmp_xmax=$2
    tmp_ymin=$3
    tmp_ymax=$4
    echo $tmp_xmin $tmp_ymin >tmp.box 
    echo $tmp_xmax $tmp_ymin >>tmp.box 
    echo $tmp_xmax $tmp_ymax >>tmp.box 
    echo $tmp_xmin $tmp_ymax >>tmp.box 
}
function getSpreadingType()
{
    rate=$1
    if [ `echo "$rate > 0" | bc` -eq 1 ] && [ `echo "$rate < 20" | bc` -eq 1 ]; then 
        type_sp=Ultraslow
        type_sp_index=0
    elif [ `echo "$rate >= 20" | bc` -eq 1 ] && [ `echo "$rate < 50" | bc` -eq 1 ]; then 
        type_sp=Slow
        type_sp_index=1
    elif [ `echo "$rate >= 50" | bc` -eq 1 ] && [ `echo "$rate < 80" | bc` -eq 1 ]; then 
        type_sp=Intermediate
        type_sp_index=2
    elif [ `echo "$rate >= 80" | bc` -eq 1 ]; then 
        type_sp=Fast
        type_sp_index=3
    fi
}
function figure3()
{
    gmt gmtset FONT_LABEL=7p
    figw=7
    figh=6
    gmt begin Figure3_v1 pdf 
        xmin=0
        xmax=150
        ymin=-6
        ymax=25
        # Fig. a spreading rate - dRMBA
        gmt basemap -JX${figw}c/${figh}c -R0/${xmax}/${ymin}/${ymax} -Baf -Bxa+l"Full spreading rate (mm yr@+-1@+)" -Bya+l"@%12%\104@%%RMBA@-TF-FZ@- (mGal)"
        echo "(a)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1
        # fill different spreading rate region
        
        fs="5p,Helvetica-Bold"
        # ultra-slow 
        h_box=`echo $ymin+1.5 |bc`
        bot_box=`echo $ymin+0.15 |bc`
        y0=`echo "${ymin}+0.75" |bc`
        makeBox 0 20 $ymin $ymax  && gmt plot tmp.box -Ggray@100 -W0p,white@100
        makeBox 0.5 20 $bot_box $h_box  && gmt plot tmp.box -G${color_ridges[0]} -W0.5p,white
        echo 10 $y0 "Ultraslow" | gmt text -F+f${fs},white=0.1p,white
        # slow
        makeBox 20 50 $ymin $ymax && gmt plot tmp.box -Ggray@80 -W0p,white@100
        makeBox 20 50 $bot_box $h_box  && gmt plot tmp.box -G${color_ridges[1]} -W0.5p,white
        echo 35 $y0 "Slow" | gmt text -F+f${fs},white=0.1p,white
        # intermediate
        makeBox 50 80 $ymin $ymax && gmt plot tmp.box -Ggray@100 -W0p,white@100
        makeBox 50 80 $bot_box $h_box  && gmt plot tmp.box -G${color_ridges[2]} -W0.5p,white
        echo 65 $y0 "Intermediate" | gmt text -F+f${fs},white=0.1p,white
        # fast
        makeBox 80 ${xmax} $ymin $ymax && gmt plot tmp.box -Ggray@80 -W0p,white@100
        makeBox 80 `echo ${xmax}-0.5 |bc` $bot_box $h_box  && gmt plot tmp.box -G${color_ridges[3]} -W0.5p,white
        echo 115 $y0 "Fast" | gmt text -F+f${fs},white=0.1p,white
        # plot y=0 line
        echo $xmin 0 >tmp.xy  && echo $xmax 0 >>tmp.xy && gmt plot tmp.xy -W1p,black,-
        # plot dRMBA points
        echo "dRMBA -- spreading rate"
        for i in {0..10};
        do 
            dataname=${TFs[$i]}
            TFName=${TFNames[$i]}
            rate=`cat ${dataname}/input/spreadingRate.txt`
            meanRMBA=${dataname}/Results_FFT/averageRMBA.txt
            data_dRMBA=`awk 'NR>1{if(NR==2){rmba_FZ1=$1}else{rmba_OTF+=$1; nOTF+=1;}} END{{rmba_OTF-=$1; nOTF-=1; rmba_FZ2=$1; meanRMBA_OTF=rmba_OTF/nOTF; dRMBA=meanRMBA_OTF-(rmba_FZ1+rmba_FZ2)/2.0; err1=meanRMBA_OTF-rmba_FZ1-dRMBA;err1=sqrt(err1*err1)};print dRMBA, err1}' $meanRMBA`
            
            # dRMBA=`awk 'NR==2{print $1-($2+$3)/2}' $meanRMBA`
            # dRMBA_l=`awk 'NR==2{print $4-$2}' $meanRMBA`
            # dRMBA_r=`awk 'NR==2{print $5-$3}' $meanRMBA`
            # dRMBA_low=`echo $dRMBA_l $dRMBA_r | awk '{if($1>$2) { print $2} else {print $1}}'`
            # dRMBA_high=`echo $dRMBA_l $dRMBA_r | awk '{if($1>$2) { print $1} else {print $2}}'`
            # errbar=`echo $dRMBA $dRMBA_low $dRMBA_high | awk '{{err1=$2-$1; err2=$3-$1; err_low=err1; err_high=err2; if(err_low>0)err_low=0;}; print err_low, err_high}'`
            # dRMBA=$dRMBA_high
            # echo $errbar
            dRMBA=`echo $data_dRMBA | awk '{print $1}'`
            errbar=`echo $data_dRMBA | awk '{print $2}'`
            mc=black
            getSpreadingType $rate
            mc=${color_ridges[$type_sp_index]}
            # plot
            echo $rate $dRMBA $errbar | gmt plot -Sc0.2c -G$mc -W0.5p,black # -Ey+p0.5p #+a
            echo ${dataname} $rate $dRMBA $errbar
            # echo $rate $dRMBA_l | gmt plot -Sc0.2c -Ggreen -W0.5p,black
            # echo $rate $dRMBA_r | gmt plot -Sc0.2c -Gblue -W0.5p,black
            justify=MC
            dxdy="0 -1"
            if [ "$TFName" == "Oceanographer" ]; then 
                justify=ML 
                dxdy="4 0"
            fi
            if [ "$TFName" == "Marathon" ]; then 
                dxdy="0 1"
            fi 
            # echo 50 20 djls ds | gmt text -F+j${justify}+f5p,Helvetica
            echo $rate $dRMBA $dxdy | awk -v TFname="$TFName" '{print $1+$3, $2+$4, TFname}' | gmt text -F+j${justify}+f5p,Helvetica
        done
        # fig. b AO. - dRMBA
        xmin=0
        xmax=30
        move_x=`echo $figw | awk '{print $1+0.5}'`
        gmt basemap -JX${figw}c/${figh}c -R0/${xmax}/${ymin}/${ymax} -Baf -BwSEn -Bxa+l"Age offset (Myr)" -Bya+l"@%12%\104@%%RMBA@-TF-FZ@- (mGal)" -X${move_x}c
        echo "(b)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1
        echo $xmin 0 >tmp.xy  && echo $xmax 0 >>tmp.xy && gmt plot tmp.xy -W1p,black,-
        echo "dRMBA -- transform age offset"
        rm dRMBA_OA.txt && touch dRMBA_OA.txt
        for i in {0..10};
        do 
            dataname=${TFs[$i]}
            TFName=${TFNames[$i]}
            rate=`cat ${dataname}/input/spreadingRate.txt`
            len=`cat ${dataname}/input/length.txt`
            AO=`echo $rate $len | awk '{print $2/$1*2}'`
            meanRMBA=${dataname}/Results_FFT/averageRMBA.txt
            data_dRMBA=`awk 'NR>1{if(NR==2){rmba_FZ1=$1}else{rmba_OTF+=$1; nOTF+=1;}} END{{rmba_OTF-=$1; nOTF-=1; rmba_FZ2=$1; meanRMBA_OTF=rmba_OTF/nOTF; dRMBA=meanRMBA_OTF-(rmba_FZ1+rmba_FZ2)/2.0; err1=meanRMBA_OTF-rmba_FZ1-dRMBA;err1=sqrt(err1*err1)};print dRMBA, err1}' $meanRMBA`
        
            # dRMBA=`awk 'NR==2{print $1-($2+$3)/2}' $meanRMBA`
            # dRMBA_l=`awk 'NR==2{print $4-$2}' $meanRMBA`
            # dRMBA_r=`awk 'NR==2{print $5-$3}' $meanRMBA`
            # dRMBA_low=`echo $dRMBA_l $dRMBA_r | awk '{if($1>$2) { print $2} else {print $1}}'`
            # dRMBA_high=`echo $dRMBA_l $dRMBA_r | awk '{if($1>$2) { print $1} else {print $2}}'`
            # errbar=`echo $dRMBA $dRMBA_low $dRMBA_high | awk '{{err1=$2-$1; err2=$3-$1; err_low=err1; err_high=err2; if(err_low>0)err_low=0;}; print err_low, err_high}'`
            # dRMBA=$dRMBA_high
            dRMBA=`echo $data_dRMBA | awk '{print $1}'`
            errbar=`echo $data_dRMBA | awk '{print $2}'`
            # plot
            mc=black
            getSpreadingType $rate
            mc=${color_ridges[$type_sp_index]}
            echo $AO $dRMBA $errbar | gmt plot -Sc0.2c -G$mc -W0.5p,black # -Ey+p0.5p #+a
            echo ${dataname} $len $rate $AO $dRMBA $errbar
            echo $AO $dRMBA $errbar>>dRMBA_OA.txt #then use Matlab to make linear fitting
            # # echo $rate $dRMBA_l | gmt plot -Sc0.2c -Ggreen -W0.5p,black
            # # echo $rate $dRMBA_r | gmt plot -Sc0.2c -Gblue -W0.5p,black
            justify=MC
            dxdy="0 -1"
            if [ "$TFName" == "Clipperton" ] || [ "$TFName" == "Marathon" ] || [ "$TFName" == "Chile Ridge" ]; then 
                justify=ML 
                dxdy="0.4 0"
            fi
            echo $AO $dRMBA $dxdy | awk -v TFname="$TFName" '{print $1+$3, $2+$4, TFname}' | gmt text -F+j${justify}+f5p,Helvetica
        done
        # # curve fitting of dRMBA-OA (without Marion and Atlantis II): p=[4.48 0.80]
        # echo "0 0.8048" >tmp.txt
        # echo "0.7019 3.9896" >>tmp.txt
        # echo "5.1605 23.9236" >>tmp.txt
        # gmt plot tmp.txt -W1p,gray,-
    gmt end show 
    rm tmp.box
}
function figure3_v2()
{
    gmt gmtset PS_CHAR_ENCODING=Standard+
    error_bar=$1
    figname=Figure3
    if [ "$error_bar" == "ErrorBar" ]; then 
        figname=Figure3_ErrorBar
    fi
    gmt gmtset FONT_LABEL=7p
    figw=8
    figh=6
    minAO_log=0.25
    maxAO_log=1.5
    gmt makecpt -Cbasecpt_AO.cpt -T0.05/1.5/0.01 --COLOR_BACKGROUND=white@100 --COLOR_FOREGROUND=white@100  > tmp_AO.cpt
    echo "#log10(AO) a AO" > label_cb_AO.txt
    for AO in 2 4 6 8 10 20 30 
    do 
        echo $AO | awk '{print log($1)/log(10), "a", $1}' >> label_cb_AO.txt
    done
    gmt begin $figname pdf 
        xmin=0
        xmax=150
        ymin=-6
        ymax=25
        # Fig. a spreading rate - dRMBA
        gmt basemap -JX${figw}c/${figh}c -R${xmin}/${xmax}/${ymin}/${ymax} -Baf -Bxa+l"Full spreading rate (mm yr@+-1@+)" -Bya+l"@%12%\104@%%RMBA@-TF-FZ@- (mGal)"
        # echo "(a)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1
        # fill different spreading rate region
        
        fs="5p,Helvetica-Bold"
        # ultra-slow 
        h_box=`echo $ymin+1.5 |bc`
        bot_box=`echo $ymin+0.15 |bc`
        y0=`echo "${ymin}+0.75" |bc`
        # see Beaulieu et al.(2015). Fig. 1 https://doi.org/10.1016/j.dsr2.2015.05.001
        # makeBox 0 20 $ymin $ymax  && gmt plot tmp.box -Ggray@100 -W0p,white@100
        makeBox 0.5 20 $bot_box $h_box  && gmt plot tmp.box -W1p,white -Glightgray@50
        echo 10 $y0 "Ultraslow" | gmt text -F+f${fs},${color_ridges[0]} #=0.1p,white
        # slow
        # makeBox 20 50 $ymin $ymax && gmt plot tmp.box -Ggray@80 -W0p,white@100
        makeBox 20 55 $bot_box $h_box  && gmt plot tmp.box -W1p,white -Glightgray@50
        echo 35 $y0 "Slow" | gmt text -F+f${fs},${color_ridges[1]} #=0.1p,white
        # intermediate
        # makeBox 50 80 $ymin $ymax && gmt plot tmp.box -Ggray@100 -W0p,white@100
        makeBox 55 80 $bot_box $h_box  && gmt plot tmp.box -W1p,white  -Glightgray@50
        echo 67.5 $y0 "Intermediate" | gmt text -F+f${fs},${color_ridges[2]} #=0.1p,white
        # fast
        # makeBox 80 ${xmax} $ymin $ymax && gmt plot tmp.box -Ggray@80 -W0p,white@100
        makeBox 80 `echo ${xmax}-0.5 |bc` $bot_box $h_box  && gmt plot tmp.box  -W1p,white -Glightgray@50
        echo 115 $y0 "Fast" | gmt text -F+f${fs},${color_ridges[3]} #=0.1p,white
        # plot y=0 line
        echo $xmin 0 >tmp.xy  && echo $xmax 0 >>tmp.xy && gmt plot tmp.xy -W0.5p,black,.
        # linear fit spreading rate vs dRMBA, dRMBA+err, dRMBA-err
        # area between upper and lower boundary
        gmt plot coefidence_interval.txt -L -Glightgreen@10

        # # gmt plot tmp_dw.txt -W2p,white
        awk '{print $1, $2}' fitline.txt | gmt plot -W0.5p,black,-
        # echo "75 4.8 R@+2@+=0.30" | gmt text -F+f5p,Helvetica,+a-20
        # echo "75 4.8 R@+2@+=0.30" | gmt text -F+f5p,Helvetica,+a-20
        # echo "75 7.5 (-0.17@%12%\261@%%0.01)x + (21@%12%\261@%%0.85)" | gmt text -F+f5p,Helvetica,lightgreen,+a-30
        # plot dRMBA points
        echo "dRMBA -- spreading rate"
        echo "#Spreading rate(mm/yr),age off set (Myr),dRMBA(mGal),error(mGal)" > dRMBA_spreadingRate.txt
        echo "#Spreading rate(mm/yr),age off set (Myr),dRMBA(mGal),error(mGal)" > dRMBA_spreadingRate_fit.txt
        for i in {0..10};
        do 
            dataname=${TFs[$i]}
            TFName=${TFNames[$i]}
            rate=`cat ${dataname}/input/spreadingRate.txt`
            len=`cat ${dataname}/input/length.txt`
            AO=`echo $rate $len | awk '{print $2/$1*2}'`
            meanRMBA=${dataname}/Results_FFT/averageRMBA.txt
            data_dRMBA=`awk 'NR>1{if(NR==2){rmba_FZ1=$1}else{rmba_OTF+=$1; nOTF+=1;}} END{{rmba_OTF-=$1; nOTF-=1; rmba_FZ2=$1; meanRMBA_OTF=rmba_OTF/nOTF; dRMBA=meanRMBA_OTF-(rmba_FZ1+rmba_FZ2)/2.0; err1=meanRMBA_OTF-rmba_FZ1-dRMBA;err1=sqrt(err1*err1)};print dRMBA, err1}' $meanRMBA`
            
            # dRMBA=`awk 'NR==2{print $1-($2+$3)/2}' $meanRMBA`
            # dRMBA_l=`awk 'NR==2{print $4-$2}' $meanRMBA`
            # dRMBA_r=`awk 'NR==2{print $5-$3}' $meanRMBA`
            # dRMBA_low=`echo $dRMBA_l $dRMBA_r | awk '{if($1>$2) { print $2} else {print $1}}'`
            # dRMBA_high=`echo $dRMBA_l $dRMBA_r | awk '{if($1>$2) { print $1} else {print $2}}'`
            # errbar=`echo $dRMBA $dRMBA_low $dRMBA_high | awk '{{err1=$2-$1; err2=$3-$1; err_low=err1; err_high=err2; if(err_low>0)err_low=0;}; print err_low, err_high}'`
            # dRMBA=$dRMBA_high
            # echo $errbar
            dRMBA=`echo $data_dRMBA | awk '{print $1}'`
            errbar=`echo $data_dRMBA | awk '{print $2}'`
            mc=black
            getSpreadingType $rate
            mc=${color_ridges[$type_sp_index]}
            # plot
            name_ridge=${Ridge_names[${Ridge_index[$i]}]}  
            symbol_ridge=${symbols_ridges[${Ridge_index[$i]}]}
            echo "> -S${symbol_ridge}6p ">>dRMBA_spreadingRate.txt
            echo $rate $AO $dRMBA $errbar $TFName $name_ridge | awk '{print $1, $3, log($2)/log(10), $4, $5}' >>dRMBA_spreadingRate.txt
            
            echo $rate "&" $AO "&" $dRMBA $errbar $TFName >>dRMBA_spreadingRate_fit.txt
            # gmt plot -Sc0.2c -G$mc -W0.5p,black # -Ey+p0.5p #+a
            echo ${dataname} $rate $AO $dRMBA $errbar
            # echo $rate $dRMBA_l | gmt plot -Sc0.2c -Ggreen -W0.5p,black
            # echo $rate $dRMBA_r | gmt plot -Sc0.2c -Gblue -W0.5p,black
            justify=MC
            dxdy="0 -1"
            if [ "$TFName" == "Oceanographer" ]; then 
                justify=TL 
                dxdy="1 -0.2"
            fi
            if [ "$TFName" == "Marathon" ]; then 
                dxdy="0 1"
            fi 
            # # echo 50 20 djls ds | gmt text -F+j${justify}+f5p,Helvetica
            # echo $rate $dRMBA $dxdy | awk -v TFname="$TFName" '{print $1+$3, $2+$4, TFname}' | \
            # gmt text -F+j${justify}+f5p,Helvetica # -F+f9p,Helvetica-Bold,black=1p,white+a0
        done
        if [ "$error_bar" == "ErrorBar" ]; then 
            # awk 'NR>1{print $1, $3, log($2)/log(10), $4}' dRMBA_spreadingRate.txt | \
            gmt plot dRMBA_spreadingRate.txt -Sc0.2c -Ctmp_AO.cpt -W0.5p,black -Ey+p0.5p #+a
            gmt plot dRMBA_spreadingRate.txt -Sc0.2c -Ctmp_AO.cpt -W0.5p,black
        else 
            # awk 'NR>1{print $1, $3, log($2)/log(10), $4}' dRMBA_spreadingRate.txt | \
            gmt plot dRMBA_spreadingRate.txt -Sc0.2c -Ctmp_AO.cpt -W0.5p,black #-Ey+p0.5p #+a
        fi
        justify_TFNames=(ML MR MR ML BL ML BL ML MR ML MR)
        dx_TFNames=(3 -3 -1 3 0 3 3 3 -3 3 -3)
        dy_TFNames=(0 0 -1 0 -1.2 0 0 0 0 0 0)
        for i in {0..10};
        do 
            dataname=${TFs[$i]}
            TFName=${TFNames[$i]}
            justify_TF=${justify_TFNames[$i]}
            dxdy="${dx_TFNames[$i]} ${dy_TFNames[$i]}"
            rate=`cat ${dataname}/input/spreadingRate.txt`
            len=`cat ${dataname}/input/length.txt`
            AO=`echo $rate $len | awk '{print $2/$1*2}'`
            meanRMBA=${dataname}/Results_FFT/averageRMBA.txt
            data_dRMBA=`awk 'NR>1{if(NR==2){rmba_FZ1=$1}else{rmba_OTF+=$1; nOTF+=1;}} END{{rmba_OTF-=$1; nOTF-=1; rmba_FZ2=$1; meanRMBA_OTF=rmba_OTF/nOTF; dRMBA=meanRMBA_OTF-(rmba_FZ1+rmba_FZ2)/2.0; err1=meanRMBA_OTF-rmba_FZ1-dRMBA;err1=sqrt(err1*err1)};print dRMBA, err1}' $meanRMBA`
            
            # dRMBA=`awk 'NR==2{print $1-($2+$3)/2}' $meanRMBA`
            # dRMBA_l=`awk 'NR==2{print $4-$2}' $meanRMBA`
            # dRMBA_r=`awk 'NR==2{print $5-$3}' $meanRMBA`
            # dRMBA_low=`echo $dRMBA_l $dRMBA_r | awk '{if($1>$2) { print $2} else {print $1}}'`
            # dRMBA_high=`echo $dRMBA_l $dRMBA_r | awk '{if($1>$2) { print $1} else {print $2}}'`
            # errbar=`echo $dRMBA $dRMBA_low $dRMBA_high | awk '{{err1=$2-$1; err2=$3-$1; err_low=err1; err_high=err2; if(err_low>0)err_low=0;}; print err_low, err_high}'`
            # dRMBA=$dRMBA_high
            # echo $errbar
            dRMBA=`echo $data_dRMBA | awk '{print $1}'`
            errbar=`echo $data_dRMBA | awk '{print $2}'`
            mc=black
            getSpreadingType $rate
            mc=${color_ridges[$type_sp_index]}
            # plot
            name_ridge=${Ridge_names[${Ridge_index[$i]}]}  
            symbol_ridge=${symbols_ridges[${Ridge_index[$i]}]}
            echo "> -S${symbol_ridge}6p ">>dRMBA_spreadingRate.txt
            echo $rate $AO $dRMBA $errbar $TFName $name_ridge | awk '{print $1, $3, log($2)/log(10), $4, $5}' >>dRMBA_spreadingRate.txt
            
            # echo $rate $AO $dRMBA $errbar >>dRMBA_spreadingRate_fit.txt
            # gmt plot -Sc0.2c -G$mc -W0.5p,black # -Ey+p0.5p #+a
            echo ${dataname} $rate $AO $dRMBA $errbar
            # echo $rate $dRMBA_l | gmt plot -Sc0.2c -Ggreen -W0.5p,black
            # echo $rate $dRMBA_r | gmt plot -Sc0.2c -Gblue -W0.5p,black
            # echo 50 20 djls ds | gmt text -F+j${justify}+f5p,Helvetica
            echo $rate $dRMBA "$dxdy" | awk -v TFname="$TFName" -v justify=${justify_TF} '{print $1+$3, $2+$4, "5p,Helvetica 0 ", justify, TFname}' | \
            gmt text  -Dj0p -F+f+a+j #-F+f5p,Helvetica-Bold,black=0.2p,white #-F+j${justify}+f5p,Helvetica #
        done
        gmt colorbar -Ctmp_AO.cpt -DJTR+o-3.1c/-1.2c+w2.9c/0.15c+h -Bpxclabel_cb_AO.txt+l"Age offset (Myr)" -G${minAO_log}/${maxAO_log} --FONT_LABEL=6p --FONT_ANNOT_PRIMARY=6p --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1.5p --MAP_LABEL_OFFSET=2p --MAP_ANNOT_OFFSET_PRIMARY=2p
        # legend
        edgewidth_symbol=0.5p
        offset_symbol_text=0.3c
        linespacing=1.5p
gmt legend -DjRT+w3c+o0.1c+l${linespacing} -F+p0.4p+gbeige+swhite@100 --FONT_ANNOT_PRIMARY=5p <<- EOF
N 3
S 0.1c ${symbols_ridges[0]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[0]}
S 0.1c ${symbols_ridges[1]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[1]}
S 0.1c ${symbols_ridges[2]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[2]}
S 0.1c ${symbols_ridges[3]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[3]}
S 0.1c ${symbols_ridges[4]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[4]}
S 0.1c ${symbols_ridges[5]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[5]}
EOF
    gmt end show 
    rm tmp.box
}
# dRMBA: TF - Ridge
function figure3_dRMBA_T2R()
{
    gmt gmtset PS_CHAR_ENCODING=Standard+
    error_bar=$1
    figname=Figure3_TF_Ridge
    if [ "$error_bar" == "ErrorBar" ]; then 
        figname=Figure3_TF_Ridge_ErrorBar
    fi
    gmt gmtset FONT_LABEL=7p
    figw=8
    figh=6
    minAO_log=0.25
    maxAO_log=1.5
    gmt makecpt -Cbasecpt_AO.cpt -T0.05/1.5/0.01 --COLOR_BACKGROUND=white@100 --COLOR_FOREGROUND=white@100  > tmp_AO.cpt
    echo "#log10(AO) a AO" > label_cb_AO.txt
    for AO in 2 4 6 8 10 20 30 
    do 
        echo $AO | awk '{print log($1)/log(10), "a", $1}' >> label_cb_AO.txt
    done
    gmt begin $figname pdf 
        xmin=0
        xmax=150
        ymin=-30
        ymax=30
        # Fig. a spreading rate - dRMBA
        gmt basemap -JX${figw}c/${figh}c -R${xmin}/${xmax}/${ymin}/${ymax} -Baf -Bxa+l"Full spreading rate (mm yr@+-1@+)" -Bya+l"@%12%\104@%%RMBA@-T-R@- (mGal)"
        # echo "(a)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1
        # fill different spreading rate region
        
        fs="5p,Helvetica-Bold"
        # ultra-slow 
        h_box=`echo $ymin+3.5 |bc`
        bot_box=`echo $ymin+0.15 |bc`
        y0=`echo "${ymin}+1.75" |bc`
        # see Beaulieu et al.(2015). Fig. 1 https://doi.org/10.1016/j.dsr2.2015.05.001
        # makeBox 0 20 $ymin $ymax  && gmt plot tmp.box -Ggray@100 -W0p,white@100
        makeBox 0.5 20 $bot_box $h_box  && gmt plot tmp.box -W1p,white -Glightgray@50
        echo 10 $y0 "Ultraslow" | gmt text -F+f${fs},${color_ridges[0]} #=0.1p,white
        # slow
        # makeBox 20 50 $ymin $ymax && gmt plot tmp.box -Ggray@80 -W0p,white@100
        makeBox 20 55 $bot_box $h_box  && gmt plot tmp.box -W1p,white -Glightgray@50
        echo 35 $y0 "Slow" | gmt text -F+f${fs},${color_ridges[1]} #=0.1p,white
        # intermediate
        # makeBox 50 80 $ymin $ymax && gmt plot tmp.box -Ggray@100 -W0p,white@100
        makeBox 55 80 $bot_box $h_box  && gmt plot tmp.box -W1p,white  -Glightgray@50
        echo 67.5 $y0 "Intermediate" | gmt text -F+f${fs},${color_ridges[2]} #=0.1p,white
        # fast
        # makeBox 80 ${xmax} $ymin $ymax && gmt plot tmp.box -Ggray@80 -W0p,white@100
        makeBox 80 `echo ${xmax}-0.5 |bc` $bot_box $h_box  && gmt plot tmp.box  -W1p,white -Glightgray@50
        echo 115 $y0 "Fast" | gmt text -F+f${fs},${color_ridges[3]} #=0.1p,white
        # plot y=0 line
        echo $xmin 0 >tmp.xy  && echo $xmax 0 >>tmp.xy && gmt plot tmp.xy -W0.5p,black,.
        # linear fit spreading rate vs dRMBA, dRMBA+err, dRMBA-err
        # area between upper and lower boundary
        # gmt plot coefidence_interval.txt -L -Glightgreen@10

#         # # gmt plot tmp_dw.txt -W2p,white
#         awk '{print $1, $2}' fitline.txt | gmt plot -W0.5p,black,-
#         # echo "75 4.8 R@+2@+=0.30" | gmt text -F+f5p,Helvetica,+a-20
#         # echo "75 4.8 R@+2@+=0.30" | gmt text -F+f5p,Helvetica,+a-20
#         # echo "75 7.5 (-0.17@%12%\261@%%0.01)x + (21@%12%\261@%%0.85)" | gmt text -F+f5p,Helvetica,lightgreen,+a-30
#         # plot dRMBA points
#         echo "dRMBA -- spreading rate"
        echo "#Spreading rate(mm/yr),age off set (Myr),dRMBA(mGal),error(mGal)" > dRMBA_TR_spreadingRate.txt
        for i in {0..10};
        do 
            dataname=${TFs[$i]}
            TFName=${TFNames[$i]}
            rate=`cat ${dataname}/input/spreadingRate.txt`
            len=`cat ${dataname}/input/length.txt`
            AO=`echo $rate $len | awk '{print $2/$1*2}'`
            meanRMBA_OTF=`awk 'NR>1{if(NR==2){rmba_FZ1=$1}else{rmba_OTF+=$1; nOTF+=1;}} END{{rmba_OTF-=$1; nOTF-=1; rmba_FZ2=$1; meanRMBA_OTF=rmba_OTF/nOTF; dRMBA=meanRMBA_OTF-(rmba_FZ1+rmba_FZ2)/2.0; err1=meanRMBA_OTF-rmba_FZ1-dRMBA;err1=sqrt(err1*err1)};print meanRMBA_OTF}' ${dataname}/Results_FFT/averageRMBA.txt`
            RMBA_Ridge1=`awk 'NR==2{print $1}' ${dataname}/Results_FFT/averageRMBA_Ridge.txt`
            RMBA_Ridge2=`awk 'NR==3{print $1}' ${dataname}/Results_FFT/averageRMBA_Ridge.txt`
            # echo $dataname $meanRMBA_OTF $RMBA_Ridge1 $RMBA_Ridge2
            dRMBA=`echo $meanRMBA_OTF $RMBA_Ridge1 $RMBA_Ridge2 | awk '{print $1-($2+$3)/2.0}'`
            errbar=`echo $meanRMBA_OTF $RMBA_Ridge1 $dRMBA | awk '{print sqrt(($1-$2-$3)*($1-$2-$3))}'`
            mc=black
            getSpreadingType $rate
            mc=${color_ridges[$type_sp_index]}
            # plot
            name_ridge=${Ridge_names[${Ridge_index[$i]}]}  
            symbol_ridge=${symbols_ridges[${Ridge_index[$i]}]}
            echo "> -S${symbol_ridge}6p ">>dRMBA_TR_spreadingRate.txt
            echo $rate $AO $dRMBA $errbar $TFName $name_ridge | awk '{print $1, $3, log($2)/log(10), $4, $5}' >>dRMBA_TR_spreadingRate.txt
            echo ${dataname} $rate $AO $dRMBA $errbar
        done
        if [ "$error_bar" == "ErrorBar" ]; then 
            # awk 'NR>1{print $1, $3, log($2)/log(10), $4}' dRMBA_TR_spreadingRate.txt | \
            gmt plot dRMBA_TR_spreadingRate.txt -Sc0.2c -Ctmp_AO.cpt -W0.5p,black -Ey+p0.5p #+a
            gmt plot dRMBA_TR_spreadingRate.txt -Sc0.2c -Ctmp_AO.cpt -W0.5p,black
        else 
            # awk 'NR>1{print $1, $3, log($2)/log(10), $4}' dRMBA_TR_spreadingRate.txt | \
            gmt plot dRMBA_TR_spreadingRate.txt -Sc0.2c -Ctmp_AO.cpt -W0.5p,black #-Ey+p0.5p #+a
        fi
        justify_TFNames=(ML MR ML ML BL ML BL BC MR ML MR)
        dx_TFNames=(3 -3 3 3 -4 3 3 0 -3 3 -3)
        dy_TFNames=(0 0 0 0 -2 0 0 -4.5 0 0 0)
        for i in {0..10};
        do 
            dataname=${TFs[$i]}
            TFName=${TFNames[$i]}
            justify_TF=${justify_TFNames[$i]}
            dxdy="${dx_TFNames[$i]} ${dy_TFNames[$i]}"
            rate=`cat ${dataname}/input/spreadingRate.txt`
            len=`cat ${dataname}/input/length.txt`
            AO=`echo $rate $len | awk '{print $2/$1*2}'`
            meanRMBA_OTF=`awk 'NR>1{if(NR==2){rmba_FZ1=$1}else{rmba_OTF+=$1; nOTF+=1;}} END{{rmba_OTF-=$1; nOTF-=1; rmba_FZ2=$1; meanRMBA_OTF=rmba_OTF/nOTF; dRMBA=meanRMBA_OTF-(rmba_FZ1+rmba_FZ2)/2.0; err1=meanRMBA_OTF-rmba_FZ1-dRMBA;err1=sqrt(err1*err1)};print meanRMBA_OTF}' ${dataname}/Results_FFT/averageRMBA.txt`
            RMBA_Ridge1=`awk 'NR==2{print $1}' ${dataname}/Results_FFT/averageRMBA_Ridge.txt`
            RMBA_Ridge2=`awk 'NR==3{print $1}' ${dataname}/Results_FFT/averageRMBA_Ridge.txt`
            # echo $dataname $meanRMBA_OTF $RMBA_Ridge1 $RMBA_Ridge2
            dRMBA=`echo $meanRMBA_OTF $RMBA_Ridge1 $RMBA_Ridge2 | awk '{print $1-($2+$3)/2.0}'`
            errbar=`echo $meanRMBA_OTF $RMBA_Ridge1 $dRMBA | awk '{print sqrt(($1-$2-$3)*($1-$2-$3))}'`
            mc=black
            getSpreadingType $rate
            mc=${color_ridges[$type_sp_index]}
            # plot
            name_ridge=${Ridge_names[${Ridge_index[$i]}]}  
            symbol_ridge=${symbols_ridges[${Ridge_index[$i]}]}
            echo "> -S${symbol_ridge}6p ">>dRMBA_TR_spreadingRate.txt
            echo $rate $AO $dRMBA $errbar $TFName $name_ridge | awk '{print $1, $3, log($2)/log(10), $4, $5}' >>dRMBA_TR_spreadingRate.txt
            
            # echo $rate $AO $dRMBA $errbar >>dRMBA_TR_spreadingRate_fit.txt
            # gmt plot -Sc0.2c -G$mc -W0.5p,black # -Ey+p0.5p #+a
            echo ${dataname} $rate $AO $dRMBA $errbar
            # echo $rate $dRMBA_l | gmt plot -Sc0.2c -Ggreen -W0.5p,black
            # echo $rate $dRMBA_r | gmt plot -Sc0.2c -Gblue -W0.5p,black
            # echo 50 20 djls ds | gmt text -F+j${justify}+f5p,Helvetica
            echo $rate $dRMBA "$dxdy" | awk -v TFname="$TFName" -v justify=${justify_TF} '{print $1+$3, $2+$4, "5p,Helvetica 0 ", justify, TFname}' | \
            gmt text  -Dj0p -F+f+a+j #-F+f5p,Helvetica-Bold,black=0.2p,white #-F+j${justify}+f5p,Helvetica #
        done
        gmt colorbar -Ctmp_AO.cpt -DJTR+o-3.1c/-1.2c+w2.9c/0.15c+h -Bpxclabel_cb_AO.txt+l"Age offset (Myr)" -G${minAO_log}/${maxAO_log} --FONT_LABEL=6p --FONT_ANNOT_PRIMARY=6p --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1.5p --MAP_LABEL_OFFSET=2p --MAP_ANNOT_OFFSET_PRIMARY=2p
        # legend
        edgewidth_symbol=0.5p
        offset_symbol_text=0.3c
        linespacing=1.5p
gmt legend -DjRT+w3c+o0.1c+l${linespacing} -F+p0.4p+gbeige+swhite@100 --FONT_ANNOT_PRIMARY=5p <<- EOF
N 3
S 0.1c ${symbols_ridges[0]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[0]}
S 0.1c ${symbols_ridges[1]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[1]}
S 0.1c ${symbols_ridges[2]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[2]}
S 0.1c ${symbols_ridges[3]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[3]}
S 0.1c ${symbols_ridges[4]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[4]}
S 0.1c ${symbols_ridges[5]} 5p - ${edgewidth_symbol} ${offset_symbol_text} ${Ridge_names[5]}
EOF
    gmt end show 
    rm tmp.box
}
# get vertical slice from thermal model result (vtu file)
function getSlice_thermalModel()
{
    lon1=$1
    lat1=$2
    lon2=$3
    lat2=$4
    lonlat1="$lon1 $lat1"
    lonlat2="$lon2 $lat2"
    thermal_VTU=`ls ../geodynamic-models_3D/${dataname}/IsoVisco/*.vtu | awk '{print $1}'`
    echo "Thermal model result: " $thermal_VTU
    lon0=`awk 'NR==1{print $1}' ${OTF}` 
    lat0=`awk 'NR==1{print $2}' ${OTF}` 
    x0=0
    y0=0
    cat $OTF | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 >${dataname}_OTF_UTM.txt
    # awk 'NR==2{print "length of control line: "sqrt($1*$1+$2*$2)}' input/${dataname}_OTF_UTM.txt
    # cat $ROTF | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 >tmp_ROTF_UTM.txt
    # calculate cos(theta) and sin(theta)
    costheta=`awk 'NR==2{print $1/sqrt($1*$1+$2*$2)}' ${dataname}_OTF_UTM.txt`
    sintheta=`awk 'NR==2{print $2/sqrt($1*$1+$2*$2)}' ${dataname}_OTF_UTM.txt`
    angle_rot=`echo $sintheta $costheta | awk '{printf "%.1f", atan2($1, $2)/3.141592653*180}'`
    echo "Rotation angle: "${angle_rot}" degree"
    # endpoints lonlat to xy
    xy1=`echo $lonlat1 | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 | awk '{print $1, $2}'`
    xy2=`echo $lonlat2 | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 | awk '{print $1, $2}'`
    # rot xy to thermal model CS
    xy1=`echo $xy1 | awk -v co=$costheta -v si=$sintheta '{print co*$1+si*$2, -si*$1+co*$2}'`
    xy2=`echo $xy2 | awk -v co=$costheta -v si=$sintheta '{print co*$1+si*$2, -si*$1+co*$2}'`
    # call pvpython script
    slicefile_suffix=${lon1}_${lat1}_${lon2}_${lat2}
    ./getSlice_VTU.py $thermal_VTU $xy1 $xy2 $slicefile_suffix
    # # get saved slice file name
    slice=`ls ../geodynamic-models_3D/${dataname}/IsoVisco/*${slicefile_suffix}.csv | awk '{print $1}'`
    # # 1. rotate back x,y
    # awk -v co=$costheta -v si=$sintheta -F ',' 'NR>1{{x=$2; y=$3; z=$4; T=$1;}; print co*x-si*y, si*x+co*y, x, y, z, T}' $slice >slice_rot.xyz
    # # 2. xy to lonlat
    # echo "#lon,lat,(x,y,z)[m], T[deg.C]" >$slice
    # cat slice_rot.xyz | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 -I -f "%.8E" | awk '{print $1, $2, $3, $4, $5, $6-273}' >>$slice
    
    # plot 
    figwidth=14
    faa=${dataname}/input/${dataname}_faa.nc
    gmt begin ../geodynamic-models_3D/${dataname}/IsoVisco/slice_${slicefile_suffix} pdf 
        gmt basemap -JM${figwidth}c -R$faa -Ba 
        plotControlTransformFault
        # input endpoints
        echo $lonlat1 >tmp.xy
        echo $lonlat2 >>tmp.xy 
        gmt plot tmp.xy -W1p,red
        # extracted profile
        awk '{print $1, $2}' $slice | gmt plot -Sc0.01c -W0.1p,blue
        # plot plrofile
        T=$(gmt info -T100+c5 $slice)
        gmt makecpt -Crainbow $T
        distance=`echo $xy1 $xy2 | awk '{print sqrt(($3-$1)*($3-$1) + ($4-$2)*($4-$2))/1000}'`
        gmt basemap -JX6c/-2c -R0/$distance/0/20 -Ba -Bxa+l"Distance from left endpoint (km)" -Bya+l"Depth bsf (km)" -X1.5c -Y1c
        x0y0=($xy1)
        awk -v x0=${x0y0[0]} -v y0=${x0y0[1]} 'NR>1{print sqrt(($3-x0)*($3-x0) + ($4-y0)*($4-y0))/1000, 100-$5/1000, $6}' $slice >tmp_T.xyz
        gmt contour tmp_T.xyz -C -I 
        clevels="200,400,600,800"
        gmt contour tmp_T.xyz -C$clevels -Wthin -A${clevels}+f7p,Helvetica,white
        gmt colorbar -DJCT+o0/0.1c -Bxaf+l"T (deg.C)" --FONT_LABEL=9p
    gmt end show
    rm tmp* slice_rot.xyz ${dataname}_OTF_UTM.txt
}
function extractProperties_thermalModel_lonlat()
{
    ext=dat
    loc_lonlat=$1
    ASPECT_model_name=$2
    if [ ! -f ${loc_lonlat}.${ext} ]; then 
        echo "File doesn't exist: ${loc_lonlat}.${ext}"
        exit
    fi 
    thermal_VTU=`ls ../geodynamic-models_3D/${dataname}/${ASPECT_model_name}/*.vtu | awk '{print $1}'`
    echo "Thermal model result: " $thermal_VTU
    lon0=`awk 'NR==1{print $1}' ${OTF}` 
    lat0=`awk 'NR==1{print $2}' ${OTF}` 
    x0=0
    y0=0
    cat $OTF | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 >${dataname}_OTF_UTM.txt
    # calculate cos(theta) and sin(theta)
    costheta=`awk 'NR==2{print $1/sqrt($1*$1+$2*$2)}' ${dataname}_OTF_UTM.txt`
    sintheta=`awk 'NR==2{print $2/sqrt($1*$1+$2*$2)}' ${dataname}_OTF_UTM.txt`
    angle_rot=`echo $sintheta $costheta | awk '{printf "%.1f", atan2($1, $2)/3.141592653*180}'`
    echo "Rotation angle: "${angle_rot}" degree"
    # lonlat to xy
    awk 'NR>1{print $2, $1, $3}' ${loc_lonlat}.${ext} | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 | awk '{print $1, $2, $3}' > tmp_UTM.xyz
    # rot xy to thermal model CS
    awk '{print $1, $2, $3}' tmp_UTM.xyz | awk -v co=$costheta -v si=$sintheta -v z_sf=$z_seafloor_ASEPCTmodel '{print co*$1+si*$2, -si*$1+co*$2, z_sf-$3*1000}' >tmp.xyz
    # call pvpython script
    ./py_VTKtrack.py $thermal_VTU tmp.xyz
    
    # 1. rotate back x,y
    awk -v co=$costheta -v si=$sintheta '{{x=$1; y=$2; z=$3; T=$4;}; print co*x-si*y, si*x+co*y, x, y, z, T}' tmp_T.xyz >tmp_T_rot.xyz
    # 2. xy to lonlat
    echo "#lon,lat,(x,y,z)[m], T[C]" >${loc_lonlat}_${ASPECT_model_name}_T.${ext}
    cat tmp_T_rot.xyz | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon0 +lat_0=$lat0 +x_0=$x0 +y_0=$y0 -I -f "%.8E" | awk '{print $1, $2, $3, $4, $5, $6-273}' >>${loc_lonlat}_${ASPECT_model_name}_T.${ext}
    # plot
    range_T=$(gmt info -T50+c5 ${loc_lonlat}_${ASPECT_model_name}_T.${ext})
    gmt makecpt -Crainbow $range_T >tmp_T.cpt
    figwidth=14
    gmt begin ${dataname}/figures/EQ_T_${ASPECT_model_name} pdf 
        gmt basemap -JM${figwidth}c -R$faa -Ba 
        plotControlTransformFault
        awk '{print $1, $2, $6}' ${loc_lonlat}_${ASPECT_model_name}_T.${ext} | gmt plot -Sc0.1c -Ctmp_T.cpt
        gmt colorbar -DJCB+o0/-2c -Bxaf+l"T (deg.C)" -Ctmp_T.cpt
        # plot vertical profile
        gmt basemap -JX8c/-3c -R0/360/0/20 -Ba -Bxa+l"Distance from left RTI (km)" -Bya+l"Depth bsf (km)" -X6c -Y7.5c
        echo "0 7.5" >tmp.7.5
        echo "360 7.5" >>tmp.7.5
        gmt plot tmp.7.5 -W0.5p,purple,- 
        awk 'NR>1{print sqrt($3*$3 + $4*$4)/1000, 100-$5/1000, $6}' ${loc_lonlat}_${ASPECT_model_name}_T.${ext} | \
        gmt plot -Sc0.1c -Ctmp_T.cpt -W0.3p,black
    gmt end show
    # rm tmp*  ${dataname}_OTF_UTM.txt
}
function figure1()
{
    TFNames=(Marathon     Marion      AtlantisII "Marie Celeste"           Oceanographer Atlantis     Vlamingh  "CR-TF 39@%12%\260@%%S"          Garrett Gofar Clipperton)
    loc_TFNames=(BL       BR          TL          TL                        TL            TL           BR       TR                              BL      TL     BR)
    TFs=(    MAR_Marathon SWIR_Marion AtlantisII MarieCeleste Oceanographer MAR_Atlantis Vlamingh Chile_Ridge Garrett Gofar Clipperton)
    offset=(4.1c/5.3c     5.1c/0.1c      11.9c/2.4c 9.4c/6.4c    4.1c/7.3c     7.2c/4.8c    9.6c/0.1c 0.3c/1.6c  4.6c/3.3c 0.3c/4.4c 0.3c/5.8c)
    arrows=(6.4/5.05/6.1/5.3 9.41/2.17/7.1/2.17 10.4/2.8/11.9/2.8 10.8/3.6/11.4/6.4 6.8/6.15/6/7.3 6.5/5.9/7.2/5.9 11.3/2.35/11.3/1.8 4.6/2.5/2.3/2.5 3.85/3.75/4.6/3.75 4/4.2/2.35/5 4.1/4.9/2.35/6)
    figwidth=14
    fontsize=7 #pt
    proj=-JQ-30/37.5/${figwidth}c #-JR0/${figwidth}c #-JK0/${figwidth}c #-JH0/${figwidth}c
    figheight=`gmt_get_map_height -Rg  $proj`
    gmt begin Figure1 png E900
        gmt gmtset MAP_FRAME_WIDTH=1p
        etopo=${globalDataPath}/etopo/ETOPO10.nc
        # gmt grd2cpt $etopo -Cgray #-Ic
        gmt basemap -Rd $proj -A0 
        # gmt grdgradient $etopo -A30 -Nt2.5 -Qc -Gtmp_etopo.grad
        # gmt grdimage ${etopo} -Itmp_etopo.grad -Cgray
        gmt coast -Ggray -Swhite -Da -W0.1p,gray@100
    gmt end

    gmt begin Figure1 pdf 
        gmt gmtset MAP_FRAME_WIDTH=1p
        etopo=${globalDataPath}/etopo/ETOPO10.nc
        # gmt grd2cpt $etopo -Cgray #-Ic
        gmt basemap -Rd $proj -A0 
        gmt grdgradient $etopo -A30 -Nt2.5 -Qc -Gtmp_etopo.grad
        gmt grdimage ${etopo} -Itmp_etopo.grad -Cgray
        gmt coast -Ggray@50 -Swhite -A1000 -Dl -Bf60 -W0.1p,gray@20
        # use a png mask to get rid of "ghost" grid lines which are caused by patches edges of coast datasets
        gmt image Figure1.png -Dx0c/0c+w${figwidth}c -Ggray+t
        # spreading ridges
        # gmt plot ${globalDataPath}/MOR/Muller_etal_AREPS_2016_Ridges.xy -W0.5p,black
        # https://tos.org/oceanography/assets/docs/20-1_snow.pdf fig.1
        for ridge in `ls ${globalDataPath}/MOR/*.txt`
        do 
            lc=black
            if [ `awk 'NR==1{print $2}' $ridge ` == Ultraslow ]; then 
                lc=${color_ridges[0]}
            elif [ `awk 'NR==1{print $2}' $ridge ` == Slow ]; then 
                lc=${color_ridges[1]}
            elif [ `awk 'NR==1{print $2}' $ridge ` == Intermediate ]; then 
                lc=${color_ridges[2]}
            elif [ `awk 'NR==1{print $2}' $ridge ` == Fast ]; then 
                lc=${color_ridges[3]}
            else
                echo $ridge
            fi 
            gmt plot $ridge -W1p,$lc
        done 
        # location of TFs
        for i in {0..10};
        do 
            dataname=${TFs[$i]}
            rmba=${dataname}/Results_FFT/${dataname}_rmba.nc 
            # plot
            range_rmba=`gmt_get_gridregion $rmba`
            xylim=`echo $range_rmba | awk -F '[/]' '{print $1, $2, $3, $4}'`
            makeBox $xylim
            gmt plot tmp.box -Gblack@50
        done
        for i in {0..10};
        do 
            dataname=${TFs[$i]}
            TFName=${TFNames[$i]}
            loc_TFName=${loc_TFNames[$i]}
            rmba=${dataname}/Results_FFT/${dataname}_rmba.nc 
            ROTF=${dataname}/input/${dataname}_ridge_otf.txt
            OTF=${dataname}/input/${dataname}_OTF.txt
            rate=`cat ${dataname}/input/spreadingRate.txt`
            # mask
            bathy_ship=${dataname}/input/${dataname}_ship.nc
            grid_mask=${dataname}/Results_FFT/${dataname}_mask.nc 
            if [ ! -f $grid_mask ]; then 
                gmt grdmath $bathy_ship $bathy_ship NAN = NAN.nc 
                gmt grdmath $bathy_ship NAN.nc XOR 99999 ADD = $grid_mask
            fi 
            # inset 
            min_rmba=-40
            max_rmba=30
            gmt makecpt -C${masterCPT_grav}+h -T${min_rmba}/${max_rmba}
            w_inset=2
            if [ "${dataname}" == "Oceanographer" ]; then
                w_inset=2.5
            fi
            if [ "${dataname}" == "Garrett" ]; then
                w_inset=2.5
            fi
            figw=$w_inset
            figh=`gmt_get_map_height -R$rmba  -JM${figw}c`
            # getSpreadingType $rate
            color_mor=black #${color_ridges[$type_sp_index]}
            gmt inset begin -DjBL+w${figw}c/${figh}c+o${offset[$i]} -F+gwhite+p1p,$color_mor # -F+gwhite+p1p+c0.1c+s
                gmt basemap -JM${figw}c -R$rmba -Ba -Bwsen --MAP_FRAME_PEN=thick,red@100
                gmt grdimage $rmba #-I${rmba}.grad 
                gmt grdimage $grid_mask -C$cpt_bathy -Q #-t50
                RTI_ms=0.1
                plotControlTransformFault
                if [ "${TFName}" == "Marie Celeste" ]; then
                    echo "Marie" | gmt text -F+c${loc_TFName}+f${fontsize}p -Dj0.1c/0.1c -Gwhite@20 
                    echo "Celeste" | gmt text -F+c${loc_TFName}+f${fontsize}p -Dj0.1c/0.4c -Gwhite@20 
                else
                    echo $TFName | gmt text -F+c${loc_TFName}+f${fontsize}p -Dj0.1c/0.1c -Gwhite@20 
                fi
            gmt inset end
        done
        # # plot help lines
        # gmt basemap -JX${figwidth}/${figheight} -R0/$figwidth/0/$figheight  -Bsg0.2 --MAP_GRID_PEN_SECONDARY=0.5p,green@50
        # gmt basemap -Ba1f1g1 --MAP_GRID_PEN_PRIMARY=1p,blue@50
        gmt basemap -JX${figwidth}/${figheight} -R0/$figwidth/0/$figheight -B0 -Bwsen --MAP_FRAME_PEN=thick,red@100
        for i in {0..10};
        do 
            dataname=${TFs[$i]}
            rate=`cat ${dataname}/input/spreadingRate.txt`
            getSpreadingType $rate
            color_arrow=${color_ridges[$type_sp_index]}
            echo ${arrows[$i]} | awk -F '[/]' '{print $1, $2, $3, $4}' | \
            gmt plot -Sv0.2c+s+e+g$color_arrow -W0.5p,$color_arrow -Gblack
        done
        gmt colorbar -DjBL+o0.4c/0.2c+w2.5c/0.1c+h+m -Bx+l"RMBA" -By+l"mGal" -F+gbeige+p0.2p --FONT_LABEL=${fontsize}p --FONT_ANNOT_PRIMARY=${fontsize}p --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
gmt legend -DjLT+w3.6c+o0.02c -F+p0.4p+gbeige+swhite@100 --FONT_ANNOT_PRIMARY=${fontsize}p <<- EOF
H ${fontsize}p Spreading rate 
D 0.2p
N 2
S 0.1c - 0.3c - 2p,magenta 0.4c Ultraslow
S 0.5c - 0.3c - 2p,red 0.8c Slow
S 0.1c - 0.3c - 2p,blue 0.4c Intermediate
S 0.5c - 0.3c - 2p,green4 0.8c Fast
EOF
    gmt end show
    
    rm tmp* NAN.nc Figure1.png
}
function plotControlTransformFault_fig2()
{
    text_OTF=($1)
    text_FZ=$2
    text_NTO=$3
    plotBox=$4

    OTF=${dataname}/input/${dataname}_OTF.txt
    ROTF=${dataname}/input/${dataname}_ridge_otf.txt
    ridgeIndex=${dataname}/input/ridge_index.txt
    c_OTF=red 
    c_FZ=blue
    

    awk 'NR>=1 && NR<=2 {print}' $ROTF | gmt plot -W2p,white
    awk 'NR>=1 && NR<=2 {print}' $ROTF | gmt plot -W0.5p,black,-
    awk 'NR>=3 && NR<=4 {print}' $ROTF | gmt plot -W2p,white
    awk 'NR>=3 && NR<=4 {print}' $ROTF | gmt plot -W0.5p,black,-
    awk 'NR>=5 && NR<=6 {print}' $ROTF | gmt plot -W2p,white
    awk 'NR>=5 && NR<=6 {print}' $ROTF | gmt plot -W0.5p,black,-
    awk 'NR>=7 && NR<=8 {print}' $ROTF | gmt plot -W2p,white
    awk 'NR>=7 && NR<=8 {print}' $ROTF | gmt plot -W0.5p,black,-
    gmt plot ${OTF} -W2p,white
    gmt plot ${OTF} -W1p,black
    gmt plot ${OTF} -Sd${RTI_ms}c -Gwhite -W0.5p,black
    if [ ! -z $plotBox ]; then 
        # fracture zone box
        gmt plot ${dataname}/input/averageBox_FZ1.lonlat -L -Gwhite@100 -W0.3p,$c_FZ  
        gmt plot ${dataname}/input/averageBox_FZ2.lonlat -L -Gwhite@100 -W0.3p,$c_FZ 
        gmt plot ${dataname}/input/averageBox_OTF1.lonlat -L -Gwhite@100 -W0.3p,$c_OTF 
    fi
    if [ ! -z $text_OTF ]; then 
        echo -42.4 30.35 ${text_OTF[0]} | gmt text -F+f9p,Helvetica-Bold,black=1p,white+a0
        echo -42.4 30.35 ${text_OTF[0]} | gmt text -F+f9p,Helvetica-Bold,black=0.2p,${c_OTF}+a0
        echo -42.4 30.2 ${text_OTF[1]} | gmt text -F+f9p,Helvetica-Bold,black=1p,white+a0
        echo -42.4 30.2 ${text_OTF[1]} | gmt text -F+f9p,Helvetica-Bold,black=0.2p,${c_OTF}+a0
    fi
    if [ ! -z "$text_FZ" ]; then 
        rot=-11
        echo -43.3 30.1 $text_FZ | gmt text -F+f9p,Helvetica-Bold,white=1p,white+a${rot}
        echo -43.3 30.1 $text_FZ | gmt text -F+f9p,Helvetica-Bold,black=0.2p,${c_FZ}+a${rot}
        echo -41.5 29.78 $text_FZ | gmt text -F+f9p,Helvetica-Bold,black=1p,white+a${rot}
        echo -41.5 29.78 $text_FZ | gmt text -F+f9p,Helvetica-Bold,black=0.2p,${c_FZ}+a${rot}
    fi
    if [ ! -z "$text_NTO" ]; then 
        rot=-9
        echo -42.8 29.25 $text_NTO | gmt text -F+f9p,Helvetica-Bold,black=1p,white+a${rot}
        echo -42.8 29.25 $text_NTO | gmt text -F+f9p,Helvetica-Bold,black=0.2p,white+a${rot}
    fi
}
function figure2
{
    gmt gmtset FONT_ANNOT_PRIMARY=7p
    gmt gmtset FONT_LABEL=7p
    dataname=MAR_Atlantis
    bathy=${dataname}/input/$dataname.nc
    bathy_ship=${dataname}/input/${dataname}_ship.nc
    bathy_large=${dataname}/input/${dataname}_large_${dataSource}.nc  
    mba=${dataname}/Results_FFT/${dataname}_mba.nc
    rmba=${dataname}/Results_FFT/${dataname}_rmba.nc
    moho=${dataname}/Results_FFT/${dataname}_moho.nc
    figwidth=7
    figheight=`gmt_get_map_height -R$rmba  -JM${figwidth}c`
    move_x=10
    move_y=10
    dx=0.5
    dy=1.5
    gmt gmtset MAP_FRAME_WIDTH=2p
    cpt_bathy=bathy.cpt
    gmt grd2cpt $bathy -Cbasecpt_bathy.cpt -Z >$cpt_bathy
    afg_grav=a20f5
    contours_grav="-40,-30,-20,-10,0,10,20,30,40"
    gmt begin Figure2 pdf
        # bathymetry
            gmt basemap -R$rmba -JM${figwidth}c -Ba -BWsen -X${move_x}c -Y${move_y}c
            gmt grdgradient $bathy -A30 -Nt1 -Qc -G${bathy}.grad
            # get min max of gravity
            # makecpt_grd $bathy
            gmt grdimage $bathy_large -C$cpt_bathy
            gmt grdimage $bathy -C$cpt_bathy -I${bathy}.grad -Q
            gmt colorbar -DJCB+o0/0.25c -Bxaf+l"Seafloor depth" -By+l"km" -C$cpt_bathy -W-0.001 --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            RTI_ms=0.2
            plotControlTransformFault_fig2 "Transform fault" "Fracture zone" "NTO"
            echo "(A)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # MBA
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -Ba -Bwsen -X${move_x}c #-Y${move_y}c
            makecpt_grd_basecpt $mba basecpt_MBA.cpt
            gmt grdimage $mba #-I${rmba}.grad
            # gmt grdcontour $mba -C$contours_grav
            gmt colorbar -DJCB+o0/0.25c -Bx${afg_grav}+l"MBA" -By+l"mGal" -G${data_min}/${data_max} --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault_fig2
            echo "(B)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # Thermal correction
            move_x=`echo "- $figwidth - $dx" | bc`
            move_y=`echo "-$figheight - $dy" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -Ba -BWSen -X${move_x}c -Y${move_y}c
            modelname=Visco-plastic
            grav_therm=${dataname}/grav_Thermal/grav_${modelname}.nc
            grav_therm_minusMeanValue=grav_therm.nc
            # mean value of the thermal gravity
            meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
            gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
            # gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
            makecpt_grd_basecpt $grav_therm_minusMeanValue basecpt_thermal.cpt
            gmt grdimage $grav_therm_minusMeanValue # -I${grav_therm_minusMeanValue}.grad
            # gmt grdcontour $grav_therm_minusMeanValue -C$contours_grav
            gmt colorbar -DJCB+o0/0.7c -Bx${afg_grav}+l"Gravity (@%12%\144@%%g@-T@-) of thermal structure" -By+l"mGal" -G${data_min}/${data_max} --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault_fig2
            echo "(C)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # RMBA
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -Ba -BwSen -X${move_x}c #-Y${move_y}c
            makecpt_grd $rmba
            gmt grdimage $rmba #-I${rmba}.grad
            # gmt grdcontour $rmba -C$contours_grav
            gmt colorbar -DJCB+o0/0.7c -Bx${afg_grav}+l"RMBA" -By+l"mGal" -G${data_min}/${data_max} --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault_fig2  "Thin crust" "Thick crust" "" "True"
            echo "(D)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
    gmt end show
    rm *.nc *.grad
}
function plotControlTransformFault_SI()
{
    plot_averageBox=$1
    plot_averageText=$2
    OTF=${dataname}/input/${dataname}_OTF.txt
    ROTF=${dataname}/input/${dataname}_ridge_otf.txt
    ridgeIndex=${dataname}/input/ridge_index.txt
    c_OTF=red 
    c_FZ=blue
    nRidge=`awk 'END{print NR}' $ridgeIndex`
    for (( i=1; i<=$nRidge; i++ ));
    do 
        start_ridge=`awk -v row=$i 'NR==row{print $1}' $ridgeIndex `
        end_ridge=`awk -v row=$i 'NR==row{print $2}' $ridgeIndex `
        awk -v row1=$start_ridge -v row2=$end_ridge 'NR>=row1 && NR<=row2 {print}' $ROTF | gmt plot -W2p,white
        awk -v row1=$start_ridge -v row2=$end_ridge 'NR>=row1 && NR<=row2 {print}' $ROTF | gmt plot -W0.5p,black,-
    done
    gmt plot ${OTF} -W2p,white
    gmt plot ${OTF} -W1p,black
    gmt plot ${OTF} -Sd${RTI_ms}c -Gwhite -W0.5p,black
    # boxnames=(averageBox_OTF averageBox_lFZ averageBox_rFZ averageBox_lOTF averageBox_rOTF)
    if [ "$plot_averageBox" == "AverageBoxes" ]; then 
        file_boxinfo=${dataname}/input/averageRMA_boxinfo.txt
        nbox=`awk 'END{print NR}' $file_boxinfo`
        for (( i=2; i<=$nbox; i++ ));
        do 
            boxname=averageBox_`awk -v row=$i 'NR==row{print $1}' $file_boxinfo`
            ind_OTF=`awk -v row=$i 'NR==row{print $2}' $file_boxinfo`
            loc_box=`awk -v row=$i 'NR==row{print $3}' $file_boxinfo`
            l_box=`awk -v row=$i 'NR==row{print $4}' $file_boxinfo`
            w_box=`awk -v row=$i 'NR==row{print $5}' $file_boxinfo`
            rot_angle=`awk -v row=$i 'NR==row{print $6}' $file_boxinfo`
            boxcolor=`awk -v row=$i 'NR==row{print $7}' $file_boxinfo`
            gmt plot ${dataname}/input/${boxname}.lonlat -W0.5p,$boxcolor -L -Gwhite@50
            if [ "$plot_averageText" == "ValueText" ]; then 
                getAverageBox_TF_FZ $boxname $ind_OTF $loc_box $l_box $w_box $rot_angle
                # plot text of mean rmba
                averageValue=`awk -v row=$i 'NR==row{print $1}' ${dataname}/Results_FFT/averageRMBA.txt`
                echo $lonClatC_box $averageValue | gmt text -F+f7p,Helvetica-Bold,$boxcolor=0.1p,black+a$angle_rot
            fi 
        done
    fi 
}
function SupplementFigure1_11
{
    dataname=$1
    gmt gmtset FONT_ANNOT_PRIMARY=7p
    gmt gmtset FONT_LABEL=7p
    # gmt gmtset MAP_ANNOT_OBLIQUE=32
    bathy=${dataname}/input/$dataname.nc
    bathy_ship=${dataname}/input/${dataname}_ship.nc
    bathy_large=${dataname}/input/${dataname}_large_${dataSource}.nc  
    mba=${dataname}/Results_FFT/${dataname}_mba.nc
    rmba=${dataname}/Results_FFT/${dataname}_rmba.nc
    moho=${dataname}/Results_FFT/${dataname}_moho.nc
    figwidth=8
    figheight=`gmt_get_map_height -R$rmba  -JM${figwidth}c`
    move_x=10
    move_y=10
    dx=0.5
    dy=1.5
    gmt gmtset MAP_FRAME_WIDTH=0.7p
    cpt_bathy=bathy.cpt
    gmt grd2cpt $bathy -Cbasecpt_bathy.cpt -Z >$cpt_bathy
    afg_grav=a20f5
    contours_grav="-40,-30,-20,-10,0,10,20,30,40"
    loc_legend=BR
    if [ "$dataname" == "Chile_Ridge" ]; then 
        loc_legend=TR 
    elif [ "$dataname" == "AtlantisII" ]; then 
        loc_legend=TR 
    elif [ "$dataname" == "MAR_Marathon" ]; then 
        loc_legend=TR 
    fi
    gmt begin SupplementFigures/SI_${dataname}  pdf
        # bathymetry
            gmt basemap -R$rmba -JM${figwidth}c -B0 -BWSen -X${move_x}c -Y${move_y}c
            gmt grdgradient $bathy -A30 -Nt1 -Qc -G${bathy}.grad
            # get min max of gravity
            # makecpt_grd $bathy
            gmt grdimage $bathy_large -C$cpt_bathy -t50
            gmt grdimage $bathy_ship -C$cpt_bathy -I${bathy}.grad -Q
            gmt colorbar -DJCB+o0/0.7c -Bxa1f0.2+l"Seafloor depth" -By+l"km" -C$cpt_bathy -W-0.001 --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            RTI_ms=0.2
            plotControlTransformFault_SI "AverageBoxes"
            gmt basemap -BWSen -Bafg -Lj${loc_legend}+o0.5c/0.5c+w50k+f+u
            echo "(A)" | gmt text -F+cTL+f10p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # MBA
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -B0 -X${move_x}c #-Y${move_y}c
            makecpt_grd_basecpt $mba basecpt_MBA.cpt
            gmt grdimage $mba #-I${rmba}.grad
            # gmt grdcontour $mba -C$contours_grav
            gmt colorbar -DJCB+o0/0.7c -Bx${afg_grav}+l"MBA" -By+l"mGal" -G${data_min}/${data_max} --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault_SI
            gmt basemap -BwseN -Bafg #-LjBR+o0.5c/0.5c+w50k+f+u
            echo "(B)" | gmt text -F+cTL+f10p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # Thermal correction
            # move_x=`echo "- $figwidth - $dx" | bc`
            move_y=`echo "-$figheight - $dy" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -B0 -BwSen -X${move_x}c #-Y${move_y}c
            modelname=Visco-plastic
            grav_therm=${dataname}/grav_Thermal/grav_${modelname}.nc
            grav_therm_minusMeanValue=grav_therm.nc
            # mean value of the thermal gravity
            meanGrav=`gmt grdmath $grav_therm MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
            gmt grdmath $grav_therm $meanGrav SUB = $grav_therm_minusMeanValue
            # gmt grdgradient $grav_therm_minusMeanValue -A30 -Nt0.6 -Qc -G${grav_therm_minusMeanValue}.grad
            makecpt_grd_basecpt $grav_therm_minusMeanValue basecpt_thermal.cpt
            gmt grdimage $grav_therm_minusMeanValue # -I${grav_therm_minusMeanValue}.grad
            # gmt grdcontour $grav_therm_minusMeanValue -C$contours_grav
            gmt colorbar -DJCB+o0/0.7c -Bx${afg_grav}+l"Gravity (@%12%\144@%%g@-T@-) of thermal structure" -By+l"mGal" -G${data_min}/${data_max} --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault_SI
            gmt basemap -BwSen -Bafg
            echo "(C)" | gmt text -F+cTL+f10p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # RMBA
            # move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -B0 -BwSEn -X${move_x}c #-Y${move_y}c
            makecpt_grd $rmba
            gmt grdimage $rmba #-I${rmba}.grad
            # gmt grdcontour $rmba -C$contours_grav
            gmt colorbar -DJCB+o0/0.7c -Bx${afg_grav}+l"RMBA" -By+l"mGal" -G${data_min}/${data_max} --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault_SI  "AverageBoxes" "ValueText"
            gmt basemap -BwsEN -Bafg
            echo "(D)" | gmt text -F+cTL+f10p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
            data_dRMBA=`awk 'NR>1{if(NR==2){rmba_FZ1=$1}else{rmba_OTF+=$1; nOTF+=1;}} END{{rmba_OTF-=$1; nOTF-=1; rmba_FZ2=$1; meanRMBA_OTF=rmba_OTF/nOTF; dRMBA=meanRMBA_OTF-(rmba_FZ1+rmba_FZ2)/2.0; err1=meanRMBA_OTF-rmba_FZ1-dRMBA;err1=sqrt(err1*err1)};print 0, dRMBA, err1}' ${dataname}/Results_FFT/averageRMBA.txt`
            echo $data_dRMBA | awk '{print "@%12%\104@%%RMBA@-TF-FZ@- = "$2" @%12%\261@%% "$3" mGal"}' | gmt text -F+c${loc_legend}+f10p,Helvetica-Bold -Dj0.05c/0.3c -Gwhite@20
    gmt end show
    rm *.nc *.grad
}
function plotBathymetry_differentTypesSpreading
{
    figwidth=7
    # figheight=`gmt_get_map_height -R$rmba  -JM${figwidth}c`
    move_x=1
    move_y=1
    
    for i in 1 5 6 10 4 ;
    do 
        dataName=${TFs[$i]}
        bathy=${dataName}/input/${dataName}_ship.nc
        cpt_bathy=bathy.cpt
        gmt grd2cpt $bathy -Cbasecpt_bathy.cpt -Z >$cpt_bathy
        # plot
        gmt begin bathymetry_diffTypesSpread/${dataName} pdf
            gmt basemap -R$bathy -JM${figwidth}c -Ba -BWsen -X${move_x}c -Y${move_y}c
            gmt grdgradient $bathy -A30 -Nt1 -Qc -G${bathy}.grad
            # get min max of gravity
            # makecpt_grd $bathy
            # gmt grdimage $bathy_large -C$cpt_bathy
            gmt grdimage $bathy -C$cpt_bathy -I${bathy}.grad -Q
            gmt colorbar -DJCB+o0/0.25c -Bxaf+l"Seafloor depth" -By+l"km" -C$cpt_bathy -W-0.001 --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
        gmt end show
    done
}
#  call functions

# 1. Prepare data
PrepareData

# # 2. transform coordinates
# waterDepth
# transformCoordinates ${OTF} 0 0 `echo $meanWaterDepth $z_seafloor_ASEPCTmodel | awk '{printf "%.0f", sqrt($1*$1)+$2}'`

# 2. calculate MBA
# FAA2MBA_FFT
# # # 2.2 plot profiles
# # # plotProfiles input/${dataname}_faa.nc FAA
# plotProfiles ${dataname}/Results_FFT/${dataname}_mba.nc MBA FFT "MBA profiles: FFT approach, -W${W}"
# plotProfiles ${dataname}/Results_FFT/${dataname}_rmba.nc RMBA FFT "RMBA profiles: FFT approach, -W${W}"

# # 2.2 plot ASPECT result slice
# xy2lonlat_vtu_ASPECT ${OTF} 0 0 ${dataname}/slice/slice_stressxx_top

# # 3. RMBA to Moho
# RMBA2Moho 1.5
# RMBA2Moho_v2 0 800
# plotMohoInversion ${dataname}/Results_FFT/tmp_Moho/Moho_800.nc

# 4. extract and calculate deltaRMBA
# calAverageRMBA
# 4.1 plot delta RMBA between TF and FZ
# figure3
# figure3_v2
# figure3_v2 ErrorBar
# figure3_dRMBA_T2R ErrorBar
# figure1
# figure2
# for i in {0..10};
# do 
#     dataname=${TFs[$i]}
#     SupplementFigure1_11 $dataname
# done
# SupplementFigure1_11 MarieCeleste

# # 5. get slice from thermal model result (vtu file), ad save to .csv file
# getSlice_thermalModel -130.7 44.6 -128.75 44.06

# # 5.1 get temperature of earthquake epercenters
# extractProperties_thermalModel_lonlat ${dataname}/input/cata No_isostatic
# extractProperties_thermalModel_lonlat ${dataname}/input/cata isostatic
# extractProperties_thermalModel_lonlat ${dataname}/input/cata IsoVisco
# extractProperties_thermalModel_lonlat ${dataname}/input/cata Visco-plastic

# FAA2MBA_Spatial
# plotProfiles ${dataname}/Results_Spatial/${dataname}_mba.nc MBA Spatial "MBA profiles: Spatial domain approach"
# plotProfiles ${dataname}/Results_Spatial/${dataname}_rmba.nc RMBA Spatial "RMBA profiles: Spatial domain approach"


# # 3. convert gravity of thermal effect to grd format
# gravity_xyz2grd grav_Thermal/grav_Viscous-plastic

# 4. difference of gravtiy between different geodynamic models
# diff_Gravity_models

# plotBathymetry_differentTypesSpreading

rm gmt.history