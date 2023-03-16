#===== This script is used to plot the supplementary figures=====
# 03.2023---Sibiao Liu
. gmt_shell_functions.sh

#===== Input parameters' value change=====
    dataname=Pitman
    bathy_ship=../inputs/${dataname}_ship.nc
    bathy_large_ship=../inputs/${dataname}_large_ship.nc
    OTF=../inputs/${dataname}_OTF.txt
    ROTF=../inputs/${dataname}_ridge_otf.txt
    mba=Results/${dataname}_mba.nc
    faa=../inputs/${dataname}_faa.nc
    
    cpt_bathy=bathy.cpt
    gmt grd2cpt $bathy_ship -C../../../basecpt_bathy.cpt -Z >$cpt_bathy
    gmt makecpt -C../../../romaO.cpt+h -I -T-10/5/1 -Z >grav_diff.cpt
    gmt makecpt -C../../../romaO.cpt+h -I -T-5/10/1 -Z >grav_diff2.cpt
    gmt makecpt -C../../../romaO.cpt+h -I -T-50/35/20 -Z >grav_faa.cpt
    gmt makecpt -C../../../basecpt_grav.cpt+h -T-50/50/20 -Z >grav_rmba.cpt


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
        #gmt plot ${file_ROTF} -W1p,white
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
function makecpt_grd_basecpt()
{
    grdfile=$1
    basecpt=$2
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
function shiftDatabyMean()
{
    grd_data=$1
    grd_data_shift=$2
    meanValue=`gmt grdmath ${grd_data} MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
    gmt grdmath ${grd_data} ${meanValue} SUB = ${grd_data_shift}
    str_shiftData=", shift ${meanValue} mGal"
}
function gravity_xyz2grd()
{
    gravityFile=$1
    gmt grd2xyz ${faa} >${dataname}_faa.xyz
    paste ${dataname}_faa.xyz ${gravityFile}.txt  | awk '{print $1,$2,$4}' >${gravityFile}.lonlat
    gmt xyz2grd ${gravityFile}.lonlat -R$faa -G${gravityFile}.nc
    rm ${gravityFile}.lonlat ${dataname}_faa.xyz
}
# Plot the figure, including
    # 0.0 Bathymetry, 0.1 FAA, 0.2 RMBA_nlvp
    # 1.0 RMBA_isov-hsc; 0.2 RMBA_nlvp-hsc; 0.3 RMBA_nlvp-isov
    Etas=(hsc isov vp)
    m2=${Etas[2]}
    gravity_xyz2grd grav_$m2
    grav_therm=grav_$m2.nc
    fullrmba_vp=${dataname}_fullrmba_$m2.nc 
    gmt grdmath $mba $grav_therm SUB = $fullrmba_vp
    shiftDatabyMean $fullrmba_vp $fullrmba_vp
    shiftRMBA=${meanValue}

    # Plot range
    lon_min_small=`gmt grdinfo $fullrmba_vp | grep "x_max" | awk '{print $3}'`
    lon_max_small=`gmt grdinfo $fullrmba_vp | grep "x_max" | awk '{print $5}'`
    lat_min_small=`gmt grdinfo $fullrmba_vp | grep "y_max" | awk '{print $3}'`
    lat_max_small=`gmt grdinfo $fullrmba_vp | grep "y_max" | awk '{print $5}'`
    range_small=${lon_min_small}/${lon_max_small}/${lat_min_small}/${lat_max_small}

    gmt gmtset FONT_ANNOT_PRIMARY=7p
    gmt gmtset FONT_LABEL=7p
    figwidth=5
    figheight=`gmt_get_map_height -R$range_small  -JM${figwidth}c`
    move_x=10
    move_y=10
    dx=0.5
    dy=1.7
    gmt gmtset MAP_FRAME_WIDTH=2p
    afg_grav=a10f5

    gmt begin ${dataname}_rmba_supplement pdf
        # Bathymetry
            gmt basemap -R$range_small -JM${figwidth}c -Ba -BWSen -X${move_x}c -Y${move_y}c
            gmt grdgradient $bathy_large_ship -A30 -Nt1 -Qc -G${bathy_large_ship}.grad
            gmt grdimage $bathy_large_ship -C$cpt_bathy -I${bathy_large_ship}.grad -Q
            gmt colorbar -DJCB+o0/0.6c -Bxaf+l"Depth" -By+l"km" -C$cpt_bathy -W-0.001 --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(a)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # FAA
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$range_small -JM${figwidth}c -Ba -Bwsen -X${move_x}c #-Y${move_y}c
            gmt grdimage $faa -Cgrav_faa.cpt
            gmt colorbar -DJCB+o0/0.6c -Bx${afg_grav}+l"Free-air anomaly" -By+l"mGal" -Cgrav_faa.cpt --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(b)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite         
        # RMBA
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$range_small -JM${figwidth}c -Ba -Bwsen -X${move_x}c #-Y${move_y}c
            gmt grdimage $fullrmba_vp -Q -Cgrav_rmba.cpt #-Imbashift.grad
            gmt colorbar -DJCB+o0/0.6c -Bx${afg_grav}+l"RMBA, shift $shiftRMBA" -By+l"mGal" -Cgrav_rmba.cpt --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(c)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # dRMBA:isov-hsc
            move_x=`echo "-$figwidth - $dx -$figwidth - $dx" | bc`
            move_y=`echo "-$figheight - $dy" | bc`
            gmt basemap -R$range_small -JM${figwidth}c -Ba -BWSen -X${move_x}c -Y${move_y}c

            m1=${Etas[0]}
            m2=${Etas[1]}
            rmba1=Results/${dataname}_fullrmba_$m1.nc
            rmba2=Results/${dataname}_fullrmba_$m2.nc
            grav_diff=fullrmba_diff_$m2-$m1.nc
            gmt grdmath $rmba2 $rmba1 SUB = ${grav_diff}
            gmt grdimage $grav_diff -Q -Cgrav_diff.cpt
            gmt colorbar -DJCB+o0/0.6c -Bxaf+l"RMBA difference: ISOV - HSC" -By+l"mGal" -Cgrav_diff.cpt --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(d)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # dRMBA: nlvp-hsc
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$range_small -JM${figwidth}c -Ba -Bwsen -X${move_x}c #-Y${move_y}c
            m1=${Etas[0]}
            m2=${Etas[2]}
            rmba1=Results/${dataname}_fullrmba_$m1.nc
            rmba2=Results/${dataname}_fullrmba_$m2.nc
            grav_diff=fullrmba_diff_$m2-$m1.nc
            gmt grdmath $rmba2 $rmba1 SUB = ${grav_diff}
            gmt grdimage $grav_diff -Q -Cgrav_diff.cpt
            gmt colorbar -DJCB+o0/0.6c -Bxaf+l"RMBA difference: NLVP - HSC" -By+l"mGal" -Cgrav_diff.cpt --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(e)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # dRMBA: nlvp-isov
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$range_small -JM${figwidth}c -Ba -Bwsen -X${move_x}c #-Y${move_y}c
            m1=${Etas[1]}
            m2=${Etas[2]}
            rmba1=Results/${dataname}_fullrmba_$m1.nc
            rmba2=Results/${dataname}_fullrmba_$m2.nc
            grav_diff=fullrmba_diff_$m2-$m1.nc
            gmt grdmath $rmba2 $rmba1 SUB = ${grav_diff}
            gmt grdimage $grav_diff -Q -Cgrav_diff2.cpt
            gmt colorbar -DJCB+o0/0.6c -Bxaf+l"RMBA difference: NLVP - ISOV" -By+l"mGal" -Cgrav_diff2.cpt --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(f)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
    gmt end show
rm gmt.* *.nc *.cpt *.grad

# # Plot the figure, including
#     # 1.0 RMBA_isov-hsc; 0.2 RMBA_nlvp-hsc; 0.3 RMBA_nlvp-isov
#     Etas=(hsc isov vp)
#     gmt begin Results/${dataname}_supplement pdf
#       gmt subplot begin 2x3 -Fs16.5c/8.5c -M0.5c/1.3c -R$faa -JM15c -Ba1df30m -BWSen
#         # 0.0 Bathymetry, 0.1 FAA, 0.2 RMBA_nlvp
#         gmt subplot set 0,0
#             gmt grdgradient $bathy_large_ship -A30 -Nt0.6 -Qc -G${bathy_ship}.grad 
#             gmt grdimage $bathy_large_ship -BWseN -C$cpt_bathy -I${bathy_ship}.grad -Q 
#             gmt colorbar -DJCB+o0/0.6c -Bxaf+l"Bathymetry" -By+l"km" -C$cpt_bathy
#             echo "(a)" | gmt text -F+cTL+f15p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite            
#             plotControlTransformFault
#         gmt subplot set 0,1
#             gmt grdimage $faa -Cgrav_faa.cpt -BWseN
#             gmt colorbar -DJCB+o0/0.6c -Bx${afg_grav}+l"Free-air anomaly" -By+l"mGal" -Cgrav_faa.cpt
#             echo "(b)" | gmt text -F+cTL+f15p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
#             plotControlTransformFault
#         gmt subplot set 0,2
#             m1=${Etas[2]}
#             gravity_xyz2grd grav_$m1
#             grav_therm=grav_$m1.nc
#             fullrmba_vp=${dataname}_fullrmba_$m1.nc 
#             gmt grdmath $mba $grav_therm SUB = $fullrmba_vp
#             #shiftDatabyMean $fullrmba_vp
#             #shiftRMBA=${meanValue}
#             gmt grdimage $fullrmba_vp -BwseN -Cgrav_rmba.cpt
#             gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA" -By+l"mGal" -Cgrav_rmba.cpt            
#             #gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA_$m1 $shiftRMBA" -By+l"mGal" -Cgrav_rmba.cpt
#             echo "(c)" | gmt text -F+cTL+f15p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
#             plotControlTransformFault
#         gmt subplot set 1,0
#             m1=${Etas[0]}
#             m2=${Etas[1]}
#             rmba1=Results/${dataname}_rmba_$m1.nc
#             rmba2=Results/${dataname}_rmba_$m2.nc
#             grav_diff=rmba_diff_$m2-$m1.nc
#             gmt grdmath $rmba2 $rmba1 SUB = ${grav_diff}
#             gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
#             gmt grdimage $grav_diff -BWseN -Cgrav_diff.cpt 
#             gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: ISOV - HSC" -By+l"mGal" -Cgrav_diff.cpt 
#             echo "(d)" | gmt text -F+cTL+f15p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
#             plotControlTransformFault
#         gmt subplot set 1,1
#             m1=${Etas[0]}
#             m2=${Etas[2]}
#             rmba1=Results/${dataname}_rmba_$m1.nc
#             rmba2=Results/${dataname}_rmba_$m2.nc
#             grav_diff=rmba_diff_$m2-$m1.nc
#             gmt grdmath $rmba2 $rmba1 SUB = ${grav_diff}
#             gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff}.grad
#             gmt grdimage $grav_diff -BWseN -Cgrav_diff.cpt 
#             gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: NLVP - HSC" -By+l"mGal" -Cgrav_diff.cpt 
#             echo "(e)" | gmt text -F+cTL+f15p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
#             plotControlTransformFault
#         gmt subplot set 1,2
#             m1=${Etas[1]}
#             m2=${Etas[2]}
#             rmba1=Results/${dataname}_rmba_$m1.nc
#             rmba2=Results/${dataname}_rmba_$m2.nc
#             grav_diff=rmba_diff_$m2-$m1.nc
#             gmt grdmath $rmba2 $rmba1 SUB = ${grav_diff}
#             gmt grdgradient $grav_diff -A30 -Nt0.6 -Qc -G${grav_diff2}.grad
#             gmt grdimage $grav_diff -BWseN -Cgrav_diff2.cpt 
#             gmt colorbar -DJCB+o0/0.5c -Bxaf+l"RMBA difference: NLVP - ISOV" -By+l"mGal" -Cgrav_diff2.cpt 
#             echo "(f)" | gmt text -F+cTL+f15p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
#             plotControlTransformFault
#         gmt subplot end
#     gmt end show
# rm gmt.* *.nc *.cpt *.grad

