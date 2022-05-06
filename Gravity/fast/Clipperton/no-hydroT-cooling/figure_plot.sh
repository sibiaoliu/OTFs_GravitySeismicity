#===== This script is used to plot the final figure=====
. gmt_shell_functions.sh

#===== Input parameters' value change=====
    dataname=Clipperton
    bathy_ship=../inputs/${dataname}_ship.nc
    OTF=../inputs/${dataname}_OTF.txt
    ROTF=../inputs/${dataname}_ridge_otf.txt
    mba=Results/${dataname}_mba.nc
    rmba=Results/${dataname}_rmba_vp.nc
    moho=Results/${dataname}_moho_vp.nc
    gmt makecpt -C../../../basecpt_grav.cpt+h -T-30/30/10 -Z >grav_rmba.cpt

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
function makecpt_grd_moho()
{
    grdfile=$1
    masterCPT_grav=$2 #../../../basecpt_grav.cpt
    data_min=`gmt grdinfo $grdfile | grep "v_min" | awk '{printf "%.1f", $3}'`
    data_max=`gmt grdinfo $grdfile | grep "v_max" | awk '{printf "%.1f", $5}'`
    cpt_min=`echo ${data_min} | awk '{printf "%.1f", $1}'`
    cpt_max=`echo ${data_max} | awk '{printf "%.1f", $1}'`
    gmt makecpt -Iz -C${masterCPT_grav}+h -T-3.0/3.0
}
function shiftDatabyMean()
{
    grd_data=$1
    grd_data_shift=$2
    meanValue=`gmt grdmath ${grd_data} MEAN = tmp.nc && gmt grd2xyz tmp.nc | awk 'NR==1{printf "%.0f", $3}' && rm tmp.nc`
    gmt grdmath ${grd_data} ${meanValue} SUB = ${grd_data_shift}
    str_shiftData=", shift ${meanValue} mGal"
}

# Plot the figure, including
    # 1. Bathymetry, 2. MBA, 3. RMBA_vp 4. Moho_vp
    gmt gmtset FONT_ANNOT_PRIMARY=7p
    gmt gmtset FONT_LABEL=7p
    figwidth=7
    figheight=`gmt_get_map_height -R$rmba  -JM${figwidth}c`
    move_x=10
    move_y=10
    dx=0.5
    dy=1.5
    gmt gmtset MAP_FRAME_WIDTH=2p
    cpt_bathy=bathy.cpt
    gmt grd2cpt $bathy_ship -C../../../basecpt_bathy.cpt -Z >$cpt_bathy
    afg_grav=a20f5
    contours_grav="-40,-30,-20,-10,0,10,20,30,40"
    gmt begin Figure pdf
        # bathymetry
            gmt basemap -R$rmba -JM${figwidth}c -Ba -BWsen -X${move_x}c -Y${move_y}c
            gmt grdgradient $bathy_ship -A30 -Nt1 -Qc -G${bathy_ship}.grad
            # makecpt_grd $bathy
            #gmt grdimage $bathy_large_ship -C$cpt_bathy
            gmt grdimage $bathy_ship -C$cpt_bathy -I${bathy_ship}.grad -Q
            gmt colorbar -DJCB+o0/0.3c -Bxaf+l"Bathymetry" -By+l"km" -C$cpt_bathy -W-0.001 --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(A)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # MBA
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -Ba -BWSen -X${move_x}c #-Y${move_y}c
            shiftDatabyMean $mba ${mba}.shift
            makecpt_grd_basecpt ${mba}.shift ../../../basecpt_MBA.cpt
            gmt grdimage ${mba}.shift -Q
            gmt colorbar -DJCB+o0/0.3c -Bx${afg_grav}+l"MBA" -By+l"mGal" -G${data_min}/${data_max} --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(B)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # RMBA
            move_x=`echo "- $figwidth - $dx" | bc`
            move_y=`echo "-$figheight - $dy" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -Ba -BwSen -X${move_x}c -Y${move_y}c
            gmt grdimage $rmba -Cgrav_rmba.cpt
            # gmt grdcontour $rmba -C$contours_grav
            gmt colorbar -DJCB+o0/0.7c -Bx${afg_grav}+l"RMBA" -By+l"mGal" -Cgrav_rmba.cpt --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(C)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # Moho            
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -Ba -Bwsen -X${move_x}c #-Y${move_y}c
            gmt grdgradient $moho -A30 -Nt0.6 -Qc -G$moho.grad 
            makecpt_grd_moho $moho ../../../roma.cpt
            gmt grdimage $moho #-I$moho.grad
            gmt colorbar -DJCB+o0/0.7c -Bx${afg_grav}+l"Relative crustal thickness" -By+l"km" -G-3/3 --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(D)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
    gmt end show
rm gmt.history bathy.cpt Results/*.shift Results/*.grad