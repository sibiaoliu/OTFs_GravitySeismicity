#===== This script is used to plot the final figure=====
. gmt_shell_functions.sh

#===== Input parameters' value change=====
    dataname=Bouvet
    etaname=vp
    bathy_ship=../inputs/${dataname}_ship.nc
    OTF=../inputs/${dataname}_OTF.txt
    ROTF=../inputs/${dataname}_ridge_otf.txt
    mba=Results/${dataname}_mba.nc
    rmba=Results/${dataname}_rmba_${etaname}.nc
    moho=Results/${dataname}_moho_${etaname}.nc
    gmt makecpt -C../../../basecpt_MBA.cpt+h -T-120/120/20 -Z >grav_mba.cpt
    gmt makecpt -C../../../basecpt_grav.cpt+h -T-100/100/20 -Z >grav_rmba.cpt
    gmt makecpt -C../../../romaO.cpt+h -T-2/2/0.5 -Z >grav_moho.cpt

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
    # 1. Bathymetry, 2. MBA, 3. RMBA 4. Moho
    gmt gmtset FONT_ANNOT_PRIMARY=7p
    gmt gmtset FONT_LABEL=7p
    figwidth=7
    figheight=`gmt_get_map_height -R$rmba  -JM${figwidth}c`
    move_x=10
    move_y=10
    dx=0.5
    dy=2
    gmt gmtset MAP_FRAME_WIDTH=2p
    cpt_bathy=bathy.cpt
    gmt grd2cpt $bathy_ship -C../../../basecpt_bathy.cpt -Z >$cpt_bathy
    cpt_moho=moho.cpt
    gmt grd2cpt $moho -C../../../romaO.cpt -Z >$cpt_moho
    afg_grav=a10f5
    contours_grav="-40,-30,-20,-10,0,10,20,30,40"
    grid_mask=tmp_mask.nc
    alpha_mask=50
    gmt grdmath $bathy_ship $bathy_ship NAN = NAN.nc 
    gmt grdmath $bathy_ship NAN.nc XOR 99999 ADD = $grid_mask && rm NAN.nc
    gmt begin Figure pdf
        # bathymetry
            gmt basemap -R$rmba -JM${figwidth}c -Ba -BWSen -X${move_x}c -Y${move_y}c
            gmt grdgradient $bathy_ship -A30 -Nt1 -Qc -G${bathy_ship}.grad
            # makecpt_grd $bathy
            #gmt grdimage $bathy_large_ship -C$cpt_bathy
            gmt grdimage $bathy_ship -C$cpt_bathy -I${bathy_ship}.grad -Q
            if [ -f $grid_mask ]; then 
                gmt grdimage $grid_mask -Cbathy.cpt -Q -t$alpha_mask
            fi
            gmt colorbar -DJCB+o0/0.6c -Bxaf+l"Bathymetry" -By+l"km" -C$cpt_bathy -W-0.001 --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(A)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # MBA
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -Ba -Bwsen -X${move_x}c #-Y${move_y}c
            shiftDatabyMean $mba mbashift.nc
            #gmt grdgradient mbashift.nc -A10 -Nt1 -Qc -Gmbashift.grad
            gmt grdimage mbashift.nc -Q -Cgrav_mba.cpt #-Imbashift.grad
            if [ -f $grid_mask ]; then 
                gmt grdimage $grid_mask -Cbathy.cpt -Q # -t$alpha_mask
            fi
            gmt colorbar -DJCB+o0/0.6c -Bx${afg_grav}+l"MBA" -By+l"mGal" -Cgrav_mba.cpt --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(B)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # RMBA
            move_x=`echo "- $figwidth - $dx" | bc`
            move_y=`echo "-$figheight - $dy" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -Ba -Bwsen -X${move_x}c -Y${move_y}c
            gmt grdimage $rmba -Cgrav_rmba.cpt -Q
            if [ -f $grid_mask ]; then 
                gmt grdimage $grid_mask -Cbathy.cpt -Q # -t$alpha_mask
            fi
            gmt colorbar -DJCB+o0/0.6c -Bx${afg_grav}+l"RMBA" -By+l"mGal" -Cgrav_rmba.cpt --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(C)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
        # Moho            
            move_x=`echo "$figwidth + $dx" | bc`
            gmt basemap -R$rmba -JM${figwidth}c -Ba -Bwsen -X${move_x}c #-Y${move_y}c
            gmt grdgradient $moho -A30 -Nt0.6 -Qc -G$moho.grad 
            gmt grdimage $moho -Cgrav_moho.cpt -Q
            if [ -f $grid_mask ]; then 
                gmt grdimage $grid_mask -Cbathy.cpt -Q # -t$alpha_mask
            fi
            gmt colorbar -DJCB+o0/0.6c -Bx${afg_grav}+l"Relative crustal thickness" -By+l"km" -Cgrav_moho.cpt --MAP_FRAME_PEN=0.5p --MAP_TICK_LENGTH_PRIMARY=1p
            plotControlTransformFault $ROTF $OTF
            echo "(D)" | gmt text -F+cTL+f12p,Helvetica-Bold -Dj0.1c/0.1 -Gwhite
    gmt end show
rm gmt.* tmp_mask.nc mbashift.nc *.cpt Results/*.grad