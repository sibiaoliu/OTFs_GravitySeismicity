#===== Script for RMBA to Crustal thickness=====
# 15.04.2022
# Generate the rmba-column, nx, ny, longx, longy
#===== Input parameters' value change=====
process=$1
dataname=(Clipperton Discovery Garrett Gofar)
etas=(hsc isov disl vp vep)

function rmbatxt()
{
    for i in {0..3}; do
        for j in {0..4}; do
            path_rmba=${dataname[i]}/no-hydroT-cooling/Results/${dataname[i]}_fullrmba_${etas[j]}
            gmt grd2xyz ${path_rmba}.nc > ${path_rmba}.xyz
            awk '{print $3}' ${path_rmba}.xyz > ${dataname[i]}/no-hydroT-cooling/Results/${dataname[i]}_fullrmba_${etas[j]}.txt
            #rm ${path_rmba}.xyz            
        done
    done
}
function mohonc()
{
    for i in {0..3}; do
        for j in {0..4}; do           
            path_rmba=${dataname[i]}/no-hydroT-cooling/Results/${dataname[i]}_fullrmba_${etas[j]}
            path_moho=${dataname[i]}/no-hydroT-cooling/Results/${dataname[i]}_moho_${etas[j]}
            path_mohoinv=${dataname[i]}/no-hydroT-cooling/Results/${dataname[i]}_mohoinv_${etas[j]}
            paste ${path_rmba}.xyz ${path_moho}.txt  | awk '{print $1,$2,$4}' >${path_moho}.lonlat
            gmt xyz2grd ${path_moho}.lonlat -R${path_rmba}.nc -G${path_moho}.nc
            paste ${path_rmba}.xyz ${path_mohoinv}.txt  | awk '{print $1,$2,$4}' >${path_mohoinv}.lonlat
            gmt xyz2grd ${path_mohoinv}.lonlat -R${path_rmba}.nc -G${path_mohoinv}.nc
        done
    done            
}
function deletetmp()
{
    for i in {0..3}; do
        for j in {0..4}; do           
            path_rmba=${dataname[i]}/no-hydroT-cooling/Results/${dataname[i]}_fullrmba_${etas[j]}
            path_moho=${dataname[i]}/no-hydroT-cooling/Results/${dataname[i]}_moho_${etas[j]}
            path_mohoinv=${dataname[i]}/no-hydroT-cooling/Results/${dataname[i]}_mohoinv_${etas[j]}
            rm ${path_moho}.lonlat ${path_mohoinv}.lonlat ${path_rmba}.xyz 
            rm ${path_moho}.txt ${path_mohoinv}.txt ${path_rmba}.txt
        done
    done            
}
function help()
{
    echo "-) bash rmba_moho.sh rmbatxt"
    echo "-) bash rmba_moho.sh mohonc"
    echo "-) bash rmba_moho.sh deletetmp"
    exit
}
       
if [ "$process" == "rmbatxt" ]; then 
    echo "Save one column of RMBA"
    rmbatxt
elif [ "$process" == "mohonc" ]; then 
    echo "Paste mono column with lon/lat"
    mohonc
elif [ "$process" == "deletetmp" ]; then 
    echo "Delete all tmp files"
    deletetmp       
else
    help
fi