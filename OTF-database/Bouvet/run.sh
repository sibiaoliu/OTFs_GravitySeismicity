process=$1
dataname=${PWD##*/}
datasource=ship

siteFile_thermal=input/MAR_Marathon_sites_103374.xyz
outPath_thermal=grav_Thermal
if [ ! -f $siteFile_thermal ] ; then 
    echo "File does not exist: $siteFile_thermal"
    exit
fi
function thermal()
{
    thermalmodel=Visco-plastic
    vtu_rho=../../geodynamic-models_3D/MAR_Marathon/Visco-plastic/Marathon-vp-8.51Ma.vtu
    vtk2grav -p $siteFile_thermal -F density -D 3300 -i $vtu_rho -o ${outPath_thermal}/grav_${thermalmodel}.txt
}
function help()
{
    echo "-) bash run.sh waterCrust"
    echo "-) bash run.sh thermal"
    echo "-) bash run.sh download"
    exit
}
function download()
{
    path_nevado="nevado:/home/zguo/programing/otf_grevemeyer_2021/gravity/TransformFaults"
    scp ${path_nevado}/${dataname}/$outPath_thermal/*.txt $outPath_thermal
}
if [ "$process" == "waterCrust" ]; then 
    echo "Water crust layer correction"
elif [ "$process" == "download" ]; then 
    download
elif [ "$process" == "thermal" ]; then 
    echo "thermal correction calculation"
    thermal
else
    help
fi