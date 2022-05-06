globalDataPath=~/MyData/DATA/GlobalGeophysicalData
ETOPO1=${globalDataPath}/etopo/ETOPO1_Bed_g_gdal.nc
dataname=Kane #DuTroit,Bouvet,Romanche, Kane, EPR_Discovery, Chile_Ridge, MAR_Marathon, SWIR_Marion, Blanco, Oceanographer_all Oceanographer, Gofar, Garrett, Vlamingh, MAR_Atlantis, Clipperton, MarieCeleste, AtlantisII
TFs=(    MAR_Marathon SWIR_Marion AtlantisII MarieCeleste Oceanographer MAR_Atlantis Vlamingh Chile_Ridge Garrett Gofar Clipperton)
TFNames=(Marathon     Marion      AtlantisII "Marie Celeste"           Oceanographer Atlantis     Vlamingh  "CR-TF 39@%12%\260@%%S"          Garrett Gofar Clipperton)
Ridges=(MAR           SWIR        SWIR       CIR                       MAR            MAR         SEIR       CR                                   EPR    EPR   EPR)
Ridge_index=(2        0           0          1                          2              2            4         3                                     5     5     5)
Ridge_names=(SWIR CIR MAR CR SEIR EPR)
symbols_ridges=(c d i s t a)
color_ridges=(magenta   red  blue         green4 black)
type_ridges=(Ultraslow Slow Intermediate Fast  Other)
bathy=${dataname}/input/$dataname.nc
bathy_ship=${dataname}/input/${dataname}_ship.nc

dataSource=ship # Sat, ship
bathy_large=${dataname}/input/${dataname}_large_${dataSource}.nc  
lon_min_large=`gmt grdinfo $bathy | grep "x_min" | awk '{print $3-0.2}'`
lon_max_large=`gmt grdinfo $bathy | grep "x_min" | awk '{print $5+0.2}'`
lat_min_large=`gmt grdinfo $bathy | grep "y_min" | awk '{print $3-0.2}'`
lat_max_large=`gmt grdinfo $bathy | grep "y_min" | awk '{print $5+0.2}'`
range_large=${lon_min_large}/${lon_max_large}/${lat_min_large}/${lat_max_large}
# echo $lon_max_large $lat_max_large | cs2cs +proj=latlon +to +proj=tmerc +datum=potsdam +lon_0=$lon_min_large +lat_0=$lat_min_large +x_0=0 +y_0=0 -f "%.2E" | awk '{print "xy长度: "$1/1000, $2/1000, "对角线长度：", sqrt($1*$1+$2*$2)/1000}'

OTF=${dataname}/input/${dataname}_OTF.txt
ROTF=${dataname}/input/${dataname}_ridge_otf.txt
faa=${dataname}/input/${dataname}_faa.nc
rho_crust=2700
rho_water=1020
rho_mantle=3300
crust_thickness=6000
pi=3.141592653
z_seafloor_ASEPCTmodel=100000

G=6.67E-6 #will generate gravity in unit mGal
W=0 # used in gravfft
function waterDepth()
{
    meanWaterDepth=`gmt grdmath $bathy MEAN = tmp.nc && gmt grdinfo tmp.nc | grep "v_min" | awk '{print $3}' && rm tmp.nc`
    echo "Mean water depth" $meanWaterDepth
    meanMohoDepth=`echo $meanWaterDepth $crust_thickness | awk '{print $1-$2}'`
    g1=`echo $pi $G $rho_water $rho_crust $meanWaterDepth | awk '{print -2*$1*$2*($3-$4)*$5}'`
    g2=`echo $pi $G $rho_crust $rho_mantle $meanMohoDepth | awk '{print -2*$1*$2*($3-$4)*$5}'`
    # gravity difference between FFT method and "Spatial domain method"
    gconst=`echo $g1 $g2 | awk '{print $1+$2}'`
    # echo "g1="$g1 "g2="$g2 "g1+g2="$gconst
}
masterCPT_grav=${dataname}/input/basecpt_grav.cpt
basecpt_moho=basecpt_moho.cpt

# gmt makecpt -T-40/60/1 -Z >mba.cpt
# gmt makecpt -T-45/90/1 -Z >rmba.cpt
# gmt makecpt -T-5/8/1 -Z >grav_diff.cpt

if [ ! -d ${dataname} ]; then 
    mkdir ${dataname}
fi
if [ ! -d ${dataname}/input ]; then 
    mkdir ${dataname}/input
    touch ${dataname}/input/${dataname}_OTF.txt 
    touch ${dataname}/input/${dataname}_ridge_otf.txt 
    touch ${dataname}/input/${dataname}_ridgeN.txt 
    touch ${dataname}/input/${dataname}_ridgeS.txt 
    cp basecpt_grav.cpt ${dataname}/input
fi
if [ ! -d ${dataname}/grav_Thermal ]; then 
    mkdir ${dataname}/grav_Thermal
fi
if [ ! -d ${dataname}/Results_FFT ]; then 
    mkdir ${dataname}/Results_FFT
fi
# if [ ! -d ${dataname}/Results_Spatial ]; then 
#     mkdir ${dataname}/Results_Spatial
# fi