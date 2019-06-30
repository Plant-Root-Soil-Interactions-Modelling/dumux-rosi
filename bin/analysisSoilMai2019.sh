#extract Variables
mkdir extract
extractSoilResult.py $1-soil.pvd
PVDtoHF5.py ex_$1-soil.pvd

mv ex_$1-soil.h5 ./extract/

cd extract
analysisSoilMai2019.py ex_$1-soil.h5
