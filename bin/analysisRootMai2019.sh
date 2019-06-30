#extract Variables
mkdir extract
extractRootResult.py $1-root.pvd
PVDtoHF5.py ex_$1-root.pvd

mv ex_$1-root.h5 ./extract/

cd extract
analysisUptakeMai2019.py ex_$1-root.h5

