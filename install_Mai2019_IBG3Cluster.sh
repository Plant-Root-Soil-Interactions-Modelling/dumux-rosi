# One click install script dumux

# make a new folder containing everything

echo "*********************************************************************************************"
echo "(1/2) Cloning repositories. This may take a while. Make sure to be connected to the internet."
echo "*********************************************************************************************"
# the core modules
for MOD in common geometry grid localfunctions istl; do
    if [ ! -d "dune-$MOD" ]; then
        git clone -b releases/2.5 https://gitlab.dune-project.org/core/dune-$MOD.git
    else
        echo "Skip cloning dune-$MOD because the folder already exists."
        cd dune-$MOD
        git checkout releases/2.5
        cd ..
    fi
done

# uggrid
git clone -b releases/2.5 https://gitlab.dune-project.org/staging/dune-uggrid.git

# foamgrid
if [ ! -d "dune-foamgrid" ]; then
    git clone -b releases/2.5 https://gitlab.dune-project.org/extensions/dune-foamgrid.git
else
    echo "Skip cloning dune-foamgrid because the folder already exists."
    cd dune-foamgrid
    git checkout releases/2.5
    cd ..
fi

# dumux
if [ ! -d "dumux" ]; then
    git clone -b releases/2.12 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git
else
    echo "Skip cloning dumux because the folder already exists."
    cd dumux
    git checkout releases/2.12
    cd ..
fi

if [ $? -ne 0 ]; then
    echo "*********************************************************************************************"
    echo "Failed to clone the repositories."
    echo "*********************************************************************************************"
    exit $?
fi

\cp dumux-rosi/dumux/_/common/timemanager.hh dumux/dumux/common/
\cp dumux-rosi/dumux/_/implicit/model.hh dumux/dumux/implicit/
\cp dumux-rosi/dumux/_/implicit/problem.hh dumux/dumux/implicit/
\cp dumux-rosi/dumux/_/io/gridcreator.hh dumux/dumux/io/
\cp dumux-rosi/dumux/_/io/restart.hh dumux/dumux/io/

echo "*********************************************************************************************"
echo "(2/2) Configure dune modules and dumux. Build the dune libaries. This may take several minutes."
echo "*********************************************************************************************"
# run build
## TODO: UPDATE LINK
./dune-common/bin/dunecontrol --opts=dumux-rosi/optim_cluster.opts --builddir=$(pwd)/build-pubMai2019  all
#
if [ $? -ne 0 ]; then
    echo "*********************************************************************************************"
    echo "Failed to build the dune libaries."
    echo "*********************************************************************************************"
    exit $?
fi

# echo result
echo "*********************************************************************************************"
echo "Succesfully configured and built dune, dumux and ROSI"
echo "*********************************************************************************************"

export PATH="$PATH:$PWD/dumux-rosi/bin"
#correct python path
var=$(which python2)
var2="#!$var"
cd $PWD/dumux-rosi/bin
for FILE in extractRootResultV2.py PVDtoHF5_.py extractSoilResult.py filterHF5byLayersSoilTopBottom_.py filterHF5_rootTipsUptakeMai2019V4FixedY.py HF5Data.py PVDtoHF5.py; do
    sed -i '1d' $FILE
    echo $var2 | cat - $FILE > temp && mv temp $FILE
done
cd ../..

cd build-pubMai2019/dumux-rosi/rosi_multiscale/rosiHybridUG2c/
make rosiHybridUG2c

# echo result
echo "*********************************************************************************************"
echo "Succesfully compiled ROSI models"
echo "*********************************************************************************************"

