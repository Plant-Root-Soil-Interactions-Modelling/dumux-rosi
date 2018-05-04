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
    git clone -b releases/2.11 https://git.iws.uni-stuttgart.de/dumux-repositories/dumux.git
else
    echo "Skip cloning dumux because the folder already exists."
    cd dumux
    git checkout releases/2.11
    cd ..
fi

if [ $? -ne 0 ]; then
    echo "*********************************************************************************************"
    echo "Failed to clone the repositories."
    echo "*********************************************************************************************"
    exit $?
fi

echo "*********************************************************************************************"
echo "(2/2) Configure dune modules and dumux. Build the dune libaries. This may take several minutes."
echo "*********************************************************************************************"
# run build
## TODO: UPDATE LINK
./dune-common/bin/dunecontrol --opts=./dumux-rosi/optim.opts --builddir=$(pwd)/build-pubMai2018  all
#
if [ $? -ne 0 ]; then
    echo "*********************************************************************************************"
    echo "Failed to build the dune libaries."
    echo "*********************************************************************************************"
    exit $?
fi

# echo result
echo "*********************************************************************************************"
echo "Succesfully configured and built dune and dumux and dumux-rosi multiscale Mai2018."
echo "*********************************************************************************************"

cd build-pubMai2018/dumux-rosi/rosi_multiscale/1dRhizosphere/
make rosi1d
cd ../2d
make rosi2d
cd ../3d
make rosi3d
cd ../3dMultiscale
make rosiMultiscale
cd ../3dRhizoRoot
make 3dRhizoRoot
cd ..

# echo result
echo "*********************************************************************************************"
echo "Succesfully compiled models from multiscale Mai2018."
echo "*********************************************************************************************"

