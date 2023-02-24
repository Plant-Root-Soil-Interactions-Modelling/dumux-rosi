#
# 1. deletes current build
# 2. builds a non-mpi version
# 3. delete again, but keep rosi_richards_cyl
# 4. build mpi-version
# 5. copy rosi_richards_cyl to the right location
#

# run from DUMUX/dumux-rosi/
cd ..

# delete all builds
cd dumux
sudo rm build-cmake -rf
sudo rm CMakeCache.txt
cd ..

cd dumux-rosi
sudo rm build-cmake -rf
sudo rm CMakeCache.txt
cd ..

for MOD in common foamgrid geometry grid istl localfunctions pybindxi spgrid uggrid; do
cd dune-$MOD
sudo rm build-cmake -rf
sudo rm CMakeCache.txt
cd ..
done

# build without MPI
./dune-common/bin/dunecontrol --opts=dumux-rosi/cmake_nompi.opts all

# copy richards_cyl to DUMUX/
cp dumux-rosi/build-cmake/cpp/python_binding/rosi_richards_cyl*.so ./
cp dumux-rosi/build-cmake/cpp/python_binding/rosi_richardsnc_cyl*.so ./

# delete all builds
cd dumux
sudo rm build-cmake -rf
sudo rm CMakeCache.txt
cd ..

cd dumux-rosi
sudo rm build-cmake -rf
sudo rm CMakeCache.txt
cd ..

for MOD in common foamgrid geometry grid istl localfunctions pybindxi spgrid uggrid; do
cd dune-$MOD
sudo rm build-cmake -rf
sudo rm CMakeCache.txt
cd ..
done

# rebuild again with mpi
./dune-common/bin/dunecontrol --opts=dumux-rosi/cmake.opts all

# copy from DUMUX/ the non-mpi richards_cyl (overwriting the MPI version)
cp rosi_richards_cyl*.so dumux-rosi/build-cmake/cpp/python_binding/
cp rosi_richardsnc_cyl*.so dumux-rosi/build-cmake/cpp/python_binding/
