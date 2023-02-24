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
rm build-cmake -rf
rm CMakeCache.txt
cd ..

cd dumux-rosi
rm build-cmake -rf
rm CMakeCache.txt
cd ..

for MOD in common foamgrid geometry grid istl localfunctions pybindxi spgrid uggrid alugrid; do
cd dune-$MOD
rm build-cmake -rf
rm CMakeCache.txt
cd ..
done

# build without MPI
./dune-common/bin/dunecontrol --opts=dumux-rosi/cmake_nompi.opts all
# ./dune-common/bin/dunecontrol --opts=dumux-rosi/cmake_nompi_cluster.opts all


