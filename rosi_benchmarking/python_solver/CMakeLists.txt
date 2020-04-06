#
# in order to use it install module dune-pybindxi and add "dune-pybindxi (>= 0)" to dune.module, 
# the uncomment cmake commands 
#

# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")

# build pyhthon binding 
find_package (Python3 COMPONENTS Interpreter Development)
set(PYBIND11_PYTHON_VERSION 3.7.5)
set(DUNE_PYBINDXI_PYTHON_VERSION 3.7.5)

include(DunePybindxiUtils)
dune_pybindxi_add_module(dumux_rosi py_richards.cc) 
dune_pybindxi_add_module(dumux_rosi_cyl py_richards_cyl.cc) # we need to build seperately, otherwise types are redefined and linker complains 
target_link_dune_default_libraries(dumux_rosi)
target_link_dune_default_libraries(dumux_rosi_cyl)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
