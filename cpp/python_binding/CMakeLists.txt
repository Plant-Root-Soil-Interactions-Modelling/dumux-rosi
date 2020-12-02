#
# in order to use it install module dune-pybindxi and add "dune-pybindxi (>= 0)" to dune.module, 
# the uncomment cmake commands 
#


# build pyhthon binding 
# find_package (Python3 COMPONENTS Interpreter Development)
# set(PYBIND11_PYTHON_VERSION 3.7.5)
# set(DUNE_PYBINDXI_PYTHON_VERSION 3.7.5)

# we need to build seperately, otherwise types are redefined and the linker complains 
include(DunePybindxiUtils)
dune_pybindxi_add_module(rosi_richards py_richards.cc) 
dune_pybindxi_add_module(rosi_richards_cyl py_richards_cyl.cc)
dune_pybindxi_add_module(rosi_richardsnc_cyl py_richardsnc_cyl.cc) 
target_link_dune_default_libraries(rosi_richards)
target_link_dune_default_libraries(rosi_richards_cyl)
target_link_dune_default_libraries(rosi_richardsnc_cyl)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
