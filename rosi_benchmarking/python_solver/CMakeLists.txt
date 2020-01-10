#
# in order to use it install module dune-pybindxi and add "dune-pybindxi (>= 0)" to dune.module, 
# the uncomment cmake commands 
#

# create a link to the grid file and the input file in the build directory
# dune_symlink_to_source_files(FILES "grids" "python" "input")

# build pyhthon binding 
# include(DunePybindxiUtils)
# dune_pybindxi_add_module(richards_sp_solver richardsspsolver.cc) # solverbase.cc ???
# target_link_dune_default_libraries(richards_sp_solver)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
# set(CMAKE_BUILD_TYPE RelWithDebInfo)
