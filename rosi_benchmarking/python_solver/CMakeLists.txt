# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")


include(DunePybindxiUtils)
dune_pybindxi_add_module(richards_yasp_solver richardsyaspsolver.cc) # solverbase.cc ???
target_link_dune_default_libraries(richards_yasp_solver)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
