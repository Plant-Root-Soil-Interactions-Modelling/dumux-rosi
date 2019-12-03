# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "python" "input" "grids" "modelparameter" "results")

add_executable(coupled_1pnc_richards EXCLUDE_FROM_ALL coupled_stomata.cc)
target_compile_definitions(coupled_1pnc_richards PUBLIC DGF PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)