# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "input" "modelparameter" "python")

add_executable(rootsystem_1pnc EXCLUDE_FROM_ALL rootsystem_1pnc.cc)
target_compile_definitions(rootsystem_1pnc PUBLIC DGF)

add_executable(rootsystem_periodic_1pnc EXCLUDE_FROM_ALL rootsystem_periodic_1pnc.cc)
target_compile_definitions(rootsystem_periodic_1pnc PUBLIC DGF)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
