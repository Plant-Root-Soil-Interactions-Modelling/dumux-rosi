# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "input" "modelparameter" "python")

add_executable(rootsystem_periodic EXCLUDE_FROM_ALL rootsystem_periodic.cc)
target_compile_definitions(rootsystem_periodic PUBLIC DGF)

add_executable(rootsystem_rb_periodic EXCLUDE_FROM_ALL rootsystem_periodic.cc)
target_compile_definitions(rootsystem_rb_periodic PUBLIC ROOTBOX)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
