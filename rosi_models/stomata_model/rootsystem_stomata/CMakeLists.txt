# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "input" "modelparameter" "python")

add_executable(rootsystem_stomata EXCLUDE_FROM_ALL rootsystem_stomata.cc)
target_compile_definitions(rootsystem_stomata PUBLIC DGF)

add_executable(rootsystemperiodic_stomata EXCLUDE_FROM_ALL rootsystem_periodic_stomata.cc)
target_compile_definitions(rootsystemperiodic_stomata PUBLIC DGF)

add_executable(rootsystem_rb_stomata EXCLUDE_FROM_ALL rootsystem_stomata.cc)
target_compile_definitions(rootsystem_rb_stomata PUBLIC ROOTBOX)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
