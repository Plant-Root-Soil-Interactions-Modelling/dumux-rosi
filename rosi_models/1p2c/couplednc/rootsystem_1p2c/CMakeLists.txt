# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "input" "modelparameter" "python")

add_executable(rootsystem_1p2c EXCLUDE_FROM_ALL rootsystem_1p2c.cc)
target_compile_definitions(rootsystem_1p2c PUBLIC DGF)

add_executable(rootsystem_rb_1p2c EXCLUDE_FROM_ALL rootsystem_1p2c.cc)
target_compile_definitions(rootsystem_rb_1p2c PUBLIC ROOTBOX)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
