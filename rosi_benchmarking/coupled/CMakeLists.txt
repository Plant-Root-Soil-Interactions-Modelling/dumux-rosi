# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "python" "input" "grids")

add_executable(decoupled EXCLUDE_FROM_ALL decoupled.cc)

add_executable(coupled1 EXCLUDE_FROM_ALL coupled1.cc)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
