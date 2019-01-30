# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "python" "input")

add_executable(decoupled EXCLUDE_FROM_ALL decoupled.cc)

# add_executable(coupled EXCLUDE_FROM_ALL coupled.cc)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
