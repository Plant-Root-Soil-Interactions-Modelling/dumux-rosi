# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "python" "input" "grids")

add_executable(decoupled EXCLUDE_FROM_ALL decoupled.cc)

add_executable(coupled_seq EXCLUDE_FROM_ALL coupled_seq.cc)
target_compile_definitions(coupled_seq PUBLIC DGF)

add_executable(coupled_seq_rb EXCLUDE_FROM_ALL coupled_seq.cc)
target_compile_definitions(coupled_seq_rb PUBLIC ROOTBOX)

add_executable(coupled EXCLUDE_FROM_ALL coupled.cc)
target_compile_definitions(coupled PUBLIC DGF)

add_executable(coupled_rb EXCLUDE_FROM_ALL coupled.cc)
target_compile_definitions(coupled_rb PUBLIC ROOTBOX)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
