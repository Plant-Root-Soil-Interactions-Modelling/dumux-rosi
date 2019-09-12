# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "python" "input" "grids" "modelparameter")

add_executable(decoupled_stomata EXCLUDE_FROM_ALL decoupled_stomata.cc)
target_compile_definitions(decoupled_stomata PUBLIC DGF)

add_executable(coupled_seq_stomata EXCLUDE_FROM_ALL coupled_seq_stomata.cc)
target_compile_definitions(coupled_seq_stomata PUBLIC DGF)

add_executable(coupled_seq_rb_stomata EXCLUDE_FROM_ALL coupled_seq_stomata.cc)
target_compile_definitions(coupled_seq_rb_stomata PUBLIC ROOTBOX)

add_executable(coupled_stomata EXCLUDE_FROM_ALL coupled_stomata.cc)
target_compile_definitions(coupled_stomata PUBLIC DGF)

add_executable(coupled_rb_stomata EXCLUDE_FROM_ALL coupled_stomata.cc)
target_compile_definitions(coupled_rb_stomata PUBLIC ROOTBOX)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
