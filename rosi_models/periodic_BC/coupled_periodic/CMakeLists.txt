# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "python" "input" "grids" "modelparameter")

add_executable(decoupled_periodic EXCLUDE_FROM_ALL decoupled_periodic.cc)
target_compile_definitions(decoupled_periodic PUBLIC DGF)

add_executable(coupled_seq_periodic EXCLUDE_FROM_ALL coupled_seq_periodic.cc)
target_compile_definitions(coupled_seq_periodic PUBLIC DGF)

add_executable(coupled_seq_rb_periodic EXCLUDE_FROM_ALL coupled_seq_periodic.cc)
target_compile_definitions(coupled_seq_rb_periodic PUBLIC ROOTBOX)

add_executable(coupled_periodic EXCLUDE_FROM_ALL coupled_periodic.cc)
target_compile_definitions(coupled_periodic PUBLIC DGF)

add_executable(coupled_rb_periodic EXCLUDE_FROM_ALL coupled_periodic.cc)
target_compile_definitions(coupled_rb_periodic PUBLIC ROOTBOX)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
