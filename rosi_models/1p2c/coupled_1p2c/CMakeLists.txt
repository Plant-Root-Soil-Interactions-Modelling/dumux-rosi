# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "python" "input" "grids" "modelparameter")

add_executable(decoupled_1p2c EXCLUDE_FROM_ALL decoupled_1p2c.cc)
target_compile_definitions(decoupled_1p2c PUBLIC DGF)

add_executable(coupled_seq_1p2c EXCLUDE_FROM_ALL coupled_seq_1p2c.cc)
target_compile_definitions(coupled_seq_1p2c PUBLIC DGF)

add_executable(coupled_seq_rb_1p2c EXCLUDE_FROM_ALL coupled_seq_1p2c.cc)
target_compile_definitions(coupled_seq_rb_1p2c PUBLIC ROOTBOX)

add_executable(coupled_1p2c EXCLUDE_FROM_ALL coupled_1p2c.cc)
target_compile_definitions(coupled_1p2c PUBLIC DGF)

add_executable(coupled_rb_1p2c EXCLUDE_FROM_ALL coupled_1p2c.cc)
target_compile_definitions(coupled_rb_1p2c PUBLIC ROOTBOX)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
