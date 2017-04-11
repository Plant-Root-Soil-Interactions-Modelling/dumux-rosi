dune_add_test(SOURCES test_3d.cc)
dune_symlink_to_source_files(FILES "grids" "test_3d.input")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/simulation)
