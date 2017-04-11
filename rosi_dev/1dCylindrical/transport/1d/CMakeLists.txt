dune_add_test(SOURCES test_1dRhizosphere.cc)
dune_symlink_to_source_files(FILES "grids" "test_1dRhizosphere.input")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/simulation)
