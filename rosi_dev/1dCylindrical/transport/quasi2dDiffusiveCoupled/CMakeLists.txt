dune_add_test(SOURCES test_quasi2dDiffusive2_coupled.cc)
dune_symlink_to_source_files(FILES "grids" "test_quasi2dDiffusive2_coupled.input")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/simulation)
