dune_add_test(SOURCES rosi2d.cc)
dune_symlink_to_source_files(FILES "grids" "multiRootDryingLoamKb40_05.input" "multiRootDryingLoamKb05_05.input" "run2d.sh")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
