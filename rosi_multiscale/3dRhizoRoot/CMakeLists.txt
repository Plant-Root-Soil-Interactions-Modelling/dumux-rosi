dune_add_test(SOURCES 3dRhizoRoot.cc)
dune_symlink_to_source_files(FILES "grids"
"multiRootDryingLoamKb05_05_10.input")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/singleRoot)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/rootSystem)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/multiRoot)
