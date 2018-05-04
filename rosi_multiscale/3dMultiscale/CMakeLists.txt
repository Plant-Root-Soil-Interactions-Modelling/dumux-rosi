dune_add_test(SOURCES rosiMultiscale.cc)
dune_symlink_to_source_files(FILES "grids"
"singleRootDryingLoamKb40_10.input" 
"singleRootDryingLoamKb20_10.input"
"singleRootDryingLoamKb80_10.input"
"singleRootDryingLoamKb40_05.input"
"singleRootLoamKb40_10theta0.1.input"
"singleRootLoamKb40_10theta0.2.input" 
"singleRootLoamKb20_10theta0.1.input"
"singleRootLoamKb20_10theta0.2.input"
"singleRootLoamKb80_10theta0.1.input"
"singleRootLoamKb80_10theta0.2.input"
"multiRootDryingLoamKb40_05_40.input"
"multiRootDryingLoamKb05_05_10.input"
"rootSystemDryingLoamKb40_40.input"
"run_rosi3dMultiscale.sh")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/singleRoot)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/rootSystem)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/multiRoot)
