dune_add_test(SOURCES rosi1d.cc)
dune_symlink_to_source_files(FILES 
"singleRootLoamKb40_theta0.1.input" 
"singleRootLoamKb40_theta0.2.input"
"singleRootLoamKb20_theta0.1.input" 
"singleRootLoamKb20_theta0.2.input"
"singleRootLoamKb80_theta0.1.input" 
"singleRootLoamKb80_theta0.2.input"
"singleRootDryingLoamKb20.input"
"singleRootDryingLoamKb40.input"
"singleRootDryingLoamKb80.input")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
