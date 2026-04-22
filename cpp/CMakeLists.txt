# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "modelparameter")

#add_subdirectory("soil_richards")
add_subdirectory("python_binding")
