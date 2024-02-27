# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "modelparameter")

#add_subdirectory("roots_1p")
#add_subdirectory("roots_1pnc")
#add_subdirectory("soil_richards")
add_subdirectory("soil_richardsnc")
#add_subdirectory("soil_tracer")
#add_subdirectory("coupled_1p_richards")
#add_subdirectory("coupled_1pnc_richards")
#add_subdirectory("coupled_1pnc_richardsnc")
add_subdirectory("python_binding")
