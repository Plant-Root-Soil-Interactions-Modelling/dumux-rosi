dune_add_test(SOURCES test_rosiRichards2cbuffernitrate.cc)
dune_symlink_to_source_files(FILES "grids" "test_rosiRichards2cbuffernitrate.input" )
#file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
ADD_CUSTOM_COMMAND(COMMAND "ln -s $PWD/.. $PWD/src")
