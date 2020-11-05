dune_add_test(SOURCES test_rosiCitrate.cc)
dune_symlink_to_source_files(FILES "grids" "test_rosiCitrate.input")
#ADD_CUSTOM_COMMAND(COMMAND "ln -s $PWD/.. $PWD/src")
