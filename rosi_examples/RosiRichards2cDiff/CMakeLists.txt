dune_add_test(SOURCES test_rosiRichards2cDiff.cc)
dune_symlink_to_source_files(FILES "grids" "test_rosiRichards2cDiff.input")
ADD_CUSTOM_COMMAND(COMMAND "ln -s $PWD/.. $PWD/src")
