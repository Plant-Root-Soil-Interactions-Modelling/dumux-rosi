dune_add_test(SOURCES test_rosiRichards2cExud.cc)
dune_symlink_to_source_files(FILES "grids" "test_rosiRichards2cExud.input" "results.pvsm")
ADD_CUSTOM_COMMAND(COMMAND "ln -s $PWD/.. $PWD/src")
