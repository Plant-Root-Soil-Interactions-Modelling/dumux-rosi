dune_add_test(SOURCES test_ccrootsystem.cc)
dune_symlink_to_source_files(FILES "test_ccrootsystem.input" "grids")
#ADD_CUSTOM_COMMAND("ln -s $PWD/.. $PWD/src")
