dune_add_test(SOURCES test_ccrootsystemgrowth.cc)
dune_symlink_to_source_files(FILES "test_ccrootsystemgrowth.input" "grids" "$PWD/..")
#ADD_CUSTOM_COMMAND(COMMAND "ln -s $PWD/.. $PWD/src")
