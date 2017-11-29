dune_add_test(SOURCES test_rosiMucilage.cc)
dune_symlink_to_source_files(FILES "grids" "test_rosiMucilage.input")
#ADD_CUSTOM_COMMAND(COMMAND "ln -s $PWD/.. $PWD/src")
