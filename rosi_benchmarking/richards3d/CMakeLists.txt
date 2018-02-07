dune_add_test(SOURCES richards3d.cc)
dune_symlink_to_source_files(FILES "grids" "input" "richards3d.input" )

#install sources
install(FILES
richardsproblem3d.hh
richardsparams.hh
richards3d.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux-rosi/rosi_benchmarking/richards3d)
