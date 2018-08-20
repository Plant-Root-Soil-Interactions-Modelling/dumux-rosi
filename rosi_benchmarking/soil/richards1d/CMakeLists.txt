dune_add_test(SOURCES richards1d.cc)
dune_symlink_to_source_files(FILES "grids" "input" "richards1d.input" )

#install sources
install(FILES
richardsproblem1d.hh
richardsparams.hh
richards1d.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux-rosi/rosi_benchmarking/richards1d)
