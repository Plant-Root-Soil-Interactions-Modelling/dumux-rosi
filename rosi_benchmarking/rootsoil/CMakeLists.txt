dune_add_test(SOURCES richards1d.cc)
dune_symlink_to_source_files(FILES "grids" "input")

#install sources
install(FILES
richardsproblem1d.hh
richardsparams.hh
rootsproblem.hh
rootsparams.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux-rosi/rosi_benchmarking/rootsoil)
