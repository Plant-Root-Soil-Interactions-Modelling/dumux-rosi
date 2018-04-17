dune_add_test(SOURCES rootsystem.cc)
dune_symlink_to_source_files(FILES "grids" "input")

#install sources
install(FILES
rootsproblem.hh
rootsparams.hh
rootsystem.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux-rosi/rosi_benchmarking/rootsystem)
