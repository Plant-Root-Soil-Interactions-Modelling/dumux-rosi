dune_add_test(SOURCES richards3d.cc)
dune_add_test(SOURCES richards3d_ug.cc)
dune_add_test(SOURCES richards3d_steadystate.cc)
dune_symlink_to_source_files(FILES "grids" "input")

#install sources
install(FILES
richardsproblem3d.hh
richardsproblem3dug.hh
richardsproblem3d_steadystate.hh
richardsparams.hh
richards3d.cc
richards3d_ug.cc
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux-rosi/rosi_benchmarking/richards3d)
