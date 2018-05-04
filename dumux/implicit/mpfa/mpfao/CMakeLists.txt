add_subdirectory("2d")

#install headers
install(FILES
mpfaointeractionvolumecontainer.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/mpfa/mpfao)
