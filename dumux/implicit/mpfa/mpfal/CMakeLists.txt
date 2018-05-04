add_subdirectory("2d")

#install headers
install(FILES
mpfalinteractionvolumecontainer.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/mpfa/mpfal)
