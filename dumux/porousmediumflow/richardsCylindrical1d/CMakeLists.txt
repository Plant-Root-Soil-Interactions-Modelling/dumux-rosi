add_subdirectory(richards)

install(FILES
localresidual.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rosi/richardsCylindrical1d)
