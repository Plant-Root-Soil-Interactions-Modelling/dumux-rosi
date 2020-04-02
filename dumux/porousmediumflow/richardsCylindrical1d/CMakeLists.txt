add_subdirectory(richards)

install(FILES
model.hh
localresidual.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/rosi/porousmediumflow/richardsCylindrical1d)
