#
# in order to use it install module dune-pybindxi and add "dune-pybindxi (>= 0)" to dune.module, 
# the uncomment cmake commands 
#



set(HEADERS
rosi_richards.hh
)

# we need to build seperately, otherwise types are redefined and the linker complains 
# target_link_dune_default_libraries(rosi_richards)
#

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/)
