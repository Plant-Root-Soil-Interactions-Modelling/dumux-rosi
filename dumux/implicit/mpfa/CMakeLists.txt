add_subdirectory("mpfal")
add_subdirectory("mpfao")

#install headers
install(FILES
mpfaassembler.hh
mpfaelementvolumevariables.hh
mpfafvelementgeometry.hh
mpfahelper.hh
mpfaimplicitdarcyfluxvariables.hh
mpfaimplicitlocaljacobian.hh
mpfaimplicitlocalresidual.hh
mpfaimplicitmodel.hh
mpfaimplicitporousmediaproblem.hh
mpfaimplicitproblem.hh
mpfaimplicitvelocityoutput.hh
mpfalocalresidual.hh
mpfaproperties.hh
mpfapropertydefaults.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/implicit/mpfa)
