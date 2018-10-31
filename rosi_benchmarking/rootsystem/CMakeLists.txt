# create a link to the crootbox parameters folder
if(HAVE_CROOTBOX)
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CROOTBOX_INCLUDE_DIR}/modelparameter ${CMAKE_CURRENT_BINARY_DIR}/params)
endif()

# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "input" "modelparameter")

# add the test executable (if CRootBox was found)
dune_add_test(NAME rootsystem
              SOURCES rootsystem.cc
              CMAKE_GUARD HAVE_CROOTBOX)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE Debug)




