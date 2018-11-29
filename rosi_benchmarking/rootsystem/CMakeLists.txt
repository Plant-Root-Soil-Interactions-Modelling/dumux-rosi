# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "input" "modelparameter")

dune_add_test(NAME rootsystem
              SOURCES rootsystem.cc
              CMAKE_GUARD HAVE_CROOTBOX)

dune_add_test(NAME main
              SOURCES main.cc
              CMAKE_GUARD HAVE_CROOTBOX
              COMPILE_DEFINITIONS TYPETAG=TubesTestBox
              CMAKE_GUARD dune-foamgrid_FOUND
              COMMAND ./convergencetest.py
              CMD_ARGS test_1p_network1d3d_box params.input -Problem.Name test_1p_network1d3d_box)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)




