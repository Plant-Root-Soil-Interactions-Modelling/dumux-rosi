# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids")

add_executable(richards3dnc EXCLUDE_FROM_ALL richardsnc.cc)
target_compile_definitions(richards3dnc PUBLIC GRIDTYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>)

add_executable(richards1dnc EXCLUDE_FROM_ALL richardsnc.cc)
target_compile_definitions(richards1dnc PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsUGnc EXCLUDE_FROM_ALL richardsnc.cc)
target_compile_definitions(richardsUGnc PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
