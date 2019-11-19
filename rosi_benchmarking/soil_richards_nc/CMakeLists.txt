# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")

add_executable(richardsnc3d EXCLUDE_FROM_ALL richardsnc.cc)
target_compile_definitions(richardsnc3d PUBLIC GRIDTYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>)

add_executable(richardsnc1d EXCLUDE_FROM_ALL richardsnc.cc)
target_compile_definitions(richardsnc1d PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsncUG EXCLUDE_FROM_ALL richardsnc.cc)
target_compile_definitions(richardsncUG PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)