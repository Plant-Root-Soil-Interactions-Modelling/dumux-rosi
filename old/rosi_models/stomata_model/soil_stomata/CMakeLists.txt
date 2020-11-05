# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")

add_executable(richards3d_stomata EXCLUDE_FROM_ALL richards_stomata.cc)
target_compile_definitions(richards3d_stomata PUBLIC GRIDTYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>)

add_executable(richards1d_stomata EXCLUDE_FROM_ALL richards_stomata.cc)
target_compile_definitions(richards1d_stomata PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsUG_stomata EXCLUDE_FROM_ALL richards_stomata.cc)
target_compile_definitions(richardsUG_stomata PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
