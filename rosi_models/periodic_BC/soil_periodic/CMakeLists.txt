# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "benchmarks_1d" "benchmarks_3d" "benchmarks_ug")

add_executable(richards3d_periodic EXCLUDE_FROM_ALL richards_periodic.cc)
target_compile_definitions(richards3d_periodic PUBLIC GRIDTYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>)

add_executable(richards1d_periodic EXCLUDE_FROM_ALL richards_periodic.cc)
target_compile_definitions(richards1d_periodic PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsUG_periodic EXCLUDE_FROM_ALL richards_periodic.cc)
target_compile_definitions(richardsUG_periodic PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
