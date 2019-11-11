# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")

include(DunePybindxiUtils)
dune_pybindxi_add_module(richards_yasp_solver richardsyaspsolver.cc solverbase.cc)
target_link_dune_default_libraries(richards_yasp_solver)

add_executable(richards3d EXCLUDE_FROM_ALL richards.cc)
target_compile_definitions(richards3d PUBLIC GRIDTYPE=Dune::YaspGrid<3,Dune::EquidistantOffsetCoordinates<double,3>>)

add_executable(richards1d EXCLUDE_FROM_ALL richards.cc)
target_compile_definitions(richards1d PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsUG EXCLUDE_FROM_ALL richards.cc)
target_compile_definitions(richardsUG PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
