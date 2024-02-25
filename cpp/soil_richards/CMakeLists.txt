# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")

add_executable(richards3d EXCLUDE_FROM_ALL richards.cc)
target_compile_definitions(richards3d PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richards1d EXCLUDE_FROM_ALL richards.cc)
target_compile_definitions(richards1d PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richards1d_cyl EXCLUDE_FROM_ALL richards_cyl.cc)
target_compile_definitions(richards1d_cyl PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richards_alu EXCLUDE_FROM_ALL richards.cc)
target_compile_definitions(richards_alu PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>)

add_executable(richards_ug EXCLUDE_FROM_ALL richards.cc)
target_compile_definitions(richards_ug PUBLIC GRIDTYPE=Dune::UGGrid<3>)

#add_executable(test_richards_benchmark_tpfa main.cc)
# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
 