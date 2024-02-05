# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")

add_executable(richardsnc3d EXCLUDE_FROM_ALL richards1p2c.cc)
target_compile_definitions(richardsnc3d PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richardsnc1d EXCLUDE_FROM_ALL richards1p2c.cc)
target_compile_definitions(richardsnc1d PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsnc1d_cyl EXCLUDE_FROM_ALL richards1p2c_cyl.cc)
target_compile_definitions(richardsnc1d_cyl PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsnc_alu  richards1p2c.cc)
target_compile_definitions(richardsnc_alu PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>)

add_executable(richardsnc_ug richards1p2c.cc)
target_compile_definitions(richardsnc_ug PUBLIC GRIDTYPE=Dune::UGGrid<3>)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
