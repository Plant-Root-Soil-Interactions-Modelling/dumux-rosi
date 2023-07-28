# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")

add_executable(richards3c3d richards1p3c.cc)
target_compile_definitions(richards3c3d PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richards3c3dR richards1p3cR.cc)
target_compile_definitions(richards3c3dR PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richardsnc3d richards1pnc.cc)
target_compile_definitions(richardsnc3d PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richards3c1dR_cyl richards1p3cR_cyl.cc)
target_compile_definitions(richards3c1dR_cyl PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richards5c1dR_cyl richards1p5cR_cyl.cc)
target_compile_definitions(richards5c1dR_cyl PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsnc1d EXCLUDE_FROM_ALL richards1pnc.cc)
target_compile_definitions(richardsnc1d PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsnc1d_cyl EXCLUDE_FROM_ALL richards1pnc_cyl.cc)
target_compile_definitions(richardsnc1d_cyl PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richards2c3d richards1p2c.cc)
target_compile_definitions(richards2c3d PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richards2c1d EXCLUDE_FROM_ALL richards1p2c.cc)
target_compile_definitions(richards2c1d PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richards2c1d_cyl richards1p2c_cyl.cc)
target_compile_definitions(richards2c1d_cyl PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richardsnc_alu EXCLUDE_FROM_ALL richards1p2c.cc)
target_compile_definitions(richardsnc_alu PUBLIC GRIDTYPE=Dune::ALUGrid<3,3,Dune::simplex,Dune::conforming>)

add_executable(richardsnc_ug EXCLUDE_FROM_ALL richards1p2c.cc)
target_compile_definitions(richardsnc_ug PUBLIC GRIDTYPE=Dune::UGGrid<3>)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
