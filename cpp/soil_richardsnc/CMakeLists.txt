# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")

add_executable(richards3c3d EXCLUDE_FROM_ALL richards1p3c.cc)
target_compile_definitions(richards3c3d PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richards3c3dR EXCLUDE_FROM_ALL richards1p3cR.cc)
target_compile_definitions(richards3c3dR PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richardsnc3d EXCLUDE_FROM_ALL richards1pnc.cc)
target_compile_definitions(richardsnc3d PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richards3c1dR_cyl EXCLUDE_FROM_ALL richards1p3cR_cyl.cc)
target_compile_definitions(richards3c1dR_cyl PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richards5c1dR_cyl EXCLUDE_FROM_ALL richards1p5cR_cyl.cc)
target_compile_definitions(richards5c1dR_cyl PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richards10c1d_cyl richards1p10c_cyl.cc)
target_compile_definitions(richards10c1d_cyl PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

add_executable(richards10c_PyBi richards1p10c_PyBi.cc)
target_compile_definitions(richards10c_PyBi PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richards3d10c_PyBi_ads richards3d1p10c_PyBi_ads.cc)
target_compile_definitions(richards3d10c_PyBi_ads PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richards10c richards1p10c.cc)
target_compile_definitions(richards10c PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(testCell EXCLUDE_FROM_ALL testCell.cc)
target_compile_definitions(testCell PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)

add_executable(richards5c1dR_cyl_k EXCLUDE_FROM_ALL richards1p5cR_cyl_k.cc)
target_compile_definitions(richards5c1dR_cyl_k PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)

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
