# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")



set(PYBIND11_DIR ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/pybind11)


# Include Pybind11
add_subdirectory(${PYBIND11_DIR} pybind11)
add_executable(richards10c1d_cyl_pybi EXCLUDE_FROM_ALL richards1d10c_PyBi.cc)
target_link_libraries(richards10c1d_cyl_pybi PUBLIC pybind11::embed)

add_executable(richards10c3d_pybi EXCLUDE_FROM_ALL richards3d10c_PyBi.cc)
target_link_libraries(richards10c3d_pybi PUBLIC pybind11::embed)

add_executable(soil_b5_10c EXCLUDE_FROM_ALL soil_b5_10c.cc)
target_link_libraries(soil_b5_10c PUBLIC pybind11::embed)

add_executable(soil_b5_10cILU EXCLUDE_FROM_ALL soil_b5_10cILU.cc)
target_link_libraries(soil_b5_10cILU PUBLIC pybind11::embed)
add_executable(soil_b5_10cILURes EXCLUDE_FROM_ALL soil_b5_10cILURes.cc)
target_link_libraries(soil_b5_10cILURes PUBLIC pybind11::embed)
add_executable(soil_b5_10cSSORC EXCLUDE_FROM_ALL soil_b5_10cSSORC.cc)
target_link_libraries(soil_b5_10cSSORC PUBLIC pybind11::embed)
add_executable(soil_b5_10cSeq EXCLUDE_FROM_ALL soil_b5_10cSeq.cc)
target_link_libraries(soil_b5_10cSeq PUBLIC pybind11::embed)

add_executable(soil_b5 EXCLUDE_FROM_ALL soil_b5.cc)
target_link_libraries(soil_b5 PUBLIC pybind11::embed)
# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
