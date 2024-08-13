# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")

message(PROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
message(PYBIND11_DIR="${PYBIND11_DIR}")


set(PYBIND11_DIR ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/pybind11)


# Include Pybind11
add_subdirectory(${PYBIND11_DIR} pybind11)

add_executable(richards1d10c_PyBi EXCLUDE_FROM_ALL richards1d10c_PyBi.cc)
target_link_libraries(richards1d10c_PyBi PUBLIC pybind11::embed)
add_executable(richards3d10c_PyBi EXCLUDE_FROM_ALL richards3d10c_PyBi.cc)
target_link_libraries(richards3d10c_PyBi PUBLIC pybind11::embed)
# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE Debug)
