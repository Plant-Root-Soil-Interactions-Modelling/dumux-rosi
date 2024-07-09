# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "python" "input")




set(PYBIND11_DIR ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/pybind11)


# Include Pybind11
add_subdirectory(${PYBIND11_DIR} pybind11)

add_executable(richards1d10c_PyBi richards1d10c_PyBi.cc)
target_link_libraries(richards1d10c_PyBi PUBLIC pybind11::embed)
# target_compile_definitions(richards1d10c_PyBi PUBLIC GRIDTYPE=Dune::FoamGrid<1,1>)
# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
