add_subdirectory(pybind11)

pybind11_add_module(rosi_richards SHARED py_richards.cc)
pybind11_add_module(rosi_richards_cyl SHARED py_richards_cyl.cc)
pybind11_add_module(rosi_richardsnc SHARED py_richardsnc.cc)
pybind11_add_module(rosi_richardsnc_cyl SHARED py_richardsnc_cyl.cc)

set_target_properties(rosi_richards PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
set_target_properties(rosi_richards_cyl PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
set_target_properties(rosi_richardsnc PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)
set_target_properties(rosi_richardsnc_cyl PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)