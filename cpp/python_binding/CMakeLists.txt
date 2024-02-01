add_subdirectory(pybind11)

pybind11_add_module(rosi_richards SHARED py_richards.cc)
pybind11_add_module(rosi_richards_cyl SHARED py_richards_cyl.cc)
pybind11_add_module(rosi_richardsnc SHARED py_richardsnc.cc)
pybind11_add_module(rosi_richardsnc_cyl SHARED py_richardsnc_cyl.cc)

