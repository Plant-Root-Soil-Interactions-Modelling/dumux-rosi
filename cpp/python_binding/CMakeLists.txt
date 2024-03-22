
# include(DunePybindxiUtils)
# dune_pybindxi_add_module(rosi_richards py_richards.cc) 
# target_link_dune_default_libraries(rosi_richards)
#dune_pybindxi_add_module(rosi_richards_cyl py_richards_cyl.cc)
#dune_pybindxi_add_module(rosi_richardsnc py_richardsnc.cc) 
#dune_pybindxi_add_module(rosi_richardsnc_cyl py_richardsnc_cyl.cc) 
#target_link_dune_default_libraries(rosi_richards_cyl)
#target_link_dune_default_libraries(rosi_richardsnc)
#target_link_dune_default_libraries(rosi_richardsnc_cyl)


#include(CMakePrintHelpers)

#pybind11_add_module(rosi_richards_cyl SHARED py_richards_cyl.cc)
#pybind11_add_module(rosi_richardsnc SHARED py_richardsnc.cc)
#pybind11_add_module(rosi_richardsnc_cyl SHARED py_richardsnc_cyl.cc)

message("libraries *********************************************************************")
message(${dune-uggrid_LIBRARIES})
message(${DUNE_LIBS})

#message(${dune-common_LIBRARIES})
#cmake_print_variables(${dune-spgrid_LIBRARIES})
#message("libraries")

set(PYBIND11_DIR ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/pybind11)
#message(${PYBIND11_DIR})

# Include Pybind11
add_subdirectory(${PYBIND11_DIR} pybind11)
#add_subdirectory(${pybind11_dir})
#include_directories(${PYBIND11_DIR}/include)

#add_subdirectory(pybind11)

pybind11_add_module(rosi_richards SHARED py_richards.cc)
target_link_dune_default_libraries(rosi_richards)  
target_link_libraries(rosi_richards PUBLIC duneuggrid ${DUNE_LIBS}) 

pybind11_add_module(rosi_richardsnc SHARED py_richardsnc.cc)
target_link_dune_default_libraries(rosi_richardsnc)  
target_link_libraries(rosi_richardsnc PUBLIC duneuggrid ${DUNE_LIBS}) 

pybind11_add_module(rosi_richards10c SHARED py_richards10c.cc)
target_link_dune_default_libraries(rosi_richards10c)  
target_link_libraries(rosi_richards10c PUBLIC duneuggrid ${DUNE_LIBS}) 

pybind11_add_module(rosi_richards_cyl SHARED py_richards_cyl.cc)
target_link_dune_default_libraries(rosi_richards_cyl)  
target_link_libraries(rosi_richards_cyl PUBLIC ${DUNE_LIBS}) 

pybind11_add_module(rosi_richardsnc_cyl SHARED py_richardsnc_cyl.cc)
target_link_dune_default_libraries(rosi_richardsnc_cyl)  
target_link_libraries(rosi_richardsnc_cyl PUBLIC ${DUNE_LIBS}) 

pybind11_add_module(rosi_richards10c_cyl SHARED py_richards10c_cyl.cc)
target_link_dune_default_libraries(rosi_richards10c_cyl)  

# target_link_libraries(rosi_richards PUBLIC CPlantBox) # /usr/bin/ld: -lCPlantBox kann nicht gefunden werden: Datei oder Verzeichnis nicht gefunden

#target_link_libraries(rosi_richards PRIVATE ${dune-common_LIBRARIES})
#target_link_libraries(rosi_richards PRIVATE ${dune-geometry_LIBRARIES})
#target_link_libraries(rosi_richards PRIVATE ${dune-grid_LIBRARIES})
#target_link_libraries(rosi_richards PRIVATE ${dune-spgrid_LIBRARIES})
#target_link_libraries(rosi_richards PRIVATE ${dune-uggrid_LIBRARIES})
#target_link_libraries(rosi_richards PUBLIC dumux)
#target_link_libraries(rosi_richards PUBLIC dune)
#target_link_dune_default_libraries(rosi_richards_cyl)
#target_link_dune_default_libraries(rosi_richardsnc)
#target_link_dune_default_libraries(rosi_richardsnc_cyl)
