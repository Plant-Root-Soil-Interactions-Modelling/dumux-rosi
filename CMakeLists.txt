cmake_minimum_required(VERSION 3.1)
project(dumux-rosi C CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

# find_package (Python3 COMPONENTS Interpreter Development)
set(PYBIND11_PYTHON_VERSION 3.7.5)
find_package( PythonInterp 3.7.5 REQUIRED )
find_package( PythonLibs 3.7.5 REQUIRED )

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules" ${dune-common_MODULE_PATH} ) 

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()
dune_enable_all_packages()

#add_subdirectory(cmake/modules)
add_subdirectory(cpp)
add_subdirectory(experimental)

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
