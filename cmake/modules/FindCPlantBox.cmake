# .. cmake_module::
#
#    Find the CPlantBox root growth library
#
#    You may set the following variables to modify the
#    behaviour of this module:
#
#    :ref:`CPLANTBOX_ROOT`
#       Path list to search for CPlantBox.
#
#    Sets the following variables:
#
#    :code:`CPLANTBOX_FOUND`
#       True if the gstat library was found.
#
#    :code:`CPLANTBOX_INCLUDE_DIRS`
#       The path to the CPlantBox headers
#
#    :code:`CPLANTBOX_LIBRARIES`
#       The CPlantBox library path
#
# .. cmake_variable:: CPLANTBOX_ROOT
#
#   You may set this variable to have :ref:`FindCPlantBox` look
#   for the CPlantBox library in the given path before inspecting
#   system paths.
#

# look for header files, only at positions given by the user
find_path(CPLANTBOX_INCLUDE_DIR
          NAMES RootSystem.h
          PATHS "${CPLANTBOX_ROOT}"
                "${CMAKE_SOURCE_DIR}/../external/cplantbox"
                "${CMAKE_SOURCE_DIR}/../external/CPlantBox"
                "${CMAKE_SOURCE_DIR}/../cplantbox"
                "${CMAKE_SOURCE_DIR}/../CPlantBox"
          PATH_SUFFIXES "cplantbox" "CPlantBox" "SRC" "src"
          NO_DEFAULT_PATH)

# this only looks if it wasn't found before
# search in default paths
find_path(CPLANTBOX_INCLUDE_DIR
          NAMES RootSystem.h
          PATH_SUFFIXES "cplantbox" "CPlantBox" "SRC" "src")

# look for library, at position fiven by the user and some defaults
find_library(CPLANTBOX_LIBRARY
             NAMES CPlantBox cplantbox libCPlantBox CPlantBox.a cplantbox.a libCPlantBox.a
             PATHS "${CPLANTBOX_ROOT}"
                   "${CMAKE_SOURCE_DIR}/../external/cplantbox"
                   "${CMAKE_SOURCE_DIR}/../external/CPlantBox"
                   "${CMAKE_SOURCE_DIR}/../cplantbox"
                   "${CMAKE_SOURCE_DIR}/../CPlantBox"
             PATH_SUFFIXES "" "lib" "lib32" "lib64" "src/.libs" "cplantbox" "CPlantBox" "SRC" "src"
             NO_DEFAULT_PATH)

# this only looks if it wasn't found before
# search in default paths
find_library(CPLANTBOX_LIBRARY
             NAMES CPlantBox cplantbox libCPlantBox CPlantBox.a cplantbox.a libCPlantBox.a
             PATH_SUFFIXES "" "lib" "lib32" "lib64" "src/.libs" "cplantbox" "CPlantBox" "SRC" "src")

include(CMakePushCheckState)
cmake_push_check_state() # Save variables

# we need if clauses here because variable is set variable-NOTFOUND
if(CPLANTBOX_INCLUDE_DIR)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${CPLANTBOX_INCLUDE_DIR})
endif()
if(CPLANTBOX_LIBRARY)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${CPLANTBOX_LIBRARY})
endif()

# handle package arguments
# this sets the CPlantBox_FOUND and CPLANTBOX_FOUND variables
# if CPLANTBOX_INCLUDE_DIR and CPLANTBOX_LIBRARY contain valid path names
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "CPlantBox"
  DEFAULT_MSG
  CPLANTBOX_INCLUDE_DIR
  CPLANTBOX_LIBRARY
)

#restore old values
cmake_pop_check_state()

# if both headers and library are found, store results
if(CPLANTBOX_FOUND)
  set(CPLANTBOX_INCLUDE_DIRS ${CPLANTBOX_INCLUDE_DIR})
  set(CPLANTBOX_LIBRARIES    ${CPLANTBOX_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of CPlantBox succeeded:\n"
    "Include directory: ${CPLANTBOX_INCLUDE_DIRS}\n"
    "Library directory: ${CPLANTBOX_LIBRARIES}\n\n")
  set(CPLANTBOX_DUNE_COMPILE_FLAGS ${CPLANTBOX_INCLUDE_DIRS}
    CACHE STRING "Compile flags used by DUNE when compiling CPlantBox programs")
  set(CPLANTBOX_DUNE_LIBRARIES ${CPLANTBOX_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking CPlantBox programs")
else()
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining location of CPlantBox failed:\n"
    "Include directory: ${CPLANTBOX_INCLUDE_DIRS}\n"
    "Library directory: ${CPLANTBOX_LIBRARIES}\n\n")
endif()

# set HAVE_CPLANTBOX for config.h
set(HAVE_CPLANTBOX ${CPLANTBOX_FOUND})

# register all CPlantBox related flags
if(CPLANTBOX_FOUND)
  dune_register_package_flags(LIBRARIES "${CPLANTBOX_LIBRARIES}"
                              INCLUDE_DIRS "${CPLANTBOX_INCLUDE_DIRS}")
endif()

mark_as_advanced(CPLANTBOX_INCLUDE_DIRS CPLANTBOX_LIBRARIES HAVE_CPLANTBOX)

# text for feature summary
set_package_properties("CPlantBox" PROPERTIES
                       DESCRIPTION "A plant growth algorithm (https://plant-root-soil-interactions-modelling.github.io/CPlantBox/)"
                       PURPOSE "Simulating growing plants with C++/Python")

function(dumux_add_cplantbox_flags _targets)
 foreach(_target ${_targets})
   if(CPLANTBOX_INCLUDE_DIRS)
     target_include_directories(${_target} PUBLIC ${CPLANTBOX_INCLUDE_DIRS})
   endif()
   if(CPLANTBOX_LIBRARIES)
     target_link_libraries(${_target} ${CPLANTBOX_LIBRARIES})
   endif()
 endforeach()
endfunction()
