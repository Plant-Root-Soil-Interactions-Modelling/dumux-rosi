# .. cmake_module::
#
#    Find the CRootBox root growth library
#
#    You may set the following variables to modify the
#    behaviour of this module:
#
#    :ref:`CROOTBOX_ROOT`
#       Path list to search for CRootBox.
#
#    Sets the following variables:
#
#    :code:`CROOTBOX_FOUND`
#       True if the gstat library was found.
#
#    :code:`CROOTBOX_INCLUDE_DIRS`
#       The path to the CRootBox headers
#
#    :code:`CROOTBOX_LIBRARIES`
#       The CRootBox library path
#
# .. cmake_variable:: CROOTBOX_ROOT
#
#   You may set this variable to have :ref:`FindCRootBox` look
#   for the CRootBox library in the given path before inspecting
#   system paths.
#

# look for header files, only at positions given by the user
find_path(CROOTBOX_INCLUDE_DIR
          NAMES RootSystem.h
          PATHS "${CROOTBOX_ROOT}"
                "${CMAKE_SOURCE_DIR}/../external/crootbox"
                "${CMAKE_SOURCE_DIR}/../external/CRootBox"
                "${CMAKE_SOURCE_DIR}/../crootbox"
                "${CMAKE_SOURCE_DIR}/../CRootBox"
          PATH_SUFFIXES "crootbox" "CRootBox" "SRC" "src"
          NO_DEFAULT_PATH)

# this only looks if it wasn't found before
# search in default paths
find_path(CROOTBOX_INCLUDE_DIR
          NAMES RootSystem.h
          PATH_SUFFIXES "crootbox" "CRootBox" "SRC" "src")

# look for library, at position fiven by the user and some defaults
find_library(CROOTBOX_LIBRARY
             NAMES CRootBox crootbox libCRootBox CRootBox.a crootbox.a libCRootBox.a
             PATHS "${CROOTBOX_ROOT}"
                   "${CMAKE_SOURCE_DIR}/../external/crootbox"
                   "${CMAKE_SOURCE_DIR}/../external/CRootBox"
                   "${CMAKE_SOURCE_DIR}/../crootbox"
                   "${CMAKE_SOURCE_DIR}/../CRootBox"
             PATH_SUFFIXES "" "lib" "lib32" "lib64" "src/.libs" "crootbox" "CRootBox" "SRC" "src"
             NO_DEFAULT_PATH)

# this only looks if it wasn't found before
# search in default paths
find_library(CROOTBOX_LIBRARY
             NAMES CRootBox crootbox libCRootBox CRootBox.a crootbox.a libCRootBox.a
             PATH_SUFFIXES "" "lib" "lib32" "lib64" "src/.libs" "crootbox" "CRootBox" "SRC" "src")

include(CMakePushCheckState)
cmake_push_check_state() # Save variables

# we need if clauses here because variable is set variable-NOTFOUND
if(CROOTBOX_INCLUDE_DIR)
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES} ${CROOTBOX_INCLUDE_DIR})
endif()
if(CROOTBOX_LIBRARY)
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${CROOTBOX_LIBRARY})
endif()

# handle package arguments
# this sets the CRootBox_FOUND and CROOTBOX_FOUND variables
# if CROOTBOX_INCLUDE_DIR and CROOTBOX_LIBRARY contain valid path names
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "CRootBox"
  DEFAULT_MSG
  CROOTBOX_INCLUDE_DIR
  CROOTBOX_LIBRARY
)

#restore old values
cmake_pop_check_state()

# if both headers and library are found, store results
if(CROOTBOX_FOUND)
  set(CROOTBOX_INCLUDE_DIRS ${CROOTBOX_INCLUDE_DIR})
  set(CROOTBOX_LIBRARIES    ${CROOTBOX_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determining location of CRootBox succeeded:\n"
    "Include directory: ${CROOTBOX_INCLUDE_DIRS}\n"
    "Library directory: ${CROOTBOX_LIBRARIES}\n\n")
  set(CROOTBOX_DUNE_COMPILE_FLAGS ${CROOTBOX_INCLUDE_DIRS}
    CACHE STRING "Compile flags used by DUNE when compiling CRootBox programs")
  set(CROOTBOX_DUNE_LIBRARIES ${CROOTBOX_LIBRARIES}
    CACHE STRING "Libraries used by DUNE when linking CRootBox programs")
else()
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determining location of CRootBox failed:\n"
    "Include directory: ${CROOTBOX_INCLUDE_DIRS}\n"
    "Library directory: ${CROOTBOX_LIBRARIES}\n\n")
endif()

# set HAVE_CROOTBOX for config.h
set(HAVE_CROOTBOX ${CROOTBOX_FOUND})

# register all CRootBox related flags
if(CROOTBOX_FOUND)
  dune_register_package_flags(LIBRARIES "${CROOTBOX_LIBRARIES}"
                              INCLUDE_DIRS "${CROOTBOX_INCLUDE_DIRS}")
endif()

mark_as_advanced(CROOTBOX_INCLUDE_DIRS CROOTBOX_LIBRARIES HAVE_CROOTBOX)

# text for feature summary
set_package_properties("CRootBox" PROPERTIES
                       DESCRIPTION "A root growth algorithm (https://plant-root-soil-interactions-modelling.github.io/CRootBox/)"
                       PURPOSE "Simulating growing root systems of different plant species with C++/Python")

function(dumux_add_crootbox_flags _targets)
 foreach(_target ${_targets})
   if(CROOTBOX_INCLUDE_DIRS)
     target_include_directories(${_target} PUBLIC ${CROOTBOX_INCLUDE_DIRS})
   endif()
   if(CROOTBOX_LIBRARIES)
     target_link_libraries(${_target} ${CROOTBOX_LIBRARIES})
   endif()
 endforeach()
endfunction()
