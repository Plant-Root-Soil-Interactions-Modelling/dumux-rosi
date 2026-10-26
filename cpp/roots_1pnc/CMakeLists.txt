# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "grids" "input" "python")

add_executable(rootsystem_stomata2 EXCLUDE_FROM_ALL rootsystem_stomata.cc)
target_compile_definitions(rootsystem_stomata2 PUBLIC DGF)
target_include_directories(rootsystem_stomata2 PRIVATE
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/structural
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/functional
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/tinyxml2
)

add_executable(rootsystem_periodic_stomata EXCLUDE_FROM_ALL rootsystem_stomata.cc)
target_compile_definitions(rootsystem_periodic_stomata PUBLIC DGF)
target_compile_definitions(rootsystem_periodic_stomata PUBLIC ROOTS_PERIODIC)
target_include_directories(rootsystem_periodic_stomata PRIVATE
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/structural
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/functional
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/tinyxml2
)

add_executable(rootsystem_1p2c EXCLUDE_FROM_ALL rootsystem_1p2c.cc)
target_compile_definitions(rootsystem_1p2c PUBLIC DGF)
target_include_directories(rootsystem_1p2c PRIVATE
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/structural
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/functional
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/tinyxml2
)

add_executable(rootsystem_1p2c_rb EXCLUDE_FROM_ALL rootsystem_1p2c.cc)
target_compile_definitions(rootsystem_1p2c_rb PUBLIC ROOTBOX)
target_include_directories(rootsystem_1p2c_rb PRIVATE
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/structural
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/functional
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/tinyxml2
)

add_executable(rootsystem_periodic_1p2c EXCLUDE_FROM_ALL rootsystem_1p2c.cc)
target_compile_definitions(rootsystem_periodic_1p2c PUBLIC DGF)
target_compile_definitions(rootsystem_periodic_1p2c PUBLIC ROOTS_PERIODIC)
target_include_directories(rootsystem_periodic_1p2c PRIVATE
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/structural
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/functional
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/tinyxml2
)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
