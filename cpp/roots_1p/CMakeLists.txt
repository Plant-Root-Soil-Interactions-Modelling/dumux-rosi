# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "input" "python")

add_executable(rootsystem EXCLUDE_FROM_ALL rootsystem.cc)
target_compile_definitions(rootsystem PUBLIC DGF)

target_include_directories(rootsystem PRIVATE
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/structural
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/functional
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/tinyxml2
)

add_executable(rootsystem_periodic EXCLUDE_FROM_ALL rootsystem_periodic.cc)
target_compile_definitions(rootsystem_periodic PUBLIC DGF)
target_include_directories(rootsystem_periodic PRIVATE
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/structural
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/functional
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/tinyxml2
)

add_executable(rootsystem_rb EXCLUDE_FROM_ALL rootsystem.cc)
target_compile_definitions(rootsystem_rb PUBLIC ROOTBOX)
target_include_directories(rootsystem_rb PRIVATE
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/structural
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/functional
    ${PROJECT_SOURCE_DIR}/../CPlantBox/src/external/tinyxml2
)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
