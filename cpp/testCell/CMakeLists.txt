
add_executable(testCell main.cc)
target_compile_definitions(testCell PUBLIC GRIDTYPE=Dune::SPGrid<double,3>)
set(CMAKE_BUILD_TYPE Debug)
