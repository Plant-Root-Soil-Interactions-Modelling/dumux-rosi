# create a link to the grid file and the input file in the build directory
dune_symlink_to_source_files(FILES "python" "input" "grids" "modelparameter")

add_executable(coupled EXCLUDE_FROM_ALL coupled.cc)
target_compile_definitions(coupled PUBLIC DGF)

add_executable(coupled_periodic EXCLUDE_FROM_ALL coupled_periodic.cc)
target_compile_definitions(coupled_periodic PUBLIC DGF)

add_executable(coupled_ug EXCLUDE_FROM_ALL coupled.cc)
target_compile_definitions(coupled_ug PUBLIC DGF GRIDTYPE=Dune::UGGrid<3>)

add_executable(coupled_rb EXCLUDE_FROM_ALL coupled.cc)
target_compile_definitions(coupled_rb PUBLIC ROOTBOX)

# add_executable(coupled_schroeder EXCLUDE_FROM_ALL coupled_schroeder.cc ../../dumux/external/brent/brent.cpp)
# target_compile_definitions(coupled_schroeder PUBLIC DGF)

# optionally set cmake build type (Release / Debug / RelWithDebInfo)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
