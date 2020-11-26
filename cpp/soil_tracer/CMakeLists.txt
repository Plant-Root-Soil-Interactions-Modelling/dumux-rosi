dune_symlink_to_source_files(FILES "grids" "python" "input")

add_executable(test_tracer EXCLUDE_FROM_ALL main.cc)
target_compile_definitions(test_tracer PUBLIC TYPETAG=TracerTestTpfa IMPLICIT=false)

