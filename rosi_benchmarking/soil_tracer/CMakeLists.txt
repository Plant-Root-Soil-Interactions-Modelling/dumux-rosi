dune_symlink_to_source_files(FILES "params.input")

# explicit tracer tests (mass fractions)
dune_add_test(NAME soil_tracer_explicit_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=false USEMOLES=false
)
dune_add_test(NAME soil_tracer_explicit_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=false USEMOLES=false
)

# explicit tracer tests (mole fractions, should yield same result)
dune_add_test(NAME soil_tracer_explicit_tpfa_mol
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=false USEMOLES=true
              CMAKE_GUARD HAVE_UMFPACK
)

dune_add_test(NAME soil_tracer_explicit_box_mol
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=false USEMOLES=true
)

# implicit tracer tests
dune_add_test(NAME soil_tracer_implicit_tpfa
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestTpfa IMPLICIT=true USEMOLES=false
)

dune_add_test(NAME soil_tracer_implicit_box
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=TracerTestBox IMPLICIT=true USEMOLES=false
)
