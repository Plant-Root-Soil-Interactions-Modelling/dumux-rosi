add_subdirectory("richards2cbuffer")
add_subdirectory("rootmodel1p2c")
add_subdirectory("richardsmodelnc")
add_subdirectory("compositional") # overwrite to introduce buffer power
add_subdirectory("richardsnc") # overwrite to change io field from moleFraction to massFraction

add_subdirectory("richardsCylindrical1d") # overwrite to change io field from moleFraction to massFraction
