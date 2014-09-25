# This is the core Pececillo library.
add_subdirectory(lib)

# This is the basic heat conduction solver.
add_subdirectory(hc)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION include)
