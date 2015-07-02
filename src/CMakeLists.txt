# This is the core Pececillo library.
add_subdirectory(lib)

# This is the basic heat conduction solver.
add_subdirectory(hc)

# This is the basic flow solver.
add_subdirectory(flow)

# A unit test executable
add_subdirectory(unit)

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION include)
