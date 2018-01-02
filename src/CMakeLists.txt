add_subdirectory(petaca) # the petaca library
add_subdirectory(lib)    # core pececillo library
add_subdirectory(unit)   # unit tests

#install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION include)
