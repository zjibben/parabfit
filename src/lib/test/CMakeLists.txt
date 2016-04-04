project(PececilloLibTest)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

# # Add Support for OpenMP.  This needs to be present in the CMakeFile
# # that controls linking of the tests
# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
# endif()

include_directories(${Pececillo_MODULE_DIR})

add_executable(test_hypre_hybrid_type test_hypre_hybrid_type.F90)
target_link_libraries(test_hypre_hybrid_type pececillo)
add_test(hypre_hybrid_type test_hypre_hybrid_type)

#set_property(TEST hypre_hybrid_type PROPERTY LABELS UNIT)
