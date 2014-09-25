project(PececilloLibTest)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Pececillo_MODULE_DIR})

add_executable(test_hypre_hybrid_type test_hypre_hybrid_type.F90)
target_link_libraries(test_hypre_hybrid_type pececillo)
add_test(hypre_hybrid_type test_hypre_hybrid_type)

#set_property(TEST hypre_hybrid_type PROPERTY LABELS UNIT)
