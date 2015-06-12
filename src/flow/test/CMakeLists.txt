project(PececilloHCTest)

include(CTest)

# Want test modules in the build directory.
unset(CMAKE_Fortran_MODULE_DIRECTORY)

include_directories(${Pececillo_MODULE_DIR})

add_executable(test_HC_AE_solver_type test_HC_AE_solver_type.F90)
target_link_libraries(test_HC_AE_solver_type pececillo-hc)
add_test(hypre_HC_AE_solver test_HC_AE_solver_type)
