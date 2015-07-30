project(Pececillo-HC)

include_directories(${Pececillo_SOURCE_DIR}/src/lib)
include_directories(${Pececillo_MODULE_DIR} ${Petaca_MODULE_DIR})

set(SOURCE_FILES
    HC_model_type.F90
    HC_precon_type.F90
    HC_AE_solver_type.F90
    HC_norm_type.F90
    HC_idaesol_model_type.F90
    HC_solver_type.F90
    HC_sim_type.F90)

# Create a library (so we can unit test)
add_library(pececillo-hc ${SOURCE_FILES})
target_link_libraries(pececillo-hc pececillo)

# Add Support for OpenMP.  This needs to be present in the CMakeFile
# that controls linking of the final executable.
# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
# endif()


# HC driver
add_executable(pececillo-HC HC_main.F90)
target_link_libraries(pececillo-HC pececillo-hc)

add_subdirectory(test)

install(TARGETS pececillo-hc DESTINATION lib)
install(TARGETS pececillo-HC DESTINATION bin)
