project(Pececillo-flow)

include_directories(${Pececillo_SOURCE_DIR}/src/lib)
include_directories(${Pececillo_MODULE_DIR} ${Petaca_MODULE_DIR})

set(SOURCE_FILES
  NS_solver_type.F90
  flow_sim_type.F90)

# Create a library (so we can unit test)
add_library(libpececillo-flow ${SOURCE_FILES})
target_link_libraries(libpececillo-flow pececillo)

# Add Support for OpenMP.  This needs to be present in the CMakeFile
# that controls linking of the final executable.
# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
# endif()


# HC driver
add_executable(pececillo-flow flow_main.F90)
target_link_libraries(pececillo-flow libpececillo-flow)

#add_subdirectory(test)

install(TARGETS pececillo-flow DESTINATION lib)
install(TARGETS pececillo-flow DESTINATION bin)
