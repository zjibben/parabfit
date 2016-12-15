project(Pececillo-unit-tests)

include_directories(${Pececillo_SOURCE_DIR}/src/lib)
include_directories(${Pececillo_MODULE_DIR} ${Petaca_MODULE_DIR})

set(SOURCE_FILES
  analytic_surface_type_test.F90
  polygon_type_test.F90
  polyhedron_type_test.F90
  locate_plane_nd_test.F90
  surface_type_test.F90
  multimat_cell_type_test.F90
  locate_plane_test.F90
  volume_track_test.F90
  volume_track_nd_test.F90
  vof_solver_type_test.F90
  string_type_test.F90
  unit_test.F90)

# Create a library (so we can unit test)
add_library(libpececillo-unit ${SOURCE_FILES})
target_link_libraries(libpececillo-unit pececillo)

# Add Support for OpenMP.  This needs to be present in the CMakeFile
# that controls linking of the final executable.
# find_package(OpenMP)
# if (OPENMP_FOUND)
#     set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#     set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#     set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_C_FLAGS}")
# endif()


# unit test driver
add_executable(pececillo-unit unit_test.F90)
target_link_libraries(pececillo-unit libpececillo-unit)

#add_subdirectory(test)

install(TARGETS pececillo-unit DESTINATION lib)
install(TARGETS pececillo-unit DESTINATION bin)
