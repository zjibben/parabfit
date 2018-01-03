set(SRC
  vof_io.F90
  vof_init_test.F90
  normal_vector_test.F90
  interface_point_test.F90
  polygon_type_test.F90
  polyhedron_type_test.F90
  locate_plane_nd_test.F90
  surface_type_test.F90
  multimat_cell_type_test.F90
  locate_plane_test.F90
  volume_track_test.F90
  volume_track_nd_test.F90
  vof_solver_type_test.F90
  bfgs_min_class_test.F90
  curvature_tests.F90
  parabfit_main.F90)

add_library(pec_parabfit ${SRC})
add_library(pececillo::parabfit ALIAS pec_parabfit)

set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
set_target_properties(pec_parabfit PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR})
target_link_libraries(pec_parabfit pececillo::core petaca)

target_include_directories(pec_parabfit PUBLIC
    $<BUILD_INTERFACE:${LIB_MOD_DIR}>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS pec_parabfit
        EXPORT  pec_parabfit
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(DIRECTORY ${LIB_MOD_DIR} DESTINATION include)

# Paraboloid fitting program
add_executable(pececillo-parabfit parabfit_main.F90)
target_link_libraries(pececillo-parabfit pececillo::parabfit)
install(TARGETS pececillo-parabfit DESTINATION bin)
