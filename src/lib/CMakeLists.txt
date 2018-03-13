set(SRC
# General utilities
    f90_assert.F90
    kinds.F90
    string_utilities.F90
    bitfield_type.F90
    permutations.F90
    graph_type.F90
    logging_services.F90
    data_layout_type.F90
    array_utils.F90
    set_type.F90
    #string_type.F90

# Region modules
    region_class.F90
    fill_region_type.F90
    plane_region_type.F90
    sphere_region_type.F90
    halfsphere_region_type.F90
    cylinder_region_type.F90
    ellipsoid_region_type.F90
    sinusoid_region_type.F90
    region_factories.F90

# Mesh modules
    unstr_mesh_type.F90
    unstr_mesh_factory.F90
    cell_topology.F90
    cell_geometry.F90
    unstr_mesh_tools.F90
    facet_hash_type.F90
    facet_table_type.F90
    exodus_c_binding.F90
    exodus_file_type.F90
    mesh_geom_type.F90
    mesh_subset_type.F90
    mixed_cell_subset_constructor.F90

# Function modules
    scalar_func_class.F90
    scalar_func_containers.F90
    const_scalar_func_type.F90
    poly_scalar_func_type.F90
    mpoly_scalar_func_type.F90
    tabular_scalar_func_type.F90
    fptr_scalar_func_type.F90
    piecewise_scalar_func_type.F90
    scalar_func_factories.F90
    paraboloid_type.F90

# Mesh function modules
    bndry_func_class.F90
    bc_factory_type.F90
    bndry_face_func_type.F90
    unstr_mesh_func.F90

# Linear algebra modules
    csr_matrix_type.F90
    csr_precon_class.F90
    csr_precon_ssor_type.F90
    csr_precon_boomer_type.F90
    hypre_c_binding.F90
    hypre_ext.c
    fhypre.F90
    hypre_pcg_type.F90
    hypre_fgmres_type.F90
    hypre_hybrid_type.F90
    upper_packed_matrix.F90

# Mimetic discretization modules
    mfd_disc_type.F90
    mfd_diff_matrix_type.F90
    mfd_diff_precon_type.F90

# DAE integrator modules
    idaesol_type.F90
    nka_type.F90

# Flow Modules
    differential_operators_module.F90
    matl_props_type.F90
    prescribed_velocity_fields.F90
    velocity_to_faces_func.F90
    predictor_module.F90
    fischer_guess_type.F90
    body_force_module.F90
    projection_module.F90

# Vof modules
    consts.F90
    hex_types.F90
    material_geometry_type.F90
    vof_init.F90
    vof_init_ex_circle.F90
    truncate_volume_module.F90
    brent_min_class.F90
    brent_root_class.F90
    bfgs_min_class.F90
    int_norm_module.F90
    lvira_normals.F90
    curvature_hf.F90
    curvature.F90
    curvature_fit_func.F90
    locate_plane_module.F90
    locate_plane_nd_module.F90
    volume_track_module.F90
    volume_track_nd_module.F90
    flux_volume_module.F90
    plane_type.F90
    polygon_type.F90
    pure_polyhedron_type.F90
    polyhedron_type.F90
    multimat_cell_type.F90
    surface_type.F90
    interface_patch_type.F90
    vof_solver_type.F90

# Output modules
    gmvwrite.c
    gmvwrite_fortran_binding.F90
    unstr_mesh_gmv.F90
    )

set(HEADERS f90_assert.fpp)

add_library(pec_core ${SRC} ${HEADERS})
add_library(pececillo::core ALIAS pec_core)

set(LIB_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
set_target_properties(pec_core PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_MOD_DIR})

target_link_libraries(pec_core petaca hypre exodus lapack)

target_include_directories(pec_core PUBLIC
    $<BUILD_INTERFACE:${LIB_MOD_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS pec_core
        EXPORT  pec_core
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(DIRECTORY ${LIB_MOD_DIR} DESTINATION include)
install(FILES ${HEADERS} DESTINATION include)

# if(ENABLE_TESTS)
#   add_subdirectory(test)
# endif()
