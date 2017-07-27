module lvira_normals

  use kinds,  only: r8
  implicit none
  private

contains

  subroutine interface_normal_lvira(int_norm, vof, mesh, gmesh)

    use int_norm_module
    use timer_tree_type

    real(r8), intent(out) :: int_norm(:,:,:)
    real(r8), intent(in)  :: vof(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh

    call start_timer("lvira normals")

    ! get the initial guess from Youngs' method
    int_norm = interface_normal(vof, mesh, gmesh, .false.)



    call stop_timer("lvira normals")

  end subroutine interface_normal_lvira

end module lvira_normals
