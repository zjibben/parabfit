module fit_normals

  use kinds, only: r8
  use unstr_mesh_type
  use mesh_geom_type
  use logging_services
  implicit none

  public :: interface_normals_fit

contains

  subroutine interface_normals_fit(int_norm, vof, mesh, gmesh)

    use int_norm_module
    use lvira_normals

    real(r8), allocatable, intent(out) :: int_norm(:,:,:)
    real(r8), intent(in) :: vof(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh

    integer, parameter :: iter_max = 5
    integer :: iter

    ! calculate the initial guess via Youngs' method
    !int_norm = interface_normal(vof, mesh, gmesh, .false.)
    call interface_normals_lvira(int_norm, vof, mesh, gmesh)

    ! iterate
    do iter = 1,iter_max
      call interface_normals_fit_iteration(int_norm, vof, mesh, gmesh)
    end do

  end subroutine interface_normals_fit


  subroutine interface_normals_fit_iteration(int_norm, vof, mesh, gmesh)

    use consts, only: cutvof
    use hex_types, only: hex_f, hex_e
    use interface_patch_type
    use surface_type
    use multimat_cell_type

    real(r8), allocatable, intent(inout) :: int_norm(:,:,:)
    real(r8), intent(in) :: vof(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh

    type(surface) :: intrec
    type(multimat_cell) :: cell
    integer :: i, m, ierr

    m = 1 ! WARN: right now assuming 2 materials

    ! get the interface reconstructions
    do i = 1,mesh%ncell
      call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, gmesh%outnorm(:,:,i), &
          mesh%volume(i))
      if (ierr /= 0) call LS_fatal ('cell_outward_volflux failed: could not initialize cell')

      call cell%partition (vof(:,i), int_norm(:,:,i))
      call intrec%append (cell%interface_polygon(1), i)
    end do

    do i = 1,mesh%ncell
      if (vof(m,i) > 1-cutvof .or. vof(m,i) < cutvof) then
        int_norm(:,:,i) = 0
        cycle
      end if

      ! WARN: right now assuming 2 materials
      int_norm(:,1,i) = normal_from_patch(intrec%local_patch(i, gmesh, vof(1,:)), &
          0.0_r8, int_norm(:,1,i))
    end do
    int_norm(:,2,:) = -int_norm(:,1,:) ! WARN: right now assuming 2 materials

  end subroutine interface_normals_fit_iteration

end module fit_normals
