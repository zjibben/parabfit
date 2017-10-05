module lvira_normals

  use kinds, only: r8
  use bfgs_min_class
  use polyhedron_type
  use logging_services
  implicit none
  private

  public :: interface_normals_lvira, interface_normal_lvira

  type, extends(bfgs_min) :: lvira_error
    private
    integer :: ncell
    type(polyhedron), allocatable :: cell(:)
    real(r8), allocatable :: vof(:), cell_vol(:)
  contains
    procedure :: init => lvira_error_init
    procedure :: f => lvira_error_f
  end type lvira_error

contains

  subroutine interface_normals_lvira(int_norm, vof, mesh, gmesh)

    use int_norm_module
    use timer_tree_type
    use unstr_mesh_type
    use mesh_geom_type
    use consts, only: cutvof
    use curvature_hf, only: HFCell ! DEBUGGING

    real(r8), allocatable, intent(out) :: int_norm(:,:,:)
    real(r8), intent(in) :: vof(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh

    integer :: i, m, ierr

    call start_timer("lvira normals")

    ! get the initial guess from Youngs' method
    int_norm = interface_normal(vof, mesh, gmesh, .false.)
    m = 1 ! WARN: right now assuming 2 materials

    ! Dynamically schedule since only a few cells are mixed
    ! and they take the vast majority of the runtime.

    !$omp parallel do schedule(dynamic,100)
    do i = 1,mesh%ncell
      call interface_normal_lvira(int_norm(:,m,i), i, vof(m,:), mesh, gmesh)
    end do
    !$omp end parallel do

    int_norm(:,2,:) = - int_norm(:,1,:) ! WARN: right now assuming 1 material

    call stop_timer("lvira normals")

  end subroutine interface_normals_lvira

  subroutine interface_normal_lvira(int_norm, i, vof, mesh, gmesh)

    use unstr_mesh_type
    use mesh_geom_type
    use consts, only: cutvof

    real(r8), intent(inout) :: int_norm(:)
    integer, intent(in) :: i
    real(r8), intent(in) :: vof(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh

    integer :: ierr
    real(r8) :: sphn(2)
    type(lvira_error) :: norm_error

    if (vof(i) > 1-cutvof .or. vof(i) < cutvof) then
      int_norm = 0
      return
    end if

    ! convert normal to spherical coordinates
    sphn(1) = acos(int_norm(3))
    sphn(2) = atan2(int_norm(2), int_norm(1))

    ! find the normal that minimizes the lvira error function
    call norm_error%init(i, vof, mesh, gmesh)

    call norm_error%find_minimum(sphn, ierr)
    !if (ierr /= 0) call LS_Fatal("lvira error")

    ! convert spherical coordinates of normal back to physical coordinates
    int_norm(1) = sin(sphn(1))*cos(sphn(2))
    int_norm(2) = sin(sphn(1))*sin(sphn(2))
    int_norm(3) = cos(sphn(1))

  end subroutine interface_normal_lvira

  subroutine lvira_error_init(this, i, vof, mesh, gmesh)

    use hex_types, only: hex_f, hex_e
    use unstr_mesh_type
    use mesh_geom_type

    class(lvira_error), intent(out) :: this
    integer, intent(in) :: i
    real(r8), intent(in) :: vof(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh

    integer :: c, cid, ierr

    this%ncell = gmesh%caneighbor(i)%n_elements + 1

    !print *, 'nc: ',this%ncell
    allocate(this%cell(this%ncell), this%vof(this%ncell), this%cell_vol(this%ncell))

    ! initialize target cell
    call this%cell(1)%init (ierr, i, mesh, gmesh)
    ! call this%cell(1)%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, gmesh%outnorm(:,:,i), &
    !     mesh%volume(i))
    this%vof(1) = vof(i)
    this%cell_vol(1) = mesh%volume(i)

    ! initialize node neighbors
    do c = 2,this%ncell
      cid = gmesh%caneighbor(i)%elements(c-1)

      call this%cell(c)%init (ierr, cid, mesh, gmesh)
      ! call this%cell(c)%init (ierr, mesh%x(:,mesh%cnode(:,cid)), hex_f, hex_e, &
      !     gmesh%outnorm(:,:,cid), mesh%volume(cid))
      this%vof(c) = vof(cid)
      this%cell_vol(c) = mesh%volume(cid)
    end do
    !print *, 'vf: ',this%vof(1)

    this%maxitr = 100; this%line_search_max = 20; this%tol = 1e-6_r8

  end subroutine lvira_error_init

  real(r8) function lvira_error_f(this, x)

    use locate_plane_nd_module
    use plane_type
    use timer_tree_type

    class(lvira_error), intent(inout) :: this
    real(r8), intent(in) :: x(:)

    integer :: c, ierr
    real(r8) :: n(3), lvira_vof
    type(plane) :: interface_plane

    !call start_timer("lvira error")

    ! get the normal vector from the input angles in spherical coordinates
    n(1) = sin(x(1))*cos(x(2))
    n(2) = sin(x(1))*sin(x(2))
    n(3) = cos(x(1))

    ! locate the plane by matching vofs on the target cell
    interface_plane = locate_plane_nd(this%cell(1), n, &
        this%vof(1)*this%cell_vol(1), this%cell_vol(1))

    lvira_error_f = 0
    do c = 2,this%ncell
      lvira_vof = this%cell(c)%volume_behind_plane(interface_plane, ierr) / this%cell_vol(c)
      if (ierr /= 0) call LS_fatal ("could not calculate volume behind lvira plane")

      lvira_error_f = lvira_error_f + (lvira_vof - this%vof(c))**2

      ! print '(a,i4,4es18.8)', 'lv: ', c, lvira_vof, this%vof(c), &
      !     lvira_vof - this%vof(c), (lvira_vof - this%vof(c))**2
      ! if ((lvira_vof - this%vof(c))**2 > 1e-7_r8) &
      !     print '(i4,2es13.3)', c, this%vof(c), (lvira_vof - this%vof(c))**2
    end do

    ! print *, 'f: ',lvira_error_f
    ! stop
    !call stop_timer("lvira error")

  end function lvira_error_f

end module lvira_normals
