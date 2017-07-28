module lvira_normals

  use kinds, only: r8
  use bfgs_min_class
  use polyhedron_type
  use logging_services
  implicit none
  private

  public :: interface_normal_lvira

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

  subroutine interface_normal_lvira(int_norm, vof, mesh, gmesh)

    use int_norm_module
    use timer_tree_type
    use unstr_mesh_type
    use mesh_geom_type
    use consts, only: cutvof

    real(r8), allocatable, intent(out) :: int_norm(:,:,:)
    real(r8), intent(in) :: vof(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh

    integer :: i, m, ierr
    real(r8) :: sphn(2)
    type(lvira_error) :: norm_error

    integer :: ii,jj,nn
    real(r8) :: ds(2), s(2)
    real(r8), parameter :: pi = 3.141592653_r8

    call start_timer("lvira normals")

    ! get the initial guess from Youngs' method
    int_norm = interface_normal(vof, mesh, gmesh, .false.)

    m = 1 ! WARN: right now assuming 1 material
    do i = 1,mesh%ncell
      if (vof(m,i) > 1-cutvof .or. vof(m,i) < cutvof) cycle
      ! print *, i, mesh%ncell

      ! print '(a,es20.10)', 'vof:           ',vof(m,i)
      ! print *, vof(m,i) > 1-cutvof, vof(m,i) < cutvof

      ! print '(a,3es13.3)', 'int_norm_grad: ',int_norm(:,m,i)

      ! convert normal to spherical coordinates
      sphn(1) = acos(int_norm(3,m,i))
      sphn(2) = atan2(int_norm(2,m,i), int_norm(1,m,i))

      !print '(a,3es13.3)', 'external s:    ',sphn

      ! find the normal that minimizes the lvira error function
      call norm_error%init(i, vof(m,:), mesh, gmesh)

      ! ! DEBUGGING ########
      ! open (99, file="err_contour.txt")
      ! nn = 200
      ! ds(1) = pi / (nn-1)
      ! ds(2) = 2*pi / (nn-1)
      ! do ii = 1,nn
      !   do jj = 1,nn
      !     s = [(ii-1)*ds(1), (jj-1)*ds(2)]
      !     write (99,'(es15.5,a,es15.5,a,es15.5)'), s(1),', ',s(2),', ', norm_error%f(s)
      !   end do
      ! end do
      ! close(99)
      ! ! ##################


      call norm_error%find_minimum(sphn, ierr)
      !if (ierr /= 0) call lvira_error_fatal(norm_error)

      ! convert spherical coordinates of normal back to physical coordinates
      int_norm(1,m,i) = sin(sphn(1))*cos(sphn(2))
      int_norm(2,m,i) = sin(sphn(1))*sin(sphn(2))
      int_norm(3,m,i) = cos(sphn(2))
    end do

    int_norm(:,2,:) = - int_norm(:,1,:) ! WARN: right now assuming 1 material

    call stop_timer("lvira normals")

  contains

    subroutine lvira_error_fatal(norm_error)

      type(lvira_error), intent(in) :: norm_error

      print *, norm_error%numitr

      call LS_Fatal ("could not minimize lvira error")

    end subroutine lvira_error_fatal

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

    allocate(this%cell(this%ncell), this%vof(this%ncell), this%cell_vol(this%ncell))

    ! initialize target cell
    call this%cell(1)%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, mesh%volume(i), &
        gmesh%outnorm(:,:,i))
    this%vof(1) = vof(i)
    this%cell_vol(1) = mesh%volume(i)

    ! initialize node neighbors
    do c = 2,this%ncell
      cid = gmesh%caneighbor(i)%elements(c-1)

      call this%cell(c)%init (ierr, mesh%x(:,mesh%cnode(:,cid)), hex_f, hex_e, mesh%volume(cid), &
          gmesh%outnorm(:,:,cid))
      this%vof(c) = vof(cid)
      this%cell_vol(c) = mesh%volume(cid)
    end do

  end subroutine lvira_error_init

  real(r8) function lvira_error_f(this, x)

    use locate_plane_nd_module
    use plane_type

    class(lvira_error), intent(in) :: this
    real(r8), intent(in) :: x(:)

    integer :: c, ierr
    real(r8) :: n(3), lvira_vof
    type(plane) :: interface_plane

    ! get the normal vector from the input angles in spherical coordinates
    n(1) = sin(x(1))*cos(x(2))
    n(2) = sin(x(1))*sin(x(2))
    n(3) = cos(x(1))

    ! locate the plane by matching vofs on the target cell
    interface_plane = locate_plane_nd(this%cell(1), n, &
        this%vof(1)*this%cell_vol(1), this%cell_vol(1))

    lvira_error_f = 0
    do c = 1,this%ncell
      lvira_vof = this%cell(c)%volume_behind_plane(interface_plane, ierr)
      if (ierr /= 0) call LS_fatal ("could not calculate volume behind lvira plane")

      lvira_error_f = lvira_error_f + (lvira_vof - this%vof(c))**2
    end do

    ! print *, 'f: ',lvira_error_f
    ! stop

  end function lvira_error_f

end module lvira_normals
