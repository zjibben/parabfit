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
    real(r8) :: sphn(2)
    type(lvira_error) :: norm_error

    integer :: ii,jj,nn, fh
    real(r8) :: ds(2), s(2), tmp(3), tmp2(3), tmp3, smin(2), smax(2)
    real(r8), parameter :: pi = 3.141592653_r8

    call start_timer("lvira normals")

    ! get the initial guess from Youngs' method
    int_norm = interface_normal(vof, mesh, gmesh, .false.)

    m = 1 ! WARN: right now assuming 2 materials
    tmp3 = 0
    !i = 687
!!$omp parallel do private(sphn, norm_error)
    do i = 1,mesh%ncell
      if (vof(m,i) > 1-cutvof .or. vof(m,i) < cutvof) then
        int_norm(:,:,i) = 0
        cycle
      end if
      !if (any(gmesh%cneighbor(:,i)<1)) cycle ! WARN: skipping boundaries.
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
      ! sphn = [pi / 2,               -0.905718691605845_r8]
      ! int_norm(1,m,i) = sin(sphn(1))*cos(sphn(2))
      ! int_norm(2,m,i) = sin(sphn(1))*sin(sphn(2))
      ! int_norm(3,m,i) = cos(sphn(1))

      ! print '(a,2es13.3)', 'sphn:     ',sphn
      ! print '(a,3es13.3)', 'int_norm: ',int_norm(:,m,i)
      ! print '(a,es13.3)',  'err:      ',norm_error%f(sphn)
      ! print *


      ! sphn = [1.54116928196139_r8,  -0.905718691605845_r8]
      ! int_norm(1,m,i) = sin(sphn(1))*cos(sphn(2))
      ! int_norm(2,m,i) = sin(sphn(1))*sin(sphn(2))
      ! int_norm(3,m,i) = cos(sphn(1))

      ! print '(a,2es13.3)', 'sphn:     ',sphn
      ! print '(a,3es13.3)', 'int_norm: ',int_norm(:,m,i)
      ! print '(a,es13.3)',  'err:      ',norm_error%f(sphn)
      ! print *

      ! stop
      ! ! ################


      ! ! DEBUGGING ########
      ! open (newunit=fh, file="err_contour.txt")
      ! nn = 200
      ! smin = [0.0_r8, -pi]
      ! smax = [pi, pi]
      ! ! smin = [0.4_r8*pi, 0.4_r8*pi]
      ! ! smax = [0.6_r8*pi, 0.5_r8*pi]

      ! ds = (smax - smin) / (nn-1)
      ! do ii = 1,nn
      !   do jj = 1,nn
      !     s = smin + [(ii-1)*ds(1), (jj-1)*ds(2)]
      !     write (fh,'(es15.5,a,es15.5,a,es15.5)'), s(1),', ',s(2),', ', norm_error%f(s)
      !   end do
      ! end do
      ! close(fh)
      ! ! ##################
      ! print *
      ! print *, 'grad: ',norm_error%f(sphn)
      ! print *

      call norm_error%find_minimum(sphn, ierr)
      ! sphn = [1.58342508284915_r8, 1.43749657637297_r8] ! converged value
      ! print *, 'conv: ',norm_error%f(sphn)
      ! print *
      !print *, 'conv: ', sphn
      !if (ierr /= 0) call lvira_error_fatal(norm_error)

      ! convert spherical coordinates of normal back to physical coordinates
      !tmp = int_norm(:,m,i)
      int_norm(1,m,i) = sin(sphn(1))*cos(sphn(2))
      int_norm(2,m,i) = sin(sphn(1))*sin(sphn(2))
      int_norm(3,m,i) = cos(sphn(1))


      ! print *, norm2(int_norm(:,m,i) - tmp)
      ! print *, tmp
      ! print '(a,3es13.3)', 'int_norm_lv:   ', int_norm(:,m,i)

      ! call HFCell(tmp3, tmp2, vof(m,:), tmp, mesh, gmesh, i)
      ! print *, tmp2

      ! sphn(1) = acos(tmp2(3))
      ! sphn(2) = atan2(tmp2(2), tmp2(1))

      ! sphn = [1.57079632679490_r8, 1.46797267946713_r8] ! hf value
      ! !print *, sphn
      ! print *, 'hf:   ',norm_error%f(sphn)

      ! sphn = [1.57079632679490_r8, 1.43749657637297_r8] ! z-perp value
      ! print *, 'zprp: ',norm_error%f(sphn)

      tmp3 = max(tmp3,abs(int_norm(3,m,i)))

      ! if (abs(int_norm(3,m,i)) > 2.9e-2) then
      !   print *, i, abs(int_norm(3,m,i))
      ! end if

      ! print *
      ! stop
    end do
!!$omp end parallel do

    int_norm(:,2,:) = - int_norm(:,1,:) ! WARN: right now assuming 1 material

    call stop_timer("lvira normals")

    ! print *, tmp3

    ! call LS_Fatal("lvira debugging")

  contains

    subroutine lvira_error_fatal(norm_error)

      type(lvira_error), intent(in) :: norm_error

      print *, norm_error%numitr

      call LS_Fatal ("could not minimize lvira error")

    end subroutine lvira_error_fatal

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
    !if (ierr /= 0) call lvira_error_fatal(norm_error)

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
    !print *, 'vf: ',this%vof(1)

    this%maxitr = 100; this%line_search_max = 20; this%tol = 1e-6_r8

  end subroutine lvira_error_init

  real(r8) function lvira_error_f(this, x)

    use locate_plane_nd_module
    use plane_type
    use timer_tree_type

    class(lvira_error), intent(in) :: this
    real(r8), intent(in) :: x(:)

    integer :: c, ierr
    real(r8) :: n(3), lvira_vof
    type(plane) :: interface_plane

    call start_timer("lvira error")

    ! get the normal vector from the input angles in spherical coordinates
    n(1) = sin(x(1))*cos(x(2))
    n(2) = sin(x(1))*sin(x(2))
    n(3) = cos(x(1))

    ! locate the plane by matching vofs on the target cell
    interface_plane = locate_plane_nd(this%cell(1), n, &
        this%vof(1)*this%cell_vol(1), this%cell_vol(1))

    lvira_error_f = 0
    do c = 1,this%ncell
      lvira_vof = this%cell(c)%volume_behind_plane(interface_plane, ierr) / this%cell_vol(c)
      !lvira_vof = quick_truncvol(interface_plane, this%cell(c)) / this%cell_vol(c)
      if (ierr /= 0) call LS_fatal ("could not calculate volume behind lvira plane")

      lvira_error_f = lvira_error_f + (lvira_vof - this%vof(c))**2

      ! print '(a,i4,4es18.8)', 'lv: ', c, lvira_vof, this%vof(c), &
      !     lvira_vof - this%vof(c), (lvira_vof - this%vof(c))**2
      ! if ((lvira_vof - this%vof(c))**2 > 1e-7_r8) &
      !     print '(i4,2es13.3)', c, this%vof(c), (lvira_vof - this%vof(c))**2
    end do

    ! print *, 'f: ',lvira_error_f
    ! stop
    call stop_timer("lvira error")

  end function lvira_error_f

  real(r8) function quick_truncvol (interface_plane, cell)

    use consts, only: nfc
    use plane_type
    use hex_types, only: reconstruction_hex
    use truncate_volume_module


    type(plane), intent(in) :: interface_plane
    type(polyhedron), intent(in) :: cell

    integer :: f
    type(truncvol_data) :: trunc_vol(nfc)
    type(reconstruction_hex) :: reccell

    reccell%P = interface_plane
    reccell%node = cell%x

    do f = 1,nfc
      trunc_vol(f) = face_param(reccell, 'full_cell', f)
    end do
    quick_truncvol = truncate_volume(reccell, trunc_vol)

  end function quick_truncvol

end module lvira_normals
