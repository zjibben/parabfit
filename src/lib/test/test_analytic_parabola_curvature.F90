program test_analytic_parabola_curvature

  use kinds, only: r8
  use logging_services
  implicit none

  integer :: status, N, p
  real(r8) :: error, R

  R = 0.35_r8

  p = 8
  N = 64 * 2**(p-1) + 1
  call analytic_test(error, R, 1.0_r8 / N, status, verbose=.true.)
  print *
  call simulate_test(error, R, 1.0_r8 / N, status) !, verbose=.true.)
  print *

  do p = 1,6
    N = 64 * 2**(p-1) + 1
    !call analytic_test(error, R, 1.0_r8 / N, status)
    call simulate_test(error, R, 1.0_r8 / N, status)
    print '(2(a,es13.3),a)', '[',1.0_r8 / N,', ', error,'],'
  end do

  call exit(status)

contains

  subroutine analytic_test(error, R, dx, status, verbose)

    real(r8), intent(out) :: error
    real(r8), intent(in) :: R, dx
    integer, intent(out) :: status
    logical, intent(in), optional :: verbose

    real(r8) :: y0, a, b, c, curvature, volfrac(-1:1), normal(2,-1:1), rho(-1:1), x(2,-1:1), &
        rnd(2,-1:1)
    logical :: verboseh

    status = 0
    if (present(verbose)) then
      verboseh = verbose
    else
      verboseh = .false.
    end if

    ! set y0
    !y0 = R - dx/2 ! the circle lands halfway through the center cell
    y0 = floor(R / dx) * dx - dx

    if (sqrt(R**2 - (3*dx/2)**2) < y0) then
      status = 1
      return
    end if

    ! calculate volume fractions
    volfrac(0) = (Yf(dx/2, R) - Yf(-dx/2, R) - y0*dx) / dx**2
    volfrac(1) = (Yf(3*dx/2, R) - Yf(dx/2, R) - y0*dx) / dx**2
    volfrac(-1) = volfrac(1) ! symmetry

    ! calculate normal vectors
    normal(1,1) = -sqrt(R**2 - (3*dx/2)**2)
    normal(1,1) = normal(1,1) + sqrt(R**2 - (dx/2)**2)
    normal(2,1) = dx
    normal(:,1) = normal(:,1) / norm2(normal(:,1))

    normal(:,0) = [0.0_r8, 1.0_r8] ! symmetry
    normal(:,-1) = [-normal(1,1), normal(2,1)] ! symmetry

    ! ! insert error
    ! call random_number(rnd)
    ! normal = normal + rnd * 1e7 * dx**2
    ! call normalize(normal(:,-1))
    ! call normalize(normal(:,0))
    ! call normalize(normal(:,1))

    ! calculate plane offsets
    rho(-1) = normal(2,-1) * (dx*volfrac(-1) + y0) - normal(1,-1) * dx
    rho(0)  = normal(2,0) * (dx*volfrac(0) + y0)
    rho(1)  = normal(2,1) * (dx*volfrac(1) + y0) + normal(1,1) * dx

    ! calculate plane centroids
    x(:,-1) = [-dx, (rho(-1) + normal(1,-1)*dx) / normal(2,-1)]
    x(:,0)  = [0.0_r8, rho(0) / normal(2,0)]
    x(:,1)  = [dx, (rho(1) - normal(1,1)*dx) / normal(2,1)]

    ! calculate parabola
    a = x(2,0)
    b = (x(2,1) - x(2,-1)) / (2*dx)
    c = (x(2,1) - 2*x(2,0) + x(2,-1)) / (2*dx**2)

    ! calculate curvature
    curvature = 2*c / sqrt(1 + b**2)

    error = abs(curvature + 1 / R) * R

    if (verboseh) then
      print '(a,3es15.5)', 'vof:       ', volfrac

      print '(a,3es15.5)', 'normal:    ', normal(:,-1)
      print '(a,3es15.5)', 'normal:    ', normal(:,0)
      print '(a,3es15.5)', 'normal:    ', normal(:,1)

      print '(a,3es15.5)', 'rho:       ', rho

      print '(a,2es20.10)', 'x:         ', x(:,-1)
      print '(a,2es20.10)', 'x:         ', x(:,0)
      print '(a,2es20.10)', 'x:         ', x(:,1)

      print '(a,3es15.5)', 'parabola:  ', a, b, c

      print '(a,es15.5)', 'curvature: ', curvature
      print '(a,es15.5)', 'error:     ', error
    end if

  end subroutine analytic_test

  pure real(r8) function Yf(x, R)
    real(r8), intent(in) :: x, R
    Yf = (x * sqrt(R**2 - x**2) + R**2 * atan(x/sqrt(R**2 - x**2))) / 2
  end function Yf

  real(r8) function rand()
    call random_number(rand)
  end function rand

  subroutine normalize(n)
    real(r8), intent(inout) :: n(:)
    n = n / norm2(n)
  end subroutine normalize

  subroutine simulate_test(error, R, dx, status, verbose)

    use vof_init_ex_circle
    use unstr_mesh_factory
    use unstr_mesh_type
    use mesh_geom_type
    use interface_patch_type
    use multimat_cell_type
    use surface_type
    use lvira_normals
    use int_norm_module
    use hex_types, only: hex_f, hex_e

    use polygon_type
    use array_utils, only: isZero

    real(r8), intent(out) :: error
    real(r8), intent(in) :: R, dx
    integer, intent(out) :: status
    logical, intent(in), optional :: verbose

    type(unstr_mesh) :: mesh
    type(mesh_geom) :: gmesh
    type(multimat_cell) :: cell
    type(surface) :: intrec
    integer :: mesh_size, i, j, c, ierr, cells(-1:1)
    real(r8), allocatable :: vof(:,:), int_norm(:,:,:)
    real(r8) :: curvature
    logical :: verboseh

    type(polygon) :: tmp

    status = 0
    if (present(verbose)) then
      verboseh = verbose
    else
      verboseh = .false.
    end if

    ! initialize mesh
    mesh_size = nint(1 / dx)
    mesh = new_unstr_mesh ([-21 * dx / 2, dx/2, -3 * dx / 2], &
        [21 * dx / 2, 0.5_r8, 3 * dx / 2], [21,(mesh_size - 1)/2,3])
    ! mesh = new_unstr_mesh ([-0.5_r8, -0.5_r8, -3 * dx / 2], &
    !     [0.5_r8, 0.5_r8, 3 * dx / 2], [mesh_size,mesh_size,3])
    call gmesh%init (mesh)

    allocate(vof(2,mesh%ncell), int_norm(3,2,mesh%ncell))

    ! find the center cell which contains the interface
    do i = 1,mesh%ncell
      if (.not.any(gmesh%cneighbor(:,i)<1) .and. gmesh%xc(1,i) == 0 .and. &
          maxval(mesh%x(2,mesh%cnode(:,i))) > R .and. minval(mesh%x(2,mesh%cnode(:,i))) < R) then
        j = i
        exit
      else if (i == mesh%ncell) then
        call LS_fatal ("couldn't find center cell")
      end if
    end do

    ! find the two neighbors
    cells = [gmesh%cneighbor(3,j), j, gmesh%cneighbor(4,j)]

    ! initialize vof
    call vof_init_circle(mesh, R, vof)

    ! calculate normal vectors
    int_norm = interface_normal(vof, mesh, gmesh, .false.) ! initial guess
    call interface_normal_lvira(int_norm(:,1,j), j, vof(1,:), mesh, gmesh)
    ! print *, 'h0', int_norm(:,1,j)
    ! print *, 'h1', interface_normal_lvira(j, vof(1,:), mesh, gmesh)
    do c = 1,gmesh%caneighbor(j)%n_elements
      i = gmesh%caneighbor(j)%elements(c)

      call interface_normal_lvira(int_norm(:,1,i), i, vof(1,:), mesh, gmesh)
    end do
    int_norm(:,2,:) = -int_norm(:,1,:)

    !print *, 'h2', interface_normal_lvira(j, vof(1,:), mesh, gmesh)

    ! calculate interface reconstructions
    do c = 1,gmesh%caneighbor(j)%n_elements + 1
      if (c > gmesh%caneighbor(j)%n_elements) then
        i = j
      else
        i = gmesh%caneighbor(j)%elements(c)
      end if

      if (.not.isMixedCell(vof(1,i))) cycle

      call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, gmesh%outnorm(:,:,i), &
          mesh%volume(i))
      if (ierr /= 0) call LS_fatal ('could not initialize cell')

      call cell%partition (vof(:,i), int_norm(:,:,i))
      call intrec%append (cell%interface_polygons(1), i)

      ! if (i == j) then
      !   print *, vof(1,i)
      !   print *, cell%mat_poly(1)%volume() / mesh%volume(i)
      !   print *, isZero(cell%volume())
      !   tmp = cell%interface_polygon(1)
      !   call tmp%print_data()
      !   print *
      ! end if
    end do

    !print *, 'h3', interface_normal_lvira(j, vof(1,:), mesh, gmesh)

    ! calculate curvature
    curvature = abs(curvature_from_patch (intrec%local_patch(j,gmesh, vof(1,:)), &
        0.0_r8, int_norm(:,1,j), vof(1,:), mesh, gmesh, j))
    error = abs(-curvature + 1/R) * R

    if (verboseh) then
      !print *, j, mesh_size

      print *, j, cells
      print '(a,3es18.8)', 'vof:       ', vof(1,cells)

      print '(a,3es15.5)', 'normal:    ', int_norm(:2,2,cells(-1))
      print '(a,3es15.5)', 'normal:    ', int_norm(:2,2,cells(0))
      print '(a,3es15.5)', 'normal:    ', int_norm(:2,2,cells(1))

      !print '(a,3es15.5)', 'rho:       ', rho

      ! do c = -1, 1
      !   do i = 1,size(intrec%element)
      !     if (intrec%cell_id(i) == cells(c)) then
      !       print '(a,3es15.5)', 'x:         ', intrec%element(i)%centroid()
      !       exit
      !     end if
      !   end do
      ! end do

      curvature = abs(curvature_from_patch (intrec%local_patch(j,gmesh, vof(1,:)), &
        0.0_r8, int_norm(:,1,j), vof(1,:), mesh, gmesh, j, verboseh))

      print '(a,es15.5)', 'curvature: ', curvature
      print '(a,es15.5)', 'error:     ', error
    end if

  end subroutine simulate_test

  logical function isMixedCell(vof)
    use consts, only: cutvof
    real(r8), intent(in) :: vof
    isMixedCell = vof < 1-cutvof .and. vof > cutvof
  end function isMixedCell

end program test_analytic_parabola_curvature
