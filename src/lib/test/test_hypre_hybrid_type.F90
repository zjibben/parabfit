program test_hypre_hybrid_type

#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  use kinds
  use csr_matrix_type
  use hypre_hybrid_type
  use parameter_list_type
  implicit none

  real(r8), parameter :: PI = 3.1415926535897931_r8
  real(r8), parameter :: TWOPI = 6.2831853071795862_r8
  integer :: status = 0

  real(r8) :: a ! set by the tests
  integer :: nx, ny, nz ! set by the tests

  call cg_test_1
  call cg_test_2
  call gmres_test_1
  call gmres_test_2
  call exit (status)

contains

  subroutine cg_test_1

    type(csr_matrix), target :: matrix
    type(parameter_list), pointer :: params
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: x(:), b(:), u(:)
    integer :: nrow, num_itr, num_dscg_itr, num_pcg_itr, stat
    real(r8) :: maxerr, l2err, norm

    write(*,'(/,a)') 'Running CG_TEST_1'

    a = 1.0_r8
    nx = 9; ny = 9; nz = 9

    allocate(params)
    call params%set ('krylov-method', 'cg')
    call params%set ('rel-tol', 1.0e-8_r8)
    call params%set ('conv-rate-tol', 0.8_r8)
    call params%set ('max-ds-iter', 50)
    call params%set ('max-amg-iter', 10)
    call params%set ('cg-use-two-norm', .true.)
    call params%set ('logging-level', 1)
    call params%set ('print-level', 0)

    call create_matrix (matrix)
    call solver%init (matrix, params)
    call solver%setup

    nrow = matrix%nrow
    allocate(u(nrow), b(nrow), x(nrow))

    call eigenvector (1, 1, 1, u)
    b = eigenvalue(1, 1, 1) * u
    x = 0.0_r8; x(1) = 1.0_r8 ! seed the initial residual with all eigen-components
    call solver%solve (b, x, stat)
    call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
    if (stat /= 0) then
      write(*,'(a)') 'error: solver failed to converge'
      status = 1
    end if

    deallocate(params)

    maxerr = maxval(abs(u-x))
    l2err = sqrt(sum((u-x)**2))
    write(*,'(3(a,i0),a,es9.2)') 'hypre: num_itr = ', num_itr, ' (', num_dscg_itr, &
        ', ', num_pcg_itr, '), ||r||_2/||b||_2 =', norm
    write(*,'(2(a,es9.2))') 'test:  ||x-x_exact||_max =', maxerr, ', ||x-x_exact||_2 =', l2err

    if (num_dscg_itr /= 22) then
      write(*,'(a,i0)') 'error: expected 22 diagonally-scaled CG iterations; got ', num_dscg_itr
      status = 1
    end if

    if (num_pcg_itr /= 0) then
      write(*,'(a,i0)') 'error: expected 0 AMG preconditioned CG iterations; got ', num_pcg_itr
      status = 1
    end if

    if (maxerr > 4.0e-9) then
      write(*,'(a,es9.2)') 'error: expected maxerr < 4.0e-9; got', maxerr
      status = 1
    end if

  end subroutine


  subroutine cg_test_2

    type(csr_matrix), target :: matrix
    type(parameter_list), pointer :: params
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: x(:), b(:), u(:)
    integer :: nrow, num_itr, num_dscg_itr, num_pcg_itr, stat
    real(r8) :: maxerr, l2err, norm

    write(*,'(/,a)') 'Running CG_TEST_2'

    a = 1.0e-4_r8
    nx = 17; ny = 17; nz = 17

    allocate(params)
    call params%set ('krylov-method', 'cg')
    call params%set ('rel-tol', 1.0e-8_r8)
    call params%set ('conv-rate-tol', 0.6_r8)
    call params%set ('max-ds-iter', 50)
    call params%set ('max-amg-iter', 10)
    call params%set ('cg-use-two-norm', .true.)
    call params%set ('logging-level', 1)
    call params%set ('print-level', 0)

    call create_matrix (matrix)
    call solver%init (matrix, params)
    call solver%setup

    nrow = matrix%nrow
    allocate(u(nrow), b(nrow), x(nrow))

    call eigenvector (1, 1, 1, u)
    b = eigenvalue(1, 1, 1) * u
    x = 0.0_r8; x(1) = 1.0_r8 ! seed the initial residual with all eigen-components
    call solver%solve (b, x, stat)
    call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
    if (stat /= 0) then
      write(*,'(a)') 'error: solver failed to converge'
      status = 1
    end if

    deallocate(params)

    maxerr = maxval(abs(u-x))
    l2err = sqrt(sum((u-x)**2))
    write(*,'(3(a,i0),a,es9.2)') 'hypre: num_itr = ', num_itr, ' (', num_dscg_itr, &
        ', ', num_pcg_itr, '), ||r||_2/||b||_2 =', norm
    write(*,'(2(a,es9.2))') 'test:  ||x-x_exact||_max =', maxerr, ', ||x-x_exact||_2 =', l2err

    if (num_dscg_itr /= 18) then
      write(*,'(a,i0)') 'error: expected 18 diagonally-scaled CG iterations; got ', num_dscg_itr
      status = 1
    end if

    if (num_pcg_itr /= 4) then
      write(*,'(a,i0)') 'error: expected 4 AMG preconditioned CG iterations; got ', num_pcg_itr
      status = 1
    end if

    if (maxerr > 2.0e-9) then
      write(*,'(a,es9.2)') 'error: expected maxerr < 2.0e-9; got', maxerr
      status = 1
    end if

  end subroutine


  subroutine gmres_test_1

    type(csr_matrix), target :: matrix
    type(parameter_list), pointer :: params
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: x(:), b(:), u(:)
    integer :: nrow, num_itr, num_dscg_itr, num_pcg_itr, stat
    real(r8) :: maxerr, l2err, norm

    write(*,'(/,a)') 'Running GMRES_TEST_1'

    a = 1.0_r8
    nx = 9; ny = 9; nz = 9

    allocate(params)
    call params%set ('krylov-method', 'gmres')
    call params%set ('rel-tol', 1.0e-6_r8)
    call params%set ('conv-rate-tol', 0.8_r8)
    call params%set ('max-ds-iter', 50)
    call params%set ('max-amg-iter', 10)
    call params%set ('gmres-krylov-dim', 5)
    call params%set ('amg-smoothing-sweeps', 1)
    call params%set ('logging-level', 1)
    call params%set ('print-level', 1)

    call create_matrix (matrix)
    call solver%init (matrix, params)
    call solver%setup

    nrow = matrix%nrow
    allocate(u(nrow), b(nrow), x(nrow))

    call eigenvector (1, 1, 1, u)
    b = eigenvalue(1, 1, 1) * u
    x = 0.0_r8; x(1) = 1.0_r8 ! seed the initial residual with all eigen-components
    call solver%solve (b, x, stat)
    call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
    if (stat /= 0) then
      write(*,'(a)') 'error: solver failed to converge'
      status = 1
    end if

    deallocate(params)

    maxerr = maxval(abs(u-x))
    l2err = sqrt(sum((u-x)**2))
    write(*,'(3(a,i0),a,es9.2)') 'hypre: num_itr = ', num_itr, ' (', num_dscg_itr, &
        ', ', num_pcg_itr, '), ||r||_2/||b||_2 =', norm
    write(*,'(2(a,es9.2))') 'test:  ||x-x_exact||_max =', maxerr, ', ||x-x_exact||_2 =', l2err

    if (num_dscg_itr /= 19) then
      write(*,'(a,i0)') 'error: expected 19 diagonally-scaled GMRES iterations; got ', num_dscg_itr
      status = 1
    end if

    if (num_pcg_itr /= 0) then
      write(*,'(a,i0)') 'error: expected 0 AMG preconditioned GMRES iterations; got ', num_pcg_itr
      status = 1
    end if

    if (maxerr > 2.0e-6) then
      write(*,'(a,es9.2)') 'error: expected maxerr < 2.0e-6; got', maxerr
      status = 1
    end if

  end subroutine


  subroutine gmres_test_2

    type(csr_matrix), target :: matrix
    type(parameter_list), pointer :: params
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: x(:), b(:), u(:)
    integer :: nrow, num_itr, num_dscg_itr, num_pcg_itr, stat
    real(r8) :: maxerr, l2err, norm

    write(*,'(/,a)') 'Running GMRES_TEST_2'

    a = 1.0e-4_r8
    nx = 9; ny = 9; nz = 9

    allocate(params)
    call params%set ('krylov-method', 'gmres')
    call params%set ('rel-tol', 1.0e-6_r8)
    call params%set ('conv-rate-tol', 0.6_r8)
    call params%set ('max-ds-iter', 50)
    call params%set ('max-amg-iter', 10)
    call params%set ('gmres-krylov-dim', 5)
    call params%set ('amg-smoothing-sweeps', 1)
    call params%set ('logging-level', 1)
    call params%set ('print-level', 1)

    call create_matrix (matrix)
    call solver%init (matrix, params)
    call solver%setup

    nrow = matrix%nrow
    allocate(u(nrow), b(nrow), x(nrow))

    call eigenvector (1, 1, 1, u)
    b = eigenvalue(1, 1, 1) * u
    x = 0.0_r8; x(1) = 1.0_r8 ! seed the initial residual with all eigen-components
    call solver%solve (b, x, stat)
    call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
    if (stat /= 0) then
      write(*,'(a)') 'error: solver failed to converge'
      status = 1
    end if

    deallocate(params)

    maxerr = maxval(abs(u-x))
    l2err = sqrt(sum((u-x)**2))
    write(*,'(3(a,i0),a,es9.2)') 'hypre: num_itr = ', num_itr, ' (', num_dscg_itr, &
        ', ', num_pcg_itr, '), ||r||_2/||b||_2 =', norm
    write(*,'(2(a,es9.2))') 'test:  ||x-x_exact||_max =', maxerr, ', ||x-x_exact||_2 =', l2err

    if (num_dscg_itr /= 5) then
      write(*,'(a,i0)') 'error: expected 5 diagonally-scaled GMRES iterations; got ', num_dscg_itr
      status = 1
    end if

    if (num_pcg_itr /= 4) then
      write(*,'(a,i0)') 'error: expected 4 AMG preconditioned GMRES iterations; got ', num_pcg_itr
      status = 1
    end if

    if (maxerr > 1.0e-6) then
      write(*,'(a,es9.2)') 'error: expected maxerr < 1.0e-6; got', maxerr
      status = 1
    end if

  end subroutine


  subroutine create_matrix (matrix)

    type(csr_matrix), intent(out) :: matrix

    integer :: ix, iy, iz, n, j, k, ntot
    integer, allocatable :: nnbr(:,:)
    type(csr_graph), pointer :: graph

    !! Stencil neighbors of each grid point (GLOBAL).
    ntot = NX*NY*NZ
    allocate(nnbr(6,ntot))
    do iz = 0, NZ-1
      do iy = 0, NY-1
        do ix = 0, NX-1
          n = linear_index(ix,iy,iz)
          nnbr(1,n) = linear_index(ix,iy,iz-1)
          nnbr(2,n) = linear_index(ix,iy-1,iz)
          nnbr(3,n) = linear_index(ix-1,iy,iz)
          nnbr(4,n) = linear_index(ix+1,iy,iz)
          nnbr(5,n) = linear_index(ix,iy+1,iz)
          nnbr(6,n) = linear_index(ix,iy,iz+1)
          !print *, n, ':', nnbr(:,n)
        end do
      end do
    end do

    !! Create the CSR matrix.
    allocate(graph)
    call graph%init (ntot)
    do j = 1, size(nnbr,2)
      call graph%add_edge (j, j)
      call graph%add_edge (j, nnbr(:,j))
    end do
    call graph%add_complete
    call matrix%init (graph, take_graph=.true.)
    do j = 1, size(nnbr,2)
      call matrix%set (j, j, a + 6.0_r8)
      do k = 1, size(nnbr,1)
        call matrix%set (j, nnbr(k,j), -1.0_r8)
      end do
    end do

  end subroutine create_matrix

  integer function linear_index (ix, iy, iz)
    integer, intent(in) :: ix, iy, iz
    linear_index = 1 + modulo(ix,NX) + NX*(modulo(iy,NY) + NY*modulo(iz,NZ))
  end function linear_index

  function eigenvalue (kx, ky, kz) result (lambda)
    integer, intent(in) :: kx, ky, kz
    real(r8) :: lambda
    lambda = a + 4*(sin(PI*kx/NX)**2 + sin(PI*ky/NY)**2 + sin(PI*kz/NZ)**2)
  end function

  subroutine eigenvector (kx, ky, kz, u)

    integer, intent(in) :: kx, ky, kz
    real(r8), intent(out) :: u(:)

    integer :: ix, iy, iz, n
    real(r8) :: dx, dy, dz

    dx = real(kx,kind=r8) / NX
    dy = real(ky,kind=r8) / NY
    dz = real(kz,kind=r8) / NZ

    do iz = 0, NZ-1
      do iy = 0, NY-1
        do ix = 0, NX-1
          n = linear_index(ix,iy,iz)
          u(n) = cos(TWOPI*(dx*ix + dy*iy + dz*iz))
          !print *, n, ':', u(n)
        end do
      end do
    end do

  end subroutine eigenvector

end program test_hypre_hybrid_type
