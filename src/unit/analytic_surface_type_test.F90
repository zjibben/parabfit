module analytic_surface_type_test

  use kinds, only: r8
  use analytic_surface_type
  implicit none
  private

  public :: analytic_surface_test_suite

contains

  subroutine analytic_surface_test_suite ()

    print '(a)'
    print '(a)', 'ANALYTIC_SURFACE'
    print '(a)', '===================================================='

    call plane_test ()
    call parabola_test ()

    print '(a)', '===================================================='
    print '(a)'

  end subroutine analytic_surface_test_suite

  subroutine plane_test ()
    
    real(r8), allocatable :: x(:,:)
    real(r8) :: dx
    integer :: N,ind,i,j
    type(analytic_surface) :: surf

    dx = 0.1_r8
    N = 9

    allocate(x(3,N*N))
    do i = 1,N
      do j = 1,N
        ind = i + (j-1)*N
        x(1,ind) = real(i-N/2+1,r8)*dx
        x(2,ind) = real(j-N/2+1,r8)*dx
        x(3,ind) = x(1,ind) + x(2,ind)
      end do
    end do
    
    call surf%init (x)

    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature([0.0_r8,0.0_r8,0.0_r8])
    
  end subroutine plane_test

  subroutine parabola_test ()
    
    real(r8), allocatable :: x(:,:)
    real(r8) :: dx
    integer :: N,ind,i,j
    type(analytic_surface) :: surf

    dx = 0.1_r8
    N = 9

    allocate(x(3,N*N))
    do i = 1,N
      do j = 1,N
        ind = i + (j-1)*N
        x(1,ind) = real(i-N/2+1,r8)*dx
        x(2,ind) = real(j-N/2+1,r8)*dx
        x(3,ind) = x(1,ind)**2 + x(2,ind)**2 ! parabola
      end do
    end do
    
    call surf%init (x)

    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature([0.0_r8,0.0_r8,0.0_r8])

  end subroutine parabola_test
    
end module analytic_surface_type_test
