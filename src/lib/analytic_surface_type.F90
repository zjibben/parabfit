!!
!! analytic_surface
!!
!! This module defines a class for an analytic implicit surface description.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! December 2016
!!

#include "f90_assert.fpp"

module analytic_surface_type

  use kinds, only: r8
  !use consts, only: ndim
  use logging_services
  implicit none
  private

  type, public :: analytic_surface
    private
    ! perhaps use a mpoly_scalar_func_type here
    real(r8), allocatable :: coeff(:)
  contains
    procedure :: init
    procedure, private :: l
    procedure, private :: Dl
    procedure :: curvature
    procedure :: Fstr
    procedure, private :: print_uf
    procedure, private :: print_f
    generic :: write(unformatted) => print_uf
    generic :: write(formatted) => print_f
  end type analytic_surface

contains

  ! perform Taubin's method (1991) to calculate the analytic surface which fits the points x(:,:)
  ! in a least-squares sense
  subroutine init (this,x)

    use array_utils, only: outer_product,isZero, normalize
    use, intrinsic :: iso_c_binding, only: c_new_line

    class(analytic_surface), intent(out) :: this
    real(r8), intent(in) :: x(:,:)

    real(r8), allocatable :: M(:,:), N(:,:), lr(:), li(:), vr(:,:), tmpM(:,:), tmpV(:), tmpV2(:)
    real(r8) :: tmpR
    integer :: ierr,s,i

    ! Exclusively 3D for now. Adding a 2D-specific version would be trivial, though.
    ASSERT(size(x,1)==3)

    s = 10 ! number of terms for a quadratic function

    ! build the matrices for the generalized eigen-problem
    allocate(this%coeff(s), M(s,s), N(s,s))
    M = 0.0_r8; N = 0.0_r8
    do i = 1,size(x,2)
      tmpV = this%l(x(:,i))
      tmpM = this%Dl(x(:,i))
      M = M + outer_product(tmpV, tmpV)
      N = N + matmul(tmpM,transpose(tmpM))
    end do

    ! print *, c_new_line, "M:"
    ! do i = 1,s
    !   print '(10es10.2)', M(:,i)
    ! end do

    ! print *, c_new_line, "N:"
    ! do i = 1,s
    !   print '(10es10.2)', N(:,i)
    ! end do

    ! solve the generalized eigenvector problem M*vr = lr*N*vr
    ! note could run this twice, first time just getting the ideal work size (tmpV2)
    allocate(lr(s),li(s),vr(s,s), tmpV2(8*s))
    call dggev ('N','V', s, M, s, N, s, &
        lr,li, tmpV, tmpM, s, vr, s, &
        tmpV2, 8*s, ierr)
    if (ierr/=0) call LS_fatal ('failed generalized eigenvalue solve')
    lr = lr/tmpV

    ! pick the eigenvector corresponding to the smallest real eigenvalue
    this%coeff = vr(:,minloc(lr, dim=1, mask=isZero(li)))

    ! print *, c_new_line, "evs:"
    ! do i = 1,s
    !   print '(10es10.2)', lr(i), li(i)
    !   print '(10es10.2)', normalize(vr(:,i))
    !   print *, c_new_line
    ! end do

    ! print *, c_new_line, "coeffs: "
    ! print '(10es10.2)', normalize(this%coeff)

    ! tmpR = minval(lr, dim=1, mask=isZero(li))
    ! print *, tmpR
    ! print '(10es10.2)', matmul(M - tmpR*N, this%coeff)

    ! print *, c_new_line

  end subroutine init

  ! return the vector of terms (without coefficients), (1, x, y, z, x**2, etc), for a point x(:)
  function l (this,x)
    class(analytic_surface), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8), allocatable :: l(:)
    ASSERT(size(x)==3)
    l = [1.0_r8, x(1), x(2), x(3), x(1)**2, x(1)*x(2), x(1)*x(3), x(2)**2, x(2)*x(3), x(3)**2]
  end function l

  ! return the Jacobian at point x(:)
  function Dl (this,x)

    class(analytic_surface), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8), allocatable :: Dl(:,:)

    ASSERT(size(x)==3)

    Dl = reshape([&
        0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 2.0_r8*x(1), x(2),   x(3),   0.0_r8,      0.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8,      x(1),   0.0_r8, 2.0_r8*x(2), x(3),   0.0_r8, &
        0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8,      0.0_r8, x(1),   0.0_r8,      x(2),  2.0_r8*x(3)&
        ], [10,3])

  end function Dl

  ! calculate the curvature at a point x
  ! this formula is generated from a Python+Sympy script
  real(r8) function curvature (this,x)

    use consts, only: alittle
    use array_utils, only: clip

    class(analytic_surface), intent(in) :: this
    real(r8), intent(in) :: x(:)

    ASSERT(size(x)==3)

    curvature = &
        ((this%coeff(2) + 2*this%coeff(5)*x(1) + this%coeff(6)*x(2) + this%coeff(7)*x(3))**2 + &
        (this%coeff(3) + this%coeff(6)*x(1) + 2*this%coeff(8)*x(2) + this%coeff(9)*x(3))**2 + &
        (this%coeff(4) + this%coeff(7)*x(1) + this%coeff(9)*x(2) + &
        2*this%coeff(10)*x(3))**2)**(-1.5)*(-2*(this%coeff(5) + this%coeff(8) + &
        this%coeff(10))*((this%coeff(2) + 2*this%coeff(5)*x(1) + this%coeff(6)*x(2) + &
        this%coeff(7)*x(3))**2 + (this%coeff(3) + this%coeff(6)*x(1) + 2*this%coeff(8)*x(2) + &
        this%coeff(9)*x(3))**2 + (this%coeff(4) + this%coeff(7)*x(1) + this%coeff(9)*x(2) + &
        2*this%coeff(10)*x(3))**2) + (2*this%coeff(5)*(this%coeff(2) + 2*this%coeff(5)*x(1) + &
        this%coeff(6)*x(2) + this%coeff(7)*x(3)) + this%coeff(6)*(this%coeff(3) + &
        this%coeff(6)*x(1) + 2*this%coeff(8)*x(2) + this%coeff(9)*x(3)) + &
        this%coeff(7)*(this%coeff(4) + this%coeff(7)*x(1) + this%coeff(9)*x(2) + &
        2*this%coeff(10)*x(3)))*(this%coeff(2) + 2*this%coeff(5)*x(1) + this%coeff(6)*x(2) + &
        this%coeff(7)*x(3)) + (this%coeff(6)*(this%coeff(2) + 2*this%coeff(5)*x(1) + &
        this%coeff(6)*x(2) + this%coeff(7)*x(3)) + 2*this%coeff(8)*(this%coeff(3) + &
        this%coeff(6)*x(1) + 2*this%coeff(8)*x(2) + this%coeff(9)*x(3)) + &
        this%coeff(9)*(this%coeff(4) + this%coeff(7)*x(1) + this%coeff(9)*x(2) + &
        2*this%coeff(10)*x(3)))*(this%coeff(3) + this%coeff(6)*x(1) + 2*this%coeff(8)*x(2) + &
        this%coeff(9)*x(3)) + (this%coeff(7)*(this%coeff(2) + 2*this%coeff(5)*x(1) + &
        this%coeff(6)*x(2) + this%coeff(7)*x(3)) + this%coeff(9)*(this%coeff(3) + &
        this%coeff(6)*x(1) + 2*this%coeff(8)*x(2) + this%coeff(9)*x(3)) + &
        2*this%coeff(10)*(this%coeff(4) + this%coeff(7)*x(1) + this%coeff(9)*x(2) + &
        2*this%coeff(10)*x(3)))*(this%coeff(4) + this%coeff(7)*x(1) + this%coeff(9)*x(2) + &
        2*this%coeff(10)*x(3)))/2

    curvature = clip(curvature, 1e10_r8, 0.0_r8)

  end function curvature

  function Fstr (this)

    use array_utils, only: isZero

    class(analytic_surface), intent(in) :: this
    character(:), allocatable :: Fstr

    character(:), allocatable :: terms(:)
    character(32) :: term_str
    integer :: i

    terms = ['1','x','y','z','x**2','x*y','x*z','y**2','y*z','z**2']
    Fstr = ''

    do i = 1,size(this%coeff)
      if (.not.isZero(this%coeff(i),1e-5_r8)) then
        if (this%coeff(i) > 0.0_r8) then
          write(term_str,'(a,es10.3,2a)') '  + ',this%coeff(i),'*',terms(i)
        else
          write(term_str,'(a,es10.3,2a)') '  - ',abs(this%coeff(i)),'*',terms(i)
        end if
        Fstr = Fstr // trim(term_str)
      end if
    end do

  end function Fstr

  subroutine print_f (this, unit, iotype, vlist, iostat, iomsg)

    class(analytic_surface), intent(in) :: this
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: vlist(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    write (unit, '(2a)', iostat=iostat) this%Fstr(), ' = 0'

  end subroutine print_f

  subroutine print_uf (this, unit, iostat, iomsg)

    class(analytic_surface), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    write (unit, '(2a)', iostat=iostat) this%Fstr(), ' = 0'

  end subroutine print_uf

end module analytic_surface_type
