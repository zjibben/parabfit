!!
!! paraboloid
!!
!! This module defines a class for an analytic paraboloid.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! February 2017
!!

#include "f90_assert.fpp"

module paraboloid_type

  use kinds, only: r8
  !use consts, only: ndim
  use logging_services
  implicit none
  private

  type, public :: paraboloid
    private
    ! perhaps use a mpoly_scalar_func_type here
    real(r8), allocatable :: coeff(:), cr(:)
    logical :: initialized
  contains
    procedure :: init
    procedure :: bestFit
    !procedure :: canonicalForm
    procedure, private :: taubinCoeffs
    procedure, private :: coeffsinOriginalSpace
    procedure, private :: l
    procedure, private :: llong
    procedure, private :: Dl
    procedure :: curvature
    procedure :: Fstr
    procedure, private :: print_uf
    procedure, private :: print_f
    generic :: write(unformatted) => print_uf
    generic :: write(formatted) => print_f
  end type paraboloid

  real(r8), parameter :: eccentricity_max = 1e4_r8

contains

  ! perform Taubin's method (1991) to calculate the analytic surface which fits the points x(:,:)
  ! in a least-squares sense
  subroutine init (this, coeff, coeff_rotated)
    class(paraboloid), intent(out) :: this
    real(r8), intent(in) :: coeff(:), coeff_rotated(:)
    this%coeff = coeff
    this%cr = coeff_rotated
    this%initialized = .true.
  end subroutine init

  subroutine bestFit (this, x)

    use array_utils, only: normalize, crossProduct, isZero

    class(paraboloid), intent(out) :: this
    real(r8), intent(in) :: x(:,:)

    real(r8), allocatable :: lr(:), vr(:,:), c(:), cr(:)
    real(r8) :: xcen(3), R(3,3), normal(3), xr(size(x,dim=1), size(x,dim=2))
    integer :: i

    ASSERT(size(x,1)==3) ! Exclusively 3D for now. Adding a 2D version would be trivial, though.

    ! get the center and normal vector
    xcen = sum(x(:,1:3), dim=2) / 3.0_r8
    normal = normalize(crossProduct(x(:,2) - x(:,1), x(:,3) - x(:,1)))

    ! set up the rotation matrix
    R(3,:) = normal
    R(2,:) = normalize(crossProduct([0.0_r8,0.0_r8,1.0_r8], normal))
    R(1,:) = crossProduct(R(2,:), normal)

    ! get the set of points in the rotated and translated coordinate space
    do i = 1,size(x, dim=2)
      xr(:,i) = matmul(R, x(:,i) - xcen)
    end do
    
    ! print *, "xc: ", xcen
    ! print *, "n:  ", normal
    ! print *, "x1: ", xr(:,1)
    ! print *

    ! fit a paraboloid to xr
    call this%taubinCoeffs (lr,vr, xr)
    cr = vr(:,minloc(lr, dim=1)) ! coefficients associated with the smallest eigenvalue

    ! make sure the z coefficient is nonzero
    if (isZero(cr(4))) call LS_fatal ('z term zero')
    cr = - cr / cr(4)

    ! print *, minval(lr)
    ! print *
    ! print *, lr
    ! print *
    ! print *, cr

    ! get the coefficients in the original space
    c = this%coeffsInOriginalSpace(cr, R, xcen)

    call this%init (c, cr)
    
  end subroutine bestFit

  function coeffsinOriginalSpace (this, c, R, xcen) result(cp)

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: c(:), R(:,:), xcen(:)
    real(r8) :: cp(10)

    real(r8) :: Ar(3,3), br(3), A(3,3), b(3)

    ! extend coefficient array
    cp = 0.0_r8
    cp(1:6) = c(1:6)
    cp(8) = c(7)

    Ar = coeffs2matrix(cp)
    br = cp(1:3)

    ! convert to original space
    A = matmul(transpose(R), matmul(Ar, R))
    b = matmul(transpose(R), br) - 2.0_r8 * matmul(A, xcen)
    cp(1) = dot_product(xcen, matmul(A,xcen)) - dot_product(matmul(transpose(R),br), xcen) + cp(1)

    cp = matvec2coeffs(A, b, cp(1))
    
  end function coeffsinOriginalSpace

  function coeffs2matrix (c)
    real(r8), intent(in) :: c(:)
    real(r8) :: coeffs2matrix(3,3)
    coeffs2matrix = reshape([&
        c(5), c(6)/2.0_r8, c(7)/2.0_r8,&
        c(6)/2.0_r8, c(8), c(9)/2.0_r8,&
        c(7)/2.0_r8, c(9)/2.0_r8, c(10)], [3,3])
  end function coeffs2matrix

  function matvec2coeffs (A, b, k) result(c)
    real(r8), intent(in) :: A(:,:), b(:), k
    real(r8) :: c(10)
    c(1) = k
    c(2:4) = b
    c(5)  = A(1,1)
    c(6)  = 2.0_r8 * A(1,2)
    c(7)  = 2.0_r8 * A(1,3)
    c(8)  = A(2,2)
    c(9)  = 2.0_r8 * A(2,3)
    c(10) = A(3,3)
  end function matvec2coeffs
  
  subroutine taubinCoeffs (this, evl, evc, x)
    
    use array_utils, only: outer_product, isZero
    use, intrinsic :: iso_c_binding, only: c_new_line
    
    class(paraboloid), intent(in) :: this
    real(r8), allocatable, intent(out) :: evl(:), evc(:,:)
    real(r8), intent(in) :: x(:,:)

    real(r8), allocatable :: M(:,:), N(:,:), lr(:), li(:), vr(:,:), tmpM(:,:), tmpV(:), tmpV2(:)
    real(r8) :: tmpR
    integer :: ierr,s,i,j
    
    ! Exclusively 3D for now. Adding a 2D-specific version would be trivial, though.
    ASSERT(size(x,1)==3)

    s = 7 ! number of terms for a paraboloid

    ! build the matrices for the generalized eigen-problem
    allocate(M(s,s), N(s,s))
    M = 0.0_r8; N = 0.0_r8
    do i = 1,size(x,2)
      tmpV = this%l(x(:,i))
      tmpM = this%Dl(x(:,i))
      M = M + outer_product(tmpV, tmpV)
      N = N + matmul(tmpM,transpose(tmpM))
    end do

    ! solve the generalized eigenvector problem M*vr = lr*N*vr
    ! note could run this twice, first time just getting the ideal work size (tmpV2)
    allocate(lr(s),li(s),vr(s,s), tmpV2(8*s))
    call dggev ('N','V', s, M, s, N, s, &
        lr,li, tmpV, tmpM, s, vr, s, &
        tmpV2, 8*s, ierr)
    if (ierr/=0) call LS_fatal ('failed generalized eigenvalue solve')

    ! only return real results
    !li = li/tmpV
    evl = pack(lr/tmpV, mask=isZero(li) .and. .not.isZero(tmpV))

    allocate(evc(s,size(evl)))
    j = 1
    do i = 1,s
      if (isZero(li(i)) .and. .not.isZero(tmpV(i))) then
        evc(:,j) = vr(:,i)
        j = j + 1
      end if
    end do

    ! print *, c_new_line, "M:"
    ! do i = 1,s
    !   print '(10es10.2)', M(:,i)
    ! end do
    
    ! print *, c_new_line, "N:"
    ! do i = 1,s
    !   print '(10es10.2)', N(:,i)
    ! end do
    
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

  end subroutine taubinCoeffs

  ! return the vector of terms (without coefficients), (1, x, y, z, x**2, etc), for a point x(:)
  function l (this,x)
    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8), allocatable :: l(:)
    ASSERT(size(x)==3)
    l = [1.0_r8, x(1), x(2), x(3), x(1)**2, x(1)*x(2), x(2)**2]
  end function l

  function llong (this,x)
    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8), allocatable :: llong(:)
    ASSERT(size(x)==3)
    llong = [1.0_r8, x(1), x(2), x(3), x(1)**2, x(1)*x(2), x(1)*x(3), x(2)**2, x(2)*x(3), x(3)**2]
  end function llong

  ! return the Jacobian at point x(:)
  function Dl (this,x)

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8), allocatable :: Dl(:,:)

    ASSERT(size(x)==3)

    Dl = reshape([&
        0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8, 2.0_r8*x(1), x(2),   0.0_r8, &
        0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8, 0.0_r8,      x(1),   2.0_r8*x(2),&
        0.0_r8, 0.0_r8, 0.0_r8, 1.0_r8, 0.0_r8,      0.0_r8, 0.0_r8 &
        ], [7,3])

  end function Dl

  ! calculate the curvature at a point x
  real(r8) function curvature (this,x)

    use consts, only: alittle
    use array_utils, only: clip

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)

    ASSERT(size(x)==3)

    curvature = 2.0_r8 * &
        (this%cr(5) + this%cr(7) + this%cr(5)*this%cr(3)**2 + this%cr(7)*this%cr(2)**2 - &
        this%cr(6)*this%cr(2)*this%cr(3)) / (1.0_r8 + this%cr(2)**2 + this%cr(3)**2)**1.5_r8
    
    curvature = clip(curvature, 1e10_r8, 0.0_r8)

  end function curvature

  function Fstr (this)

    use array_utils, only: isZero

    class(paraboloid), intent(in) :: this
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

    class(paraboloid), intent(in) :: this
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: vlist(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    write (unit, '(2a)', iostat=iostat) this%Fstr(), ' = 0'

  end subroutine print_f

  subroutine print_uf (this, unit, iostat, iomsg)

    class(paraboloid), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    write (unit, '(2a)', iostat=iostat) this%Fstr(), ' = 0'

  end subroutine print_uf

end module paraboloid_type
