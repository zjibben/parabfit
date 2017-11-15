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
    real(r8), public :: rot(3,3), offset(3)
    logical :: initialized
    integer :: lwork
  contains
    procedure :: init
    procedure :: bestFit
    procedure :: volumetricFit
    !procedure :: canonicalForm
    procedure, private :: taubinCoeffs
    procedure, private :: coeffsinOriginalSpace
    procedure, private :: paraboloidCoeffsDirect
    procedure, private :: l
    procedure, private :: llong
    procedure, private :: l_paraboloid
    procedure, private :: Dl
    procedure, private :: optimalWorkSize
    procedure, private :: localCoords
    procedure, private :: globalCoords
    procedure :: curvature
    procedure :: normal
    procedure :: normal_average
    procedure :: curvatureQdrtc
    procedure :: Fstr
    procedure :: fstr_rot
    procedure :: point_on_surface
    procedure :: point_along_line
    procedure, private :: print_uf
    procedure, private :: print_f
#ifndef NAGFOR
    ! NAG Compiler does not support this as of v6.1
    generic :: write(unformatted) => print_uf
    generic :: write(formatted) => print_f
#endif
  end type paraboloid

  real(r8), parameter :: eccentricity_max = 1e4_r8

contains

  ! perform Taubin's method (1991) to calculate the analytic surface which fits the points x(:,:)
  ! in a least-squares sense
  subroutine init (this, coeff, coeff_rotated)

    class(paraboloid), intent(out) :: this
    real(r8), intent(in) :: coeff(:), coeff_rotated(:)

    real(r8) :: work(1), tmp1(1,1), tmp2(1), tmp3(1)
    integer :: ierr

    this%coeff = coeff
    this%cr = coeff_rotated
    this%initialized = .true.

    this%lwork = this%optimalWorkSize(6)

  end subroutine init

  subroutine bestFit (this, x, weight, normal)

    use array_utils, only: rotationMatrix, normalize, crossProduct
    use logging_services

    class(paraboloid), intent(out) :: this
    real(r8), intent(in) :: x(:,:), weight(:)
    real(r8), intent(in), optional :: normal(:)

    real(r8), allocatable :: lr(:), vr(:,:), c(:), cr(:)
    real(r8) :: xcen(3), R(3,3), normalh(3), xr(size(x,dim=1), size(x,dim=2))
    integer :: i

    ASSERT(size(x,1)==3) ! Exclusively 3D for now. Adding a 2D version would be trivial, though.
    ASSERT(size(x, dim=2) == size(weight))

    ! get the center and normal vector
    !xcen = sum(x(:,1:3), dim=2) / 3.0_r8
    xcen = x(:,1)
    if (present(normal)) then
      normalh = normal
    else
      ! this is valid only when we passed in 3 points per polygon
      normalh = normalize(crossProduct(x(:,2) - x(:,1), x(:,3) - x(:,1)))
      call LS_fatal ("disabling cross product normals for now")
    end if

    ! set up the rotation matrix
    R = rotationMatrix(normalh)

    ! get the set of points in the rotated and translated coordinate space
    do i = 1,size(x, dim=2)
      xr(:,i) = matmul(R, x(:,i) - xcen)
    end do

    ! ! fit a paraboloid to xr
    ! call this%taubinCoeffs (lr,vr, xr)
    ! cr = vr(:,minloc(lr, dim=1)) ! coefficients associated with the smallest eigenvalue

    ! ! make sure the z coefficient is nonzero
    ! if (isZero(cr(4))) call LS_fatal ('z term zero')
    ! cr = - cr / cr(4)

    cr = this%paraboloidCoeffsDirect (xr, weight)

    ! print '(a)', this%Fstr_rot(cr)
    ! do i = 1,size(xr,dim=2)
    !   print '(a,3es14.4)', 'xr: ', xr(:,i)
    ! end do
    ! print *

    ! print *, minval(lr)
    ! print *
    ! print *, lr
    ! print *
    ! print *, cr

    ! get the coefficients in the original space
    c = this%coeffsInOriginalSpace(cr, R, xcen)

    call this%init (c, cr)

    this%rot = R
    this%offset = xcen

  end subroutine bestFit

  ! subroutine volumetricFit (this, interface_reconstruction)

  !   use polygon_type
  !   use array_utils, only: outer_product, rotationMatrix
  !   use logging_services
  !   external dsysv ! lapack

  !   class(paraboloid), intent(out) :: this
  !   type(polygon), intent(in) :: interface_reconstruction(:)

  !   type(polygon) :: int_rec
  !   real(r8), allocatable :: c(:), cr(:)
  !   real(r8) :: xcen(3), normal(3), R(3,3)
  !   real(r8) :: integrals(6), reconstruction_plane_coeffs(3), xv, yv, xvn, yvn, A(6,6), b(6)
  !   integer :: i, n_rec, ierr, v, vn, ipiv(6), lwork
  !   real(r8), allocatable :: work(:)

  !   xcen = interface_reconstruction(1)%centroid2()
  !   normal = interface_reconstruction(1)%norm

  !   ! get the set of points in the rotated and translated coordinate space
  !   A = 0; b = 0
  !   do i = 1,size(interface_reconstruction)
  !     int_rec = interface_reconstruction(i)
  !     call int_rec%rotate_offset(normal, xcen)

  !     reconstruction_plane_coeffs = [-dot_product(int_rec%centroid2(), int_rec%norm), &
  !         int_rec%norm(1), int_rec%norm(2)] / (-int_rec%norm(3))

  !     integrals = 0
  !     do v = 1,int_rec%nVerts
  !       vn = modulo(v,int_rec%nVerts) + 1
  !       xv = int_rec%x(1,v)
  !       yv = int_rec%x(2,v)
  !       xvn = int_rec%x(1,vn)
  !       yvn = int_rec%x(2,vn)

  !       integrals = integrals + [&
  !           (xv*yvn - xvn*yv) / 2, &
  !           (xv + xvn)*(xv*yvn - xvn*yv) / 6, &
  !           (yv + yvn)*(xv*yvn - xvn*yv) / 6, &
  !           (xv + xvn)*(xv**2 + xvn**2)*(yvn - yv) / 12, &
  !           (yvn - yv)*(3*xv**2*yv + xv**2*yvn + 2*xv*xvn*yv + 2*xv*xvn*yvn + xvn**2*yv + 3*xvn**2*yvn)/24, &
  !           (xv - xvn)*(yv + yvn)*(yv**2 + yvn**2) / 12]
  !     end do

  !     A = A + outer_product(integrals, integrals)
  !     b = b + integrals * dot_product(integrals(:3), reconstruction_plane_coeffs)
  !   end do

  !   ! get optimal work size
  !   lwork = this%optimalWorkSize(6)
  !   allocate(work(lwork))

  !   ! solve the symmetric linear system (note lapack puts the solution in b)
  !   call dsysv ('U', 6, 1, A, 6, ipiv, b, 6, work, lwork, ierr)
  !   if (ierr /= 0) call LS_fatal ('volumetricFit: failed symmetric linear solve')

  !   ! convert to full 3d form
  !   cr = [b(1), b(2), b(3), -1.0_r8, b(4), b(5), b(6)]

  !   ! get coeffs in original space
  !   R = rotationMatrix(normal)
  !   c = this%coeffsInOriginalSpace(cr, R, xcen)

  !   ! generate paraboloid object
  !   call this%init (c, cr)
  !   this%rot = R
  !   this%offset = xcen

  ! end subroutine volumetricFit

  subroutine volumetricFit (this, interface_reconstruction_collection)

    use polygon_type
    use array_utils, only: outer_product, rotationMatrix
    use logging_services
    external dsysv ! lapack

    class(paraboloid), intent(out) :: this
    type(polygon_box), intent(in) :: interface_reconstruction_collection(:)

    type(polygon) :: int_rec
    real(r8), allocatable :: c(:), cr(:)
    real(r8) :: xcen(3), normal(3), rot(3,3)
    real(r8) :: integrals(6), reconstruction_plane_coeffs(3), xv, yv, xvn, yvn, A(6,6), b(6), &
        integral_sum(6), b_dot_sum
    integer :: i, n_rec, ierr, v, vn, ipiv(6), lwork, p, r
    real(r8), allocatable :: work(:)

    ! xcen = 0
    ! do r = 1,interface_reconstruction_collection(1)%n_elements
    !   xcen = xcen + interface_reconstruction_collection(1)%elements(r)%centroid2()
    ! end do
    ! xcen = xcen / interface_reconstruction_collection(1)%n_elements
    xcen = interface_reconstruction_collection(1)%elements(1)%centroid2()
    normal = interface_reconstruction_collection(1)%elements(1)%norm

    ! get the set of points in the rotated and translated coordinate space
    A = 0; b = 0
    do p = 1,size(interface_reconstruction_collection)

      integral_sum = 0; b_dot_sum = 0
      do r = 1,interface_reconstruction_collection(p)%n_elements
        int_rec = interface_reconstruction_collection(p)%elements(r)
        call int_rec%rotate_offset(normal, xcen)

        reconstruction_plane_coeffs = [-dot_product(int_rec%centroid2(), int_rec%norm), &
            int_rec%norm(1), int_rec%norm(2)] / (-int_rec%norm(3))

        integrals = 0
        do v = 1,int_rec%nVerts
          vn = modulo(v,int_rec%nVerts) + 1
          xv = int_rec%x(1,v)
          yv = int_rec%x(2,v)
          xvn = int_rec%x(1,vn)
          yvn = int_rec%x(2,vn)

          integrals = integrals + [&
              (xv*yvn - xvn*yv) / 2, &
              (xv + xvn)*(xv*yvn - xvn*yv) / 6, &
              (yv + yvn)*(xv*yvn - xvn*yv) / 6, &
              (xv + xvn)*(xv**2 + xvn**2)*(yvn - yv) / 12, &
              (yvn - yv)*(3*xv**2*yv + xv**2*yvn + 2*xv*xvn*yv + 2*xv*xvn*yvn + xvn**2*yv + 3*xvn**2*yvn)/24, &
              (xv - xvn)*(yv + yvn)*(yv**2 + yvn**2) / 12]
        end do

        integral_sum = integral_sum + integrals
        b_dot_sum = b_dot_sum + dot_product(integrals(:3), reconstruction_plane_coeffs)
      end do

      A = A + outer_product(integral_sum, integral_sum)
      b = b + integral_sum * b_dot_sum
    end do

    ! get optimal work size
    lwork = this%optimalWorkSize(6)
    allocate(work(lwork))

    ! solve the symmetric linear system (note lapack puts the solution in b)
    call dsysv ('U', 6, 1, A, 6, ipiv, b, 6, work, lwork, ierr)
    if (ierr /= 0) call LS_fatal ('volumetricFit: failed symmetric linear solve')

    ! convert to full 3d form
    cr = [b(1), b(2), b(3), -1.0_r8, b(4), b(5), b(6)]

    ! get coeffs in original space
    rot = rotationMatrix(normal)
    c = this%coeffsInOriginalSpace(cr, rot, xcen)

    ! generate paraboloid object
    call this%init (c, cr)
    this%rot = rot
    this%offset = xcen

  end subroutine volumetricFit

  ! take a set of points in a rotated coordinate space
  ! and return the best fitting paraboloid
  function paraboloidCoeffsDirect (this, x, weight) result(c)

    use array_utils, only: outer_product
    external dsysv ! lapack symmetric linear solve

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:,:), weight(:)
    real(r8), allocatable :: c(:)

    integer :: i, s, ierr, lwork
    integer, allocatable :: ipiv(:)
    real(r8), allocatable :: A(:,:), b(:), tmpl(:), work(:)

    ASSERT(size(x, dim=1) == 3)
    ASSERT(size(x, dim=2) == size(weight))

    s = 6 ! number of terms for direct paraboloid z = f(x)

    ! set up lhs & rhs
    allocate(A(s,s), b(s), ipiv(s), c(s+1))
    A = 0; b = 0
    do i = 1,size(x, dim=2)
      tmpl = this%l_paraboloid(x(:,i))
      A = A + weight(i)**2 * outer_product(tmpl,tmpl)
      b = b + weight(i)**2 * x(3,i) * tmpl
    end do

    ! get optimal work size
    lwork = this%optimalWorkSize(s)
    allocate(work(lwork))

    ! solve the linear system
    ! note lapack puts the solution in b
    call dsysv ('U', s, 1, A, s, ipiv, b, s, work, lwork, ierr)
    if (ierr /= 0) call LS_fatal ('failed symmetric linear solve')

    ! convert coeffs to implicit form
    c(1:3) = b(1:3)
    c(4) = -1
    c(5:) = b(4:)

  end function paraboloidCoeffsDirect

  integer function optimalWorkSize (this, s)

    external dsysv ! lapack symmetric linear solve

    class(paraboloid), intent(in) :: this
    integer, intent(in) :: s

    integer :: ipiv(s), ierr
    real(r8) :: A(s,s), b(s), tmpl(s), work(1)

    call dsysv ('U', s, 1, A, s, ipiv, b, s, work, -1, ierr)
    if (ierr /= 0) call LS_error ('could not get optimal work size')

    optimalWorkSize = int(work(1))

  end function optimalWorkSize

  function coeffsInOriginalSpace (this, c, R, xcen) result(cp)

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: c(:), R(:,:), xcen(:)
    real(r8) :: cp(10)

    real(r8) :: Ar(3,3), br(3), A(3,3), b(3)

    ! extend coefficient array
    cp = 0.0_r8
    cp(1:6) = c(1:6)
    cp(8) = c(7)

    Ar = coeffs2matrix(cp)
    br = cp(2:4)

    ! convert to original space
    A = matmul(transpose(R), matmul(Ar, R))
    b = matmul(transpose(R), br) - 2 * matmul(A, xcen)
    cp(1) = dot_product(xcen, matmul(A,xcen)) - dot_product(matmul(transpose(R),br), xcen) + cp(1)

    cp = matvec2coeffs(A, b, cp(1))

  end function coeffsInOriginalSpace

  function coeffs2matrix (c)
    real(r8), intent(in) :: c(:)
    real(r8) :: coeffs2matrix(3,3)
    coeffs2matrix = reshape([&
        c(5), c(6)/2, c(7)/2,&
        c(6)/2, c(8), c(9)/2,&
        c(7)/2, c(9)/2, c(10)], [3,3])
  end function coeffs2matrix

  function matvec2coeffs (A, b, k) result(c)
    real(r8), intent(in) :: A(:,:), b(:), k
    real(r8) :: c(10)
    c(1) = k
    c(2:4) = b
    c(5)  = A(1,1)
    c(6)  = 2 * A(1,2)
    c(7)  = 2 * A(1,3)
    c(8)  = A(2,2)
    c(9)  = 2 * A(2,3)
    c(10) = A(3,3)
  end function matvec2coeffs

  subroutine taubinCoeffs (this, evl, evc, x)

    use array_utils, only: outer_product, isZero
    use, intrinsic :: iso_c_binding, only: c_new_line
    external dggev

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
    M = 0; N = 0
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

  function l_paraboloid (this, x)
    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8), allocatable :: l_paraboloid(:)
    ASSERT(size(x)==3)
    l_paraboloid = [1.0_r8, x(1), x(2), x(1)**2, x(1)*x(2), x(2)**2]
  end function l_paraboloid

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

  function point_on_surface(this, x) result(xr)

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: xr(3)

    xr = this%localCoords(x)
    xr(3) = 0
    xr(3) = sum(this%cr*this%l(xr))

    xr = this%globalCoords(xr)

  end function point_on_surface

  ! given a point and a direction, find the point on
  ! the surface on that line closest to the given point
  function point_along_line(this, xo, do) result(xr)

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: xo(:), do(:)
    real(r8) :: xr(3)

    real(r8) :: l, t, a, b, c, x(3), d(3)

    ! convert to rotated coordinate space
    x = this%localCoords(xo)
    d = matmul(this%rot, do)

    ! get quadratic coefficients
    a = this%cr(5)*d(1)**2 + this%cr(6)*d(1)*d(2) + this%cr(7)*d(2)**2
    b = this%cr(2)*d(1) + this%cr(3)*d(2) + 2*this%cr(5)*x(1)*d(1) + this%cr(6)*x(2)*d(1) &
        + this%cr(6)*x(1)*d(2) + 2*this%cr(7)*x(2)*d(2) - d(3)
    c = sum(this%cr * this%l(x))

    ! looking for the nearest point, so minimum absolute l
    l = -b + sqrt(b**2 - 4*a*c)
    t = -b - sqrt(b**2 - 4*a*c)
    if (abs(t) < abs(l)) l = t
    l = l / (2*a)

    xr = xo + l*do

  end function point_along_line

  ! calculate the curvature at a point x
  real(r8) function curvature (this,x)

    use array_utils, only: clip

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)

    real(r8) :: xr(3)

    ASSERT(size(x)==3)

    ! curvature = 2.0_r8 * &
    !     (this%cr(5) + this%cr(7) + this%cr(5)*this%cr(3)**2 + this%cr(7)*this%cr(2)**2 - &
    !     this%cr(6)*this%cr(2)*this%cr(3)) / (1.0_r8 + this%cr(2)**2 + this%cr(3)**2)**1.5_r8

    xr = this%localCoords(x)

    ! print *, curvature
    ! print *, x
    ! print *, xr

    ! curvature = 2 * (-2*(this%cr(4) + this%cr(6))*((this%cr(2) + 2*this%cr(4)*xr(1) + &
    !     this%cr(5)*xr(2))**2 + (this%cr(3) + this%cr(5)*xr(1) + 2*this%cr(6)*xr(2))**2 + 1) + &
    !     (2*this%cr(4)*(this%cr(2) + 2*this%cr(4)*xr(1) + this%cr(5)*xr(2)) + &
    !     this%cr(5)*(this%cr(3) + this%cr(5)*xr(1) + 2*this%cr(6)*xr(2)))*(this%cr(2) + &
    !     2*this%cr(4)*xr(1) + this%cr(5)*xr(2)) + (this%cr(5)*(this%cr(2) + 2*this%cr(4)*xr(1) + &
    !     this%cr(5)*xr(2)) + 2*this%cr(6)*(this%cr(3) + this%cr(5)*xr(1) + &
    !     2*this%cr(6)*xr(2)))*(this%cr(3) + this%cr(5)*xr(1) + 2*this%cr(6)*xr(2)))*((this%cr(2) + &
    !     2*this%cr(4)*xr(1) + this%cr(5)*xr(2))**2 + (this%cr(3) + this%cr(5)*xr(1) + &
    !     2*this%cr(6)*xr(2))**2 + 1)**(-1.5_r8)

    ! curvature = 2 * (-2*(this%cr(4) + this%cr(6))*((this%cr(2) + 2*this%cr(4)*xr(1) + &
    !     this%cr(5)*xr(2))**2 + (this%cr(3) + this%cr(5)*xr(1) + 2*this%cr(6)*xr(2))**2 + 1) + &
    !     (2*this%cr(4)*(this%cr(2) + 2*this%cr(4)*xr(1) + this%cr(5)*xr(2)) + &
    !     this%cr(5)*(this%cr(3) + this%cr(5)*xr(1) + 2*this%cr(6)*xr(2)))*(this%cr(2) + &
    !     2*this%cr(4)*xr(1) + this%cr(5)*xr(2)) + (this%cr(5)*(this%cr(2) + 2*this%cr(4)*xr(1) + &
    !     this%cr(5)*xr(2)) + 2*this%cr(6)*(this%cr(3) + this%cr(5)*xr(1) + &
    !     2*this%cr(6)*xr(2)))*(this%cr(3) + this%cr(5)*xr(1) + 2*this%cr(6)*xr(2)))*((this%cr(2) + &
    !     2*this%cr(4)*xr(1) + this%cr(5)*xr(2))**2 + (this%cr(3) + this%cr(5)*xr(1) + &
    !     2*this%cr(6)*xr(2))**2 + 1)**(-1.5_r8)

    curvature = (-2*(this%cr(5) + this%cr(7))*((this%cr(2) + 2*this%cr(5)*xr(1) + &
        this%cr(6)*xr(2))**2 + (this%cr(3) + this%cr(6)*xr(1) + 2*this%cr(7)*xr(2))**2 + 1) + &
        (2*this%cr(5)*(this%cr(2) + 2*this%cr(5)*xr(1) + this%cr(6)*xr(2)) + &
        this%cr(6)*(this%cr(3) + this%cr(6)*xr(1) + 2*this%cr(7)*xr(2)))*(this%cr(2) + &
        2*this%cr(5)*xr(1) + this%cr(6)*xr(2)) + (this%cr(6)*(this%cr(2) + 2*this%cr(5)*xr(1) + &
        this%cr(6)*xr(2)) + 2*this%cr(7)*(this%cr(3) + this%cr(6)*xr(1) + &
        2*this%cr(7)*xr(2)))*(this%cr(3) + this%cr(6)*xr(1) + 2*this%cr(7)*xr(2)))/((this%cr(2) + &
        2*this%cr(5)*xr(1) + this%cr(6)*xr(2))**2 + (this%cr(3) + this%cr(6)*xr(1) + &
        2*this%cr(7)*xr(2))**2 + 1)**1.5_r8

    !print *, 'c2', 2*this%cr(7)

    ! print *, curvature

    ! xr = 0
    ! curvature = 2 * (-2*(this%cr(4) + this%cr(6))*((this%cr(2) + 2*this%cr(4)*xr(1) + &
    !     this%cr(5)*xr(2))**2 + (this%cr(3) + this%cr(5)*xr(1) + 2*this%cr(6)*xr(2))**2 + 1) + &
    !     (2*this%cr(4)*(this%cr(2) + 2*this%cr(4)*xr(1) + this%cr(5)*xr(2)) + &
    !     this%cr(5)*(this%cr(3) + this%cr(5)*xr(1) + 2*this%cr(6)*xr(2)))*(this%cr(2) + &
    !     2*this%cr(4)*xr(1) + this%cr(5)*xr(2)) + (this%cr(5)*(this%cr(2) + 2*this%cr(4)*xr(1) + &
    !     this%cr(5)*xr(2)) + 2*this%cr(6)*(this%cr(3) + this%cr(5)*xr(1) + &
    !     2*this%cr(6)*xr(2)))*(this%cr(3) + this%cr(5)*xr(1) + 2*this%cr(6)*xr(2)))*((this%cr(2) + &
    !     2*this%cr(4)*xr(1) + this%cr(5)*xr(2))**2 + (this%cr(3) + this%cr(5)*xr(1) + &
    !     2*this%cr(6)*xr(2))**2 + 1)**(-1.5_r8)
    ! print *, curvature

    curvature = clip(curvature, 1e10_r8, 0.0_r8)

  end function curvature

  function normal (this, x)

    use array_utils, only: normalize

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: normal(3)

    real(r8) :: xr(3)

    ASSERT(size(x)==3)

    xr = this%localCoords(x)

    normal(1) = -(this%cr(2) + 2*this%cr(5)*xr(1) + this%cr(6)*xr(2))
    normal(2) = -(this%cr(3) + this%cr(6)*xr(1) + 2*this%cr(7)*xr(2))
    normal(3) = 1.0_r8

    normal = matmul(transpose(this%rot), normalize(normal))

  end function normal

  function normal_average (this, interface_reconstruction) result(normal)

    use array_utils, only: normalize
    use polygon_type

    class(paraboloid), intent(in) :: this
    type(polygon), intent(in) :: interface_reconstruction
    real(r8) :: normal(3)

    type(polygon) :: int_rec
    real(r8) :: xb(2)

    ! average normal
    ! get bounds
    int_rec = interface_reconstruction
    call int_rec%rotate_offset(int_rec%norm, int_rec%centroid2())
    xb = [minval(int_rec%x(2,:)), maxval(int_rec%x(2,:))] ! WARN: 2D only

    normal(1) = 0
    normal(2) = sqrt((this%cr(3) + 2*this%cr(7)*xb(1))**2 + 1)/(2*this%cr(7)) &
        -       sqrt((this%cr(3) + 2*this%cr(7)*xb(2))**2 + 1)/(2*this%cr(7))
    normal(3) = -asinh(2*(this%cr(3)/(2*this%cr(7)) + xb(1))*abs(this%cr(7)))/(2*abs(this%cr(7))) &
        +        asinh(2*(this%cr(3)/(2*this%cr(7)) + xb(2))*abs(this%cr(7)))/(2*abs(this%cr(7)))

    normal = matmul(transpose(this%rot), normalize(normal))

  end function normal_average

  function localCoords(this, x) result(xr)
    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: xr(3)
    xr = matmul(this%rot, x - this%offset)
  end function localCoords

  function globalCoords(this, xr) result(x)
    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: xr(:)
    real(r8) :: x(3)
    x = matmul(transpose(this%rot), xr) + this%offset
  end function globalCoords

  ! calculate the curvature at a point x
  real(r8) function curvatureQdrtc (this,x)

    use consts, only: alittle
    use array_utils, only: clip

    class(paraboloid), intent(in) :: this
    real(r8), intent(in) :: x(:)

    ASSERT(size(x)==3)

    curvatureQdrtc = &
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
        2*this%coeff(10)*x(3)))

    curvatureQdrtc = clip(curvatureQdrtc, 1e10_r8, 0.0_r8)

  end function curvatureQdrtc

  function Fstr (this)

    use array_utils, only: isZero

    class(paraboloid), intent(in) :: this
    character(:), allocatable :: Fstr

    character(:), allocatable :: terms(:)
    character(32) :: term_str
    integer :: i

    terms = ['1   ','x   ','y   ','z   ','x**2','x*y ','x*z ','y**2','y*z ','z**2']
    !terms = ['1','x','y','x**2','x*y','y**2']
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

  function Fstr_rot (this) result(Fstr)

    use array_utils, only: isZero

    class(paraboloid), intent(in) :: this
    character(:), allocatable :: Fstr

    character(:), allocatable :: terms(:)
    character(32) :: term_str
    integer :: i

    terms = ['1   ','x   ','y   ','z   ','x**2','x*y ','y**2']
    Fstr = ''

    do i = 1,size(this%cr)
      if (.not.isZero(this%cr(i),1e-5_r8)) then
        if (this%cr(i) > 0.0_r8) then
          write(term_str,'(a,es10.3,2a)') '  + ',this%cr(i),'*',terms(i)
        else
          write(term_str,'(a,es10.3,2a)') '  - ',abs(this%cr(i)),'*',terms(i)
        end if
        Fstr = Fstr // trim(term_str)
      end if
    end do

  end function Fstr_Rot

  subroutine print_f (this, unit, iotype, vlist, iostat, iomsg)

    class(paraboloid), intent(in) :: this
    integer, intent(in) :: unit
    character(*), intent(in) :: iotype
    integer, intent(in) :: vlist(:)
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    write (unit, '(2a)', iostat=iostat) this%Fstr_rot(), ' = 0'

  end subroutine print_f

  subroutine print_uf (this, unit, iostat, iomsg)

    class(paraboloid), intent(in) :: this
    integer, intent(in) :: unit
    integer, intent(out) :: iostat
    character(*), intent(inout) :: iomsg

    write (unit, '(2a)', iostat=iostat) this%Fstr_rot(), ' = 0'

  end subroutine print_uf

end module paraboloid_type
