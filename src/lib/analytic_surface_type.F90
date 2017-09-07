!!
!! analytic_surface
!!
!! This module defines a class for an analytic quadric surface description.
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
    logical :: initialized
  contains
    procedure :: init
    procedure :: bestFit
    procedure :: bestParaboloidFit
    procedure :: bestOneSheetFit
    procedure :: canonicalForm
    procedure, private :: goodParaboloidInSpace
    procedure, private :: taubinCoeffs
    procedure, private :: l
    procedure, private :: Dl
    procedure :: curvature
    procedure :: Fstr
    procedure, private :: print_uf
    procedure, private :: print_f
#ifndef NAGFOR
    ! NAG Compiler does not support this as of v6.1
    generic :: write(unformatted) => print_uf
    generic :: write(formatted) => print_f
#endif
  end type analytic_surface

  real(r8), parameter :: eccentricity_max = 1e4_r8

contains

  ! perform Taubin's method (1991) to calculate the analytic surface which fits the points x(:,:)
  ! in a least-squares sense
  subroutine init (this, coeff)
    class(analytic_surface), intent(out) :: this
    real(r8), intent(in) :: coeff(:)
    this%coeff = coeff
    this%initialized = .true.
  end subroutine init

  subroutine bestFit (this, x)

    use array_utils, only: isZero

    class(analytic_surface), intent(out) :: this
    real(r8), intent(in) :: x(:,:)

    real(r8), allocatable :: lr(:), vr(:,:)

    ASSERT(size(x,1)==3) ! Exclusively 3D for now. Adding a 2D version would be trivial, though.

    ! calculate the minimum-error quadric
    call this%taubinCoeffs (lr,vr, x)

    ! pick the eigenvector corresponding to the smallest real eigenvalue
    call this%init (vr(:,minloc(lr, dim=1)))

  end subroutine bestFit

  ! Andrews & Sequin (2013) algorithm for parabolic coefficients
  subroutine bestParaboloidFit (this, x)

    use array_utils, only: indexSort, normalize

    class(analytic_surface), intent(out) :: this
    real(r8), intent(in) :: x(:,:)

    real(r8) :: xcen(3)
    real(r8), allocatable :: lr(:), vr(:,:), C0(:), C1(:)
    integer, allocatable :: lr_index_sort(:)
    integer :: i,j

    ASSERT(size(x,1)==3) ! Exclusively 3D for now. Adding a 2D version would be trivial, though.

    ! calculate the Taubin basis of solutions
    call this%taubinCoeffs (lr,vr, x)

    ! get the two sets of coefficients associated with the smallest errors (eigenvalues)
    lr_index_sort = indexSort (lr)
    ! C0 = normalize(vr(:,minloc(lr, dim=1)))
    ! lr(minloc(lr, dim=1)) = huge(1.0_r8) ! nix this one to easily get the second smallest
    ! C1 = normalize(vr(:,minloc(lr, dim=1)))

    xcen = sum(x(:,1:3), dim=2) / 3.0_r8
    do j = 2,size(lr)
      C1 = normalize(vr(:,lr_index_sort(j)))
      do i = 1,j-1
        C0 = normalize(vr(:,lr_index_sort(i)))
        call this%goodParaboloidInSpace (C0, C1, xcen)
        if (this%initialized) return
      end do
    end do

    call LS_fatal ('could not find good paraboloidic fit')

  end subroutine bestParaboloidFit

  subroutine goodParaboloidInSpace (this, C0, C1, xcen)

    use array_utils, only: normalize, det, polynomial_roots, isZero, signs

    class(analytic_surface), intent(out) :: this
    real(r8), intent(in) :: C0(:), C1(:), xcen(:)

    real(r8) :: A0(3,3), A1(3,3), Atmp(3,3), cubic_coeffs(4), k, dx(3)
    real(r8), allocatable :: lr(:), vr(:,:), t_roots(:), lin(:), sqr_signs(:)
    integer :: i

    ! switch to matrix form of the coefficients
    A0 = coeffs2matrix (C0)
    A1 = coeffs2matrix (C1)

    ! construct cubic polynomial, the roots of which indicate a det=0 matrix (i.e., a paraboloid)
    cubic_coeffs = 0.0_r8

    cubic_coeffs(1) = det(A0)

    Atmp(:,1) = A1(:,1); Atmp(:,2) = A0(:,2); Atmp(:,3) = A0(:,3)
    cubic_coeffs(2) = cubic_coeffs(2) + det(Atmp)
    Atmp(:,1) = A0(:,1); Atmp(:,2) = A1(:,2); Atmp(:,3) = A0(:,3)
    cubic_coeffs(2) = cubic_coeffs(2) + det(Atmp)
    Atmp(:,1) = A0(:,1); Atmp(:,2) = A0(:,2); Atmp(:,3) = A1(:,3)
    cubic_coeffs(2) = cubic_coeffs(2) + det(Atmp)

    Atmp(:,1) = A0(:,1); Atmp(:,2) = A1(:,2); Atmp(:,3) = A1(:,3)
    cubic_coeffs(3) = cubic_coeffs(3) + det(Atmp)
    Atmp(:,1) = A1(:,1); Atmp(:,2) = A0(:,2); Atmp(:,3) = A1(:,3)
    cubic_coeffs(3) = cubic_coeffs(3) + det(Atmp)
    Atmp(:,1) = A1(:,1); Atmp(:,2) = A1(:,2); Atmp(:,3) = A0(:,3)
    cubic_coeffs(3) = cubic_coeffs(3) + det(Atmp)

    cubic_coeffs(4) = det(A1)

    t_roots = polynomial_roots(cubic_coeffs)

    !call this%init (normalize(C0 + minval(t_roots) * C1))

    ! print *, A0
    ! print *, det(A0)

    ! print *, 'c0: ', c0

    ! print *
    ! print *, 'c1: ', c1

    ! print *
    ! print *, 'tc: ', cubic_coeffs

    ! print *
    ! print *, 'roots: ',t_roots

    ! print *
    do i = 1,size(t_roots)
      dx = xcen
      call this%init (normalize(C0 + t_roots(i) * C1))
      call this%canonicalForm (lr, lin, k, dx)
      !if (isZero(k)) return ! no one-sheet quadrics of this form

      ! lr = lr/k
      ! lin = lin/k
      ! sqr_signs = signs(pack(lr, mask=.not.isZero(lr)))
      ! print *, lr
      ! print *, lin
      ! print *
      ! if (count(.not.isZero(lr)) == 2 .and. sqr_signs(1) == sqr_signs(2)) return

      ! if not too eccentric, return
      ! print *, lr
      ! print *, lin
      ! print *, k
      !if (dx(1) > 1e-4_r8) return
      if (sqrt(sum(lr**2)) < eccentricity_max) return
      !if (maxval(abs(lr)) < eccentricity_max) return
    end do

    ! do i = 1,size(x,2)
    !   print '(a,3es20.10)', 'x: ',x(:,i)
    ! end do

    this%initialized = .false.

    !call LS_fatal ('could not find paraboloidic fit in this space')

    ! do while (.not.isZero(det(coeffs2matrix(this%coeff))) .and. i <= size(t_roots))
    !   i = i + 1
    !   if (i > size(t_roots)) then
    !     i = i-1
    !     print *, A0
    !     print *
    !     print *, A1
    !     print *
    !     print *, 'root: ', t_roots(i)
    !     print '(a,4es15.5)', 'c: ', cubic_coeffs
    !     print *, 'det : ', det(A0 + t_roots(i)*A1)
    !     print *, 'det : ', det(coeffs2matrix(this%coeff))
    !     print *, 'detf: ', cubic_coeffs(1) + t_roots(i) * cubic_coeffs(2) + &
    !         t_roots(i)**2 * cubic_coeffs(3) + t_roots(i)**3 * cubic_coeffs(4)
    !     call LS_fatal ('could not find paraboloid fit')
    !   end if
    !   call this%init (normalize(C0 + t_roots(i) * C1))
    ! end do

  end subroutine goodParaboloidInSpace

  subroutine bestOneSheetFit (this, x)

    use array_utils, only: signs, isZero, diag, det

    class(analytic_surface), intent(out) :: this
    real(r8), intent(in) :: x(:,:)

    real(r8) :: k
    real(r8), allocatable :: sqr_signs(:), sqr(:), lin(:)
    integer :: ierr, N

    ASSERT(size(x,1)==3) ! Exclusively 3D for now. Adding a 2D version would be trivial, though.

    ! get the best fit
    call this%bestFit (x)

    ! get the canonical form of the best fit, so we can recognize its shape
    call this%canonicalForm (sqr, lin, k)

    ! The number of hyperbolic sheets is the number of
    ! eigenvalues (squared terms) with the opposite sign as
    ! the offset in the canonical equation. If we have a
    ! two-sheet hyperboloid, then find the best paraboloid fit.
    ! If k=0, then we have a cone. If any eigenvalues are zero
    ! or all have the same sign, then the best fit is not
    ! a hyperboloid.
    sqr_signs = signs(sqr)

    if (count(sqr_signs /= sign(1.0_r8, k)) > 1 .and. .not.isZero(k) .and. &
        .not.(any(isZero(sqr)) .or. all(sqr_signs==sqr_signs(1)))) & !then
        call this%bestParaboloidFit (x)

      ! call this%canonicalForm (evl, k)
      ! !print '(4es15.5)', evl, k
      ! print '(3es15.5)', evl/k
      ! print *, isZero(evl/k), det(coeffs2matrix(this%coeff))

      !print '(a,3(es15.5,a))', '(',evl(1)/k, ')*x + (',evl(2)/k, ')*y + (',evl(3)/k, ')*z + 1'
    !end if

    ! print '(4es15.5)', evl, k
    ! print '(3es15.5)', evl/k
    ! print *, evl_signs
    ! print *

  end subroutine bestOneSheetFit

  subroutine canonicalForm(this, sqr, lin, k, dx)

    use array_utils, only: diag, isZero, first_true_loc

    class(analytic_surface), intent(in) :: this
    real(r8), allocatable, intent(out) :: sqr(:), lin(:)
    real(r8), intent(out) :: k
    real(r8), intent(inout), optional :: dx(:)

    real(r8) :: D(3,3), evli(3), evc(3,3), tmp(1,3), work(5*3), tr(3), A(3,3), b(3), norm_fac, &
        eccentricity
    integer :: ierr, N, i

    N = 3
    if (allocated(sqr)) deallocate(sqr)
    if (allocated(lin)) deallocate(lin)
    allocate(sqr(N), lin(N))

    ! convert to matrix form
    A = coeffs2matrix(this%coeff) ! nonlinear terms
    b = this%coeff(2:4)           ! linear terms
    k = this%coeff(1)             ! constant term

    ! rotate to cancel out cross-terms (diagonalize the nonlinear terms matrix)
    call dgeev ('N','V', N, A, N, sqr, evli, tmp, 1, evc, N, work, 5*N, ierr)
    if (ierr/=0 .or. .not.all(isZero(evli))) call LS_fatal ('failed eigenvalue solve')

    D = diag(sqr)
    b = matmul(transpose(evc), b)

    ! translate to cancel out as many linear terms as possible
    tr = 0.0_r8
    do i = 1,N
      if (.not.isZero(sqr(i))) tr(i) = 0.5_r8 * b(i) / sqr(i)
    end do

    ! use a component of tr corresponding to a 0-eigenvalue to cancel out the constant term
    if (any(isZero(sqr) .and. .not.isZero(b))) then
      i = first_true_loc(isZero(sqr) .and. .not.isZero(b))
      tr(i) = (dot_product(tr, matmul(D,tr)) - dot_product(b, tr) + k) / b(i)
    end if

    k = dot_product(tr, matmul(D, tr)) - dot_product(b, tr) + k
    b = b - 2.0_r8 * matmul(tr, D)

    ! normalize the terms by the largest linear component, or the largest nonlinear component
    if (any(.not.isZero(b))) then
      norm_fac = b(first_true_loc(.not.isZero(b)))
    else
      norm_fac = sqr(first_true_loc(.not.isZero(sqr)))
    end if
    sqr = sqr / norm_fac
    lin = b / norm_fac
    k = k / norm_fac

    if (present(dx)) then
      eccentricity = sqrt(sum(sqr**2))
      dx = matmul(transpose(evc), dx) / norm_fac
      dx(1) = 2.0_r8 * sqrt(abs(dx(i))/eccentricity) * abs(norm_fac)
    end if

  end subroutine canonicalForm

  function coeffs2matrix(c)
    real(r8), intent(in) :: c(:)
    real(r8) :: coeffs2matrix(3,3)
    coeffs2matrix = reshape([&
        c(5), c(6)/2.0_r8, c(7)/2.0_r8,&
        c(6)/2.0_r8, c(8), c(9)/2.0_r8,&
        c(7)/2.0_r8, c(9)/2.0_r8, c(10)], [3,3])
  end function coeffs2matrix

  subroutine taubinCoeffs (this, evl, evc, x)

    use array_utils, only: outer_product, isZero
    use, intrinsic :: iso_c_binding, only: c_new_line

    class(analytic_surface), intent(in) :: this
    real(r8), allocatable, intent(out) :: evl(:), evc(:,:)
    real(r8), intent(in) :: x(:,:)

    real(r8), allocatable :: M(:,:), N(:,:), lr(:), li(:), vr(:,:), tmpM(:,:), tmpV(:), tmpV2(:)
    real(r8) :: tmpR
    integer :: ierr,s,i,j

    ! Exclusively 3D for now. Adding a 2D-specific version would be trivial, though.
    ASSERT(size(x,1)==3)

    s = 10 ! number of terms for a quadratic function

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
      if (isZero(li(i)) .and. .not. isZero(tmpV(i))) then
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
        2*this%coeff(10)*x(3)))

    curvature = clip(curvature, 1e10_r8, 0.0_r8)

  end function curvature

  function Fstr (this)

    use array_utils, only: isZero

    class(analytic_surface), intent(in) :: this
    character(:), allocatable :: Fstr

    character(:), allocatable :: terms(:)
    character(32) :: term_str
    integer :: i

    terms = ['1   ','x   ','y   ','z   ','x**2','x*y ','x*z ','y**2','y*z ','z**2']
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
