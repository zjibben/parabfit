!!
!! ARRAY_UTILS
!!
!! This module contains a collection of utilities for arrays.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! April 2016
!!

#include "f90_assert.fpp"

module array_utils

  use kinds,  only: r8
  use consts, only: alittle
  use logging_services
  implicit none
  private

  public :: first_true_loc,last_true_loc,xrange,reorder,insertion_sort,int2str,&
      indexSort, &
      containsPair, containsValue, containsPoint, pointIndex, signs, &
      reverse,invert,isZero, index_of, clip, minmag, &
      meanArithmetic, meanHarmonic, &
      outer_product, &
      magnitude, magnitude2, normalize, normalizeIfNonzero, projectOnto, crossProduct, &
      orthonormalBasis, rotationMatrix, &
      eliminateNoise, interpolate, &
      det, diag, &
      polynomial_roots

  ! interface append
  !   procedure append_polygon
  ! end interface append

  interface interpolate
    module procedure interpolate_r8, interpolate_r8v
  end interface interpolate

  interface reorder
    module procedure reorder_r81d, reorder_r82d
  end interface reorder

  interface insertion_sort
    module procedure insertion_sort_3r8r8, insertion_sort_ir8, insertion_sort_r8
  end interface insertion_sort

  interface reverse
    module procedure reverse_i, reverse_r8r8
  end interface reverse

  interface isZero
    module procedure isZero_r8, isZero_r8a, isZero_r8aa
  end interface isZero

  interface meanArithmetic
    module procedure meanArithmeticNoWeights, meanArithmeticWeights
  end interface meanArithmetic

contains

  ! ! append element to the end of array (increasing array size by 1)
  ! subroutine append_polygon (array, element)
  !   use polygon_type

  !   class(polygon), intent(inout) :: array(:)
  !   class(polygon), intent(in)    :: element

  !   type(polygon)                 :: tmp(size(array))
  !   integer                       :: N

  !   N = size(array)
  !   tmp = array
  !   deallocate(array)
  !   allocate(array(N+1))
  !   array(1:N) = tmp
  !   array(N+1) = element

  ! end subroutine append_polygon

  ! return the index of the last true element of a logical array
  integer function last_true_loc (mask)

    logical, intent(in) :: mask(:)

    do last_true_loc = size(mask),1,-1
      if (mask(last_true_loc)) return
    end do
    last_true_loc = 0

  end function last_true_loc

  ! return the index of the first true element of a logical array
  integer function first_true_loc (mask)

    logical, intent(in) :: mask(:)

    do first_true_loc = 1,size(mask)
      if (mask(first_true_loc)) return
    end do
    first_true_loc = 0

  end function first_true_loc

  ! create an array [a,a+1,a+2,...,b-1,b]
  ! if b<a, then    [a,a-1,a-2,...,b+1,b]
  function xrange (a,b)
    integer, intent(in) :: a,b
    integer             :: xrange(abs(b-a)+1)

    integer             :: i

    xrange(1) = a
    do i = 2,abs(b-a)+1
      xrange(i) = xrange(i-1) + sign(1,b-a)
    end do

  end function xrange

  pure logical function containsPair (pair,list)

    integer, intent(in) :: pair(:), list(:,:)

    integer :: i

    containsPair = .false.
    do i = 1,size(list, dim=2)
      containsPair = containsPair .or. &
          all(list(:,i)==pair) .or. all(list(:,i)==reverse(pair))
    end do

  end function containsPair

  pure logical function containsPoint (x,list)

    real(r8), intent(in) :: x(:), list(:,:)

    integer :: i

    containsPoint = .false.
    do i = 1,size(list, dim=2)
      containsPoint = containsPoint .or. all(isZero(x - list(:,i)))
      if (containsPoint) return
    end do

  end function containsPoint

  pure integer function pointIndex(x,list)

    real(r8), intent(in) :: x(:), list(:,:)

    do pointIndex = 1,size(list,dim=2)
      if (all(isZero(x - list(:,pointIndex)))) return
    end do
    pointIndex = -1

  end function pointIndex

  pure logical function containsValue (n, array)
    integer, intent(in) :: n, array(:)
    containsValue = index_of(n, array) > 0
  end function containsValue

  ! return an array of indicies of a sorted from smallest a to largest a
  function indexSort (a)

    real(r8), intent(in) :: a(:)
    integer :: indexSort(size(a))
    real(r8) :: a_copy(size(a))

    a_copy = a
    indexSort = xrange(1,size(a))
    call insertion_sort (indexSort, a_copy)

  end function indexSort

  subroutine insertion_sort_3r8r8 (x,key)
    real(r8), intent(inout) :: x(:,:),key(:)

    real(r8) :: tmp,tmpX(size(x,dim=1))
    integer  :: i,j

    do i = 2,size(key)
      tmp  = key(i)
      tmpX = x(:,i)
      j = i
      do while (j>1)
        ! fortran doesn't guarantee short-circuiting, so we segfault if the two checks are together
        if (key(j-1) <= tmp) exit
        key(j) = key(j-1)
        x(:,j) = x(:,j-1)
        j = j-1
      end do
      key(j) = tmp
      x(:,j) = tmpX
    end do

  end subroutine insertion_sort_3r8r8

  subroutine insertion_sort_ir8 (x,key)
    integer,  intent(inout) :: x(:)
    real(r8), intent(inout) :: key(:)

    integer  :: tmpX
    real(r8) :: tmp
    integer  :: i,j

    do i = 2,size(key)
      tmp  = key(i)
      tmpX = x(i)
      j = i
      do while (j>1)
        ! fortran doesn't guarantee short-circuiting, so we segfault if the two checks are together
        if (key(j-1)>tmp) then
          key(j) = key(j-1)
          x(j)   = x(j-1)
          j = j-1
        else
          exit
        end if
      end do
      key(j) = tmp
      x(j)   = tmpX
    end do

  end subroutine insertion_sort_ir8

  subroutine insertion_sort_r8 (key)
    real(r8), intent(inout) :: key(:)

    real(r8) :: tmp
    integer  :: i,j

    do i = 2,size(key)
      tmp = key(i)
      j = i
      do while (j>1)
        ! fortran doesn't guarantee short-circuiting, so we segfault if the two checks are together
        if (key(j-1) <= tmp) exit
        key(j) = key(j-1)
        j = j-1
      end do
      key(j) = tmp
    end do

  end subroutine insertion_sort_r8

  ! reorder an array x based on a given order
  ! element i => element order(i)
  subroutine reorder_r81d (x,order)
    real(r8), intent(inout) :: x(:)
    integer,  intent(in)    :: order(:)

    real(r8)                :: y(size(x))
    integer                 :: i

    do i = 1,size(x)
      y(order(i)) = x(i)
    end do
    x = y
  end subroutine reorder_r81d

  subroutine reorder_r82d (x,order)
    real(r8), intent(inout) :: x(:,:)
    integer,  intent(in)    :: order(:)

    real(r8)                :: y(size(x,dim=1),size(x,dim=2))
    integer                 :: i

    do i = 1,size(x,dim=2)
      y(:,order(i)) = x(:,i)
    end do
    x = y
  end subroutine reorder_r82d

  ! these functions should be put in a different module for independent utilities
  character(len=20) function int2str (k)
    integer, intent(in) :: k
    write(int2str,*) k
    int2str = adjustl(int2str)
  end function int2str

  pure logical function isZero_r8 (x,tol)
    real(r8),           intent(in) :: x
    real(r8), optional, intent(in) :: tol

    real(r8) :: tolh

    tolh = 1e4_r8*alittle
    if (present(tol)) tolh = tol

    isZero_r8 = abs(x) < tolh
  end function isZero_r8

  pure function isZero_r8a (x,tol)
    real(r8),           intent(in) :: x(:)
    real(r8), optional, intent(in) :: tol
    logical                        :: isZero_r8a(size(x))

    real(r8) :: tolh

    tolh = 1e4_r8*alittle
    if (present(tol)) tolh = tol

    isZero_r8a = abs(x) < tolh
  end function isZero_r8a

  pure function isZero_r8aa (x,tol)
    real(r8),           intent(in) :: x(:,:)
    real(r8), optional, intent(in) :: tol
    logical                        :: isZero_r8aa(size(x,dim=1),size(x,dim=2))

    real(r8) :: tolh

    tolh = 1e4_r8*alittle
    if (present(tol)) tolh = tol

    isZero_r8aa = abs(x) < tolh
  end function isZero_r8aa

  pure real(r8) function clip (x,limit,default)

    use ieee_arithmetic, only: ieee_is_nan

    real(r8), intent(in) :: x
    real(r8), optional, intent(in) :: limit,default

    real(r8) :: limith

    if (present(limit)) then
      limith = limit
    else
      limith = huge(1.0_r8)
    end if

    clip = x
    if (abs(x) > limith .or. ieee_is_nan(x)) then
      if (present(default)) then
        clip = default
      else
        clip = limith
      end if
    end if

  end function clip

  pure logical function eq (a,b)
    real(r8), intent(in) :: a,b
    eq = isZero (a-b)
  end function eq

  ! return a reversed array
  pure function reverse_r8r8 (x) result(r)
    real(r8), intent(in) :: x(:,:)
    real(r8)             :: r(size(x,dim=1),size(x,dim=2))

    integer              :: i,N

    N = size(x,dim=2)
    do i = 1,N
      r(:,i) = x(:,N-(i-1))
    end do

  end function reverse_r8r8

  pure function reverse_i (x) result(r)
    integer, intent(in) :: x(:)
    integer             :: r(size(x))

    integer             :: i,N

    N = size(x)
    do i = 1,N
      r(i) = x(N-(i-1))
    end do

  end function reverse_i

  function invert (x)
    integer, intent(in) :: x(:)
    integer             :: invert(size(x))

    integer             :: i

    do i = 1,size(x)
      invert(x(i)) = i
    end do

  end function invert

  ! returns the index of the first entry of value n in array
  ! returns -1 if not found
  pure integer function index_of (n, array)
    integer, intent(in) :: n, array(:)

    do index_of = 1,size(array)
      if (array(index_of)==n) return
    end do
    index_of = -1 ! value not found
    ! call LS_fatal ('value not found in array')

  end function index_of

  real(r8) pure function meanArithmeticNoWeights (x)
    real(r8), intent(in) :: x(:)
    meanArithmeticNoWeights = sum(x) / real(size(x), r8)
  end function meanArithmeticNoWeights

  real(r8) pure function meanArithmeticWeights (x,w)
    real(r8), intent(in) :: x(:),w(:)
    meanArithmeticWeights = sum(x*w) / sum(w)
  end function meanArithmeticWeights

  real(r8) pure function meanHarmonic (x)
    real(r8), intent(in) :: x(:)
    meanHarmonic = 0.0_r8
    if (.not.any(x==0.0_r8)) meanHarmonic = 1.0_r8 / meanArithmetic(1.0_r8/x)
  end function meanHarmonic

  real(r8) pure function magnitude (v)
    real(r8), intent(in) :: v(:)
    magnitude = sqrt(magnitude2(v))
  end function magnitude

  real(r8) pure function magnitude2 (v)
    real(r8), intent(in) :: v(:)
    magnitude2 = sum(v**2)
  end function magnitude2

  pure function normalize (v)
    real(r8), intent(in) :: v(:)
    real(r8)             :: normalize(size(v))
    real(r8) :: mag
    mag = magnitude(v)
    normalize = v / merge(mag, huge(1.0_r8), mag > 0.0_r8)
  end function normalize

  pure subroutine normalizeIfNonzero (v)
    real(r8), intent(inout) :: v(:)
    real(r8) :: mag
    mag = magnitude(v)
    if (mag > alittle) v = v/mag
  end subroutine normalizeIfNonzero

  ! zeroes out components of an array that are below a certain threshold
  pure subroutine eliminateNoise (a, limit)
    real(r8), intent(inout) :: a(:)
    real(r8), intent(in), optional :: limit
    real(r8) :: limith
    limith = merge(limit, alittle, present(limit))
    where (abs(a) < limith) a = 0.0_r8
  end subroutine eliminateNoise

  ! project vector x1 into the direction of vector x2
  pure function projectOnto (x1, x2)
    real(r8), intent(in) :: x1(:), x2(:)
    real(r8)             :: projectOnto(size(x1))
    real(r8)             :: n2(size(x1))
    n2 = normalize(x2)
    projectOnto = dot_product(x1, n2) * n2
  end function projectOnto

  pure function crossProduct (a,b)
    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: crossProduct(size(a))
    crossProduct = [&
        a(2)*b(3) - a(3)*b(2), &
        a(3)*b(1) - a(1)*b(3), &
        a(1)*b(2) - a(2)*b(1)]
  end function crossProduct

  ! given a set of vectors x, return an orthonormal basis q for the same space using Gram-Schmidt
  function orthonormalBasis (x)

    use logging_services

    real(r8), intent(in) :: x(:,:)
    real(r8), allocatable :: orthonormalBasis(:,:)

    integer :: m, n, k, lwork, ierr
    real(r8), allocatable :: Q(:,:), tau(:), work(:)

    ! copy the input variables to a matrix to be modified in-place
    m = size(x, dim=1)
    n = size(x, dim=2)
    k = min(m,n)
    lwork = 8*n
    allocate(work(lwork), tau(k))
    Q = x

    ! compute the QR factorization
    call dgeqrf (m, n, Q, m, tau, work, lwork, ierr)
    if (ierr/=0) call ls_fatal ("failed dgeqrf in orthonormalBasis")

    call dorgqr (m, m, k, Q, m, tau, work, lwork, ierr)
    if (ierr/=0) call ls_fatal ("failed dorgqr in orthonormalBasis")

    ! return the first m columns of Q as the orthonormal basis with the same span as x
    orthonormalBasis = Q(:,1:m)

    ! real(r8), intent(in) :: x(:,:)
    ! real(r8), allocatable :: orthonormalBasis(:,:)

    ! integer :: i,j,n
    ! real(r8) :: tmp, v(size(x,dim=1)), vs(size(x,dim=1),size(x,dim=2))

    ! n = 1
    ! vs(:,1) = normalize(x(:,1))

    ! do i = 2,size(x,dim=2)
    !   v = x(:,i)
    !   do j = 1,n
    !     v = v - projectOnto(v,vs(:,j))
    !   end do

    !   tmp = magnitude(v)
    !   if (.not.isZero(tmp)) then
    !     n = n+1
    !     vs(:,n) = v / tmp
    !   end if
    ! end do

    ! orthonormalBasis = vs(:,1:n)

  end function orthonormalBasis

  pure function interpolate_r8 (f, xc, xf)
    real(r8), intent(in) :: f(:),xc(:,:),xf(:)
    real(r8) :: interpolate_r8
    real(r8) :: interpolation_factor
    ! ASSERT(size(f)==2)
    ! ASSERT(size(xc,2)==2)
    interpolation_factor = dot_product(xf-xc(:,1), xc(:,2)-xc(:,1))/magnitude2(xc(:,2)-xc(:,1))
    interpolate_r8 = (1.0_r8 - interpolation_factor)*f(1) + interpolation_factor*f(2)
  end function interpolate_r8

  pure function interpolate_r8v (f, xc, xf)
    real(r8), intent(in) :: f(:,:),xc(:,:),xf(:)
    real(r8) :: interpolate_r8v(size(f,1))
    real(r8) :: interpolation_factor
    ! ASSERT(size(f,2)==2)
    ! ASSERT(size(xc,2)==2)
    interpolation_factor = dot_product(xf-xc(:,1), xc(:,2)-xc(:,1))/magnitude2(xc(:,2)-xc(:,1))
    interpolate_r8v = (1.0_r8 - interpolation_factor)*f(:,1) + interpolation_factor*f(:,2)
  end function interpolate_r8v

  pure function outer_product (a,b)

    real(r8), intent(in) :: a(:), b(:)
    real(r8) :: outer_product(size(a),size(b))

    integer :: i

    do i = 1,size(b)
      outer_product(:,i) = a(:)*b(i)
    end do

  end function outer_product

  ! return the determinant of a matrix
  real(r8) function det (A)

    real(r8), intent(in) :: A(:,:)

    real(r8), allocatable :: LU(:,:)
    integer, allocatable :: piv(:)
    integer :: s, ierr

    s = size(A, dim=1)
    LU = A
    allocate(piv(s))
    call dgetrf (s,s, LU, s, piv, ierr)
    if (ierr>0) then
      det = 0.0_r8
      return
    else if (ierr<0) then
      print *, ierr
      print *, A
      call LS_fatal ('failed determinant calculation')
    end if

    ! determinant is the product of the diagonal of U in the LU decomposition
    det = 1.0_r8
    do s = 1,size(LU,dim=1)
      det = det * LU(s,s)
      if (piv(s) /= s) det = -det
    end do

  end function det

  function diag (diagonal_elements)

    real(r8), intent(in) :: diagonal_elements(:)
    real(r8) :: diag(size(diagonal_elements),size(diagonal_elements))

    integer :: i

    diag = 0.0_r8
    do i = 1,size(diagonal_elements)
      diag(i,i) = diagonal_elements(i)
    end do

  end function diag

  ! return the roots of the polynomial sum_{i=1}^N coeffs(i) * x^(i-1)
  ! only real roots are returned
  function polynomial_roots (coeffs)

    real(r8), intent(in) :: coeffs(:)
    real(r8), allocatable :: polynomial_roots(:)

    integer :: i, N, ierr
    real(r8), allocatable :: comp_mat(:,:), lr(:), li(:), vl(:,:), vr(:,:), work(:)

    N = size(coeffs)-1

    ! build the companion matrix
    allocate(comp_mat(N,N))
    comp_mat = 0.0_r8
    do i = 1,N-1
      comp_mat(i+1,i) = 1.0_r8
    end do
    comp_mat(:,N) = - coeffs(1:N) / coeffs(N+1)

    ! the roots are the eigenvalues of the companion matrix
    allocate(lr(N), li(N), work(5*N), vl(1,N), vr(1,N))
    call dgeev ('N','N', N, comp_mat, N, lr, li, vl, 1, vr, 1, work, 5*N, ierr)
    if (ierr/=0) call LS_fatal ('failed finding polynomial roots')

    polynomial_roots = pack(lr, mask=isZero(li))

  end function polynomial_roots

  function signs (a)

    real(r8), intent(in) :: a(:)
    real(r8) :: signs(size(a))

    integer :: i

    do i = 1,size(a)
      if (isZero(a(i))) then
        signs(i) = 0.0_r8
      else
        signs(i) = sign(1.0_r8, a(i))
      end if
    end do

  end function signs

  ! generate matrix which rotates such that normal is in the z direction
  function rotationMatrix (normal)

    real(r8), intent(in) :: normal(3)
    real(r8) :: rotationMatrix(3,3)

    real(r8) :: zdir(3)

    zdir = [0.0_r8, 0.0_r8, 1.0_r8]

    if (isZero(norm2(zdir - normal))) then
      rotationMatrix = 0
      rotationMatrix(1,1) = 1
      rotationMatrix(2,2) = 1
      rotationMatrix(3,3) = 1
    else
      rotationMatrix(3,:) = normal
      rotationMatrix(2,:) = normalize(crossProduct(zdir, normal))
      rotationMatrix(1,:) = crossProduct(rotationMatrix(2,:), normal)
    end if

  end function rotationMatrix

  ! return the argument with the minimum absolute value
  real(r8) pure function minmag(a,b)
    real(r8), intent(in) :: a, b
    minmag = a
    if (abs(a) > abs(b)) minmag = b
  end function minmag

end module array_utils
