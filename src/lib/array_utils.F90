!!
!! array_utils
!!
!! This module contains a collection of utilities for arrays.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!

module array_utils

  use kinds,  only: r8
  use consts, only: alittle
  implicit none
  private

  public :: first_true_loc,last_true_loc,xrange,reorder,insertion_sort,int2str,containsPair, &
      reverse,invert,isZero, index_of, &
      meanArithmetic, meanHarmonic, &
      magnitude, magnitude2, normalize, projectOnto

  ! interface append
  !   procedure append_polygon
  ! end interface append

  interface reorder
    module procedure reorder_r81d, reorder_r82d
  end interface reorder

  interface insertion_sort
    module procedure insertion_sort_3r8r8, insertion_sort_ir8
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
        if (key(j-1)>tmp) then
          key(j) = key(j-1)
          x(:,j) = x(:,j-1)
          j = j-1
        else
          exit
        end if
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

    tolh = merge(tol, 1e4_r8*alittle, present(tol))

    isZero_r8 = abs(x) < tolh
  end function isZero_r8

  pure function isZero_r8a (x,tol)
    real(r8),           intent(in) :: x(:)
    real(r8), optional, intent(in) :: tol
    logical                        :: isZero_r8a(size(x))

    real(r8) :: tolh

    tolh = merge(tol, 1e4_r8*alittle, present(tol))    

    isZero_r8a = abs(x) < tolh
  end function isZero_r8a

  pure function isZero_r8aa (x,tol)
    real(r8),           intent(in) :: x(:,:)
    real(r8), optional, intent(in) :: tol
    logical                        :: isZero_r8aa(size(x,dim=1),size(x,dim=2))

    real(r8) :: tolh

    tolh = merge(tol, 1e4_r8*alittle, present(tol))

    isZero_r8aa = abs(x) < tolh
  end function isZero_r8aa

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
  integer function index_of (n, array)
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
    meanHarmonic = merge(1.0_r8 / meanArithmetic(1.0_r8/x), 0.0_r8, mask=.not.any(x==0.0_r8))
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
    normalize = v / magnitude(v)
  end function normalize

  ! project vector x1 into the direction of vector x2 
  pure function projectOnto (x1, x2)
    real(r8), intent(in) :: x1(:), x2(:)
    real(r8)             :: projectOnto(size(x1))
    real(r8)             :: n2(size(x1))
    n2 = normalize(x2)
    projectOnto = dot_product(x1, n2) * n2
  end function projectOnto

end module array_utils
