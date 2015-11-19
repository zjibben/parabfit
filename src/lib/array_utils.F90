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

  public :: first_true_loc,last_true_loc,xrange,reorder,insertion_sort,mag,prj,int2str,&
       reverse,invert,isZero

  ! interface append
  !   procedure append_polygon
  ! end interface append

  interface reorder
    procedure :: reorder_r81d
    procedure :: reorder_r82d
  end interface reorder

  interface insertion_sort
    procedure :: insertion_sort_3r8r8
    procedure :: insertion_sort_ir8
  end interface insertion_sort

  interface reverse
    procedure :: reverse_i
    procedure :: reverse_r8r8
  end interface reverse

  interface isZero
    procedure :: isZero_r8
    procedure :: isZero_r8a
    procedure :: isZero_r8aa
  end interface isZero

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

  pure real(r8) function dot (a,b)
    real(r8), intent(in) :: a(:),b(:)
    dot = sum(a*b)
  end function dot

  pure real(r8) function mag (a)
    real(r8), intent(in) :: a(:)
    mag = sqrt(sum(a**2))
  end function mag

  pure function prj (x,v)
    real(r8), intent(in) :: x(:),v(:)
    real(r8)             :: prj(size(x))
    prj = sum(x*v) * v / sum(v**2)
  end function prj
  
  subroutine insertion_sort_3r8r8 (x,key)
    real(r8), intent(inout) :: x(:,:),key(:)

    real(r8) :: tmp,tmpX(size(x,dim=1))
    integer  :: i,j

    do i = 2,size(key)
      tmp  = key(i)
      tmpX = x(:,i)
      j = i
      ! do while (j>1 .and. key(j-1)>tmp)
      !   key(j) = key(j-1)
      !   x(:,j) = x(:,j-1)
      !   j = j-1
      ! end do
      do while (j>1) ! fortran debug mode doesn't short-circuit, so we segfault if the two checks are together
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
      ! do while (j>1 .and. key(j-1)>tmp)
      !   key(j) = key(j-1)
      !   x(j)   = x(j-1)
      !   j = j-1
      ! end do
      do while (j>1) ! fortran debug mode doesn't short-circuit, so we segfault if the two checks are together
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
  function reverse_r8r8 (x) result(r)
    real(r8), intent(in) :: x(:,:)
    real(r8)             :: r(size(x,dim=1),size(x,dim=2))

    integer              :: i,N

    N = size(x,dim=2)
    do i = 1,N
      r(:,i) = x(:,N-(i-1))
    end do

  end function reverse_r8r8

  function reverse_i (x) result(r)
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

end module array_utils
