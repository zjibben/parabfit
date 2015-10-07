!!
!! array_utils
!!
!! This module contains a collection of utilities for arrays.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!

module array_utils
  use kinds, only: r8
  implicit none
  private

  public :: last_true_loc,xrange,reorder,insertion_sort,mag,prj

  ! interface append
  !   procedure append_polygon
  ! end interface append

  interface reorder
    procedure :: reorder_r81d
    procedure :: reorder_r82d
  end interface reorder

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
  ! this function is used in a few different modules -- it may be a good idea to bring it into one module
  function last_true_loc (mask)
    logical, intent(in) :: mask(:)
    integer             :: last_true_loc
    
    integer :: i

    do i = size(mask),1,-1
      if (mask(i)) then
        last_true_loc = i
        return
      end if
    end do
    last_true_loc = 0
    
  end function last_true_loc

  ! create an array [a,a+1,a+2,...,b-1,b]
  function xrange (a,b)
    integer, intent(in) :: a,b
    integer             :: xrange(b+1-a)

    integer             :: i

    do i = 1,b+1-a
      xrange(i) = a + (i-1)
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
  
  subroutine insertion_sort (x,key)
    real(r8), intent(inout) :: x(:,:),key(:)

    real(r8) :: tmp,tmpX(size(x,dim=1))
    integer  :: i,j

    do i = 2,size(key)
      tmp  = key(i)
      tmpX = x(:,i)
      j = i
      do while (j>1 .and. key(j-1)>tmp)
        key(j) = key(j-1)
        x(:,j) = x(:,j-1)
        j = j-1
      end do
      key(j) = tmp
      x(:,j) = tmpX
    end do

  end subroutine insertion_sort


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

end module array_utils
