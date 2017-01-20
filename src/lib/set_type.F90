!!
!! SET_TYPE
!!
!! This module defines a type for unordered sets of data.
!! It can be extended to other types, but currently is only
!! implemented for sets of integers.
!!
!! Zechariah J Jibben <zjibben@lanl.gov>
!! January 2017
!! 

module set_type

  implicit none
  private

  type, public :: set_integer
    private
    integer, public :: n_elements
    integer, allocatable, public :: elements(:)
  contains
    procedure :: add_integer_single
    procedure :: add_integer_array
    generic :: add => add_integer_single, add_integer_array
  end type set_integer

contains

  ! add a single integer to the set
  subroutine add_integer_single (this, i)

    class(set_integer), intent(inout) :: this
    integer, intent(in) :: i

    integer, allocatable :: tmp(:)

    if (allocated(this%elements)) then
      ! if this integer isn't already in the set, add it
      if (.not.any(this%elements==i)) then
        ! extend the array of elements
        tmp = this%elements
        deallocate(this%elements)
        this%n_elements = this%n_elements + 1
        allocate(this%elements(this%n_elements))
        this%elements(1:this%n_elements-1) = tmp

        ! add the new item
        this%elements(this%n_elements) = i
      end if
    else
      this%n_elements = 1
      allocate(this%elements(this%n_elements))
      this%elements(this%n_elements) = i
    end if
    
  end subroutine add_integer_single

  ! add integers from an array to the set
  ! TODO: This could be made more efficient by checking how many
  !       unique, new elements are in iv prior to adding them.
  !       This way, we only reallocate once.
  subroutine add_integer_array (this, iv)

    class(set_integer), intent(inout) :: this
    integer, intent(in) :: iv(:)

    integer :: i
    
    do i = 1,size(iv)
      call this%add (iv(i))
    end do

  end subroutine add_integer_array

end module set_type
