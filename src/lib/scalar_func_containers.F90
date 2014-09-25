!!
!! SCALAR_FUNC_CONTAINERS
!!
!! This module defines several containers for working with collections of
!! SCALAR_FUNC objects whose dynamic types may differ.  The functionality
!! provided is limited to just that required by the needs of the immediate
!! application code; this does not aim to be anything general.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  A polymorphic SCALAR_FUNC array must have a single dynamic type, and thus
!!  cannot be used to hold a collection of SCALAR_FUNC objects whose dynamic
!!  types may differ.  This module provides two derived types for addressing
!!  this problem.
!!
!!  The type SCALAR_FUNC_BOX wraps a public allocatable SCALAR_FUNC component:
!!
!!    TYPE SCALAR_FUNC_BOX
!!      CLASS(SCALAR_FUNC), ALLOCATABLE :: F
!!    END TYPE
!!
!!  Use an array of this type to hold a collection of SCALAR_FUNC objects of
!!  differing dynamic type.
!!
!!  The type SCALAR_FUNC_LIST provides a dynamic linked-list of SCALAR_FUNC
!!  objects.  Its components are private and it has a single type bound
!!  procedure APPEND(F) which appends the passed SCALAR_FUNC object F to the
!!  list.  The passed F must be allocatable, and the procedure moves the
!!  allocation from F to an internal list object; F is returned unallocated.
!!
!!  The subroutine SCALAR_FUNC_LIST_TO_BOX_ARRAY (LIST, ARRAY) converts the
!!  SCALAR_FUNC_LIST object LIST into a rank-1 SCALAR_FUNC_BOX array object.
!!  ARRAY must be allocatable; it is allocated to the correct size by the
!!  subroutine.  LIST is returned empty, the allocation for its SCALAR_FUNC
!!  objects having been moved into the returned ARRAY.
!!

#include "f90_assert.fpp"

module scalar_func_containers

  use scalar_func_class
  implicit none
  private

  public :: scalar_func ! re-export (necessary?)
  public :: scalar_func_box, scalar_func_list ! types
  public :: scalar_func_list_to_box_array ! procedure

  type :: scalar_func_list
    private
    integer :: n = 0
    type(list_func), pointer :: first => null()
  contains
    procedure :: append
    final :: delete
  end type

  type :: list_func
    class(scalar_func), allocatable :: f
    type(list_func), pointer :: next => null()
  end type

  type :: scalar_func_box
    class(scalar_func), allocatable :: f
  end type scalar_func_box

contains

  !! Final subroutine for SCALAR_FUNC_LIST objects.
  subroutine delete (this)
    type(scalar_func_list), intent(inout) :: this
    type(list_func), pointer :: rest
    do while (associated(this%first))
      rest => this%first%next
      deallocate(this%first)
      this%first => rest
    end do
  end subroutine

  subroutine append (this, f)
    class(scalar_func_list) :: this
    class(scalar_func), allocatable :: f
    type(list_func), pointer :: last
    ASSERT(allocated(f))
    if (associated(this%first)) then
      last => this%first
      do while (associated(last%next))
        last => last%next
      end do
      allocate(last%next)
      last => last%next
    else
      allocate(this%first)
      last => this%first
    end if
    call move_alloc (f, last%f)
    this%n = this%n + 1
  end subroutine append

  subroutine scalar_func_list_to_box_array (list, array)
    type(scalar_func_list) :: list
    type(scalar_func_box), allocatable, intent(out) :: array(:)
    integer :: j
    type(list_func), pointer :: rest
    allocate(array(list%n))
    do j = 1, list%n
      ASSERT(associated(list%first))
      call move_alloc (list%first%f, array(j)%f)
      rest => list%first%next
      deallocate(list%first)
      list%first => rest
    end do
    list%n = 0
    ASSERT(.not.associated(list%first))
  end subroutine scalar_func_list_to_box_array

end module scalar_func_containers
