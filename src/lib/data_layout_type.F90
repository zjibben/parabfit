!!
!! DATA_LAYOUT_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!! Frequently the degrees of freedom of a system will consist of a collection
!! of separately identifiable variables (e.g., cell temperatures, cell
!! enthalpies, face temperatures, boundary face radiosities, etc.) that need
!! to be packed into a single contiguous array of unknowns.  This module
!! provides a simple means of managing the packed layout of such data within
!! a contiguous array and providing access to the data segments of the array.
!!
!! PROGRAMMING INTERFACE
!!
!! This module defines derived data type DATA_LAYOUT (with private components)
!! with the following methods as type-bound procedures.  An instance of this
!! type is defined by one or more calls to ALLOC culminating in a call to
!! ALLOC_COMPLETE.
!!
!! No defined assignment is provided; do not use instances of this type in an
!! assignment statement unless you really know what the default assignment is
!! doing.
!!
!!  ALLOC(SIZE) allocates a contiguous segment of length SIZE in the layout
!!    and returns the integer handle that is used to access this segment.
!!
!!  ALLOC_COMPLETE() finalizes the layout after all the desired calls to ALLOC
!!    have been made.  Once called, no further calls to ALLOC are permitted.
!!
!! Once the layout has been defined by calls to ALLOC and ALLOC_COMPLETE, the
!! following methods are available.
!!
!!  DATA_PTR(ARRAY, SEGID) returns a pointer to the segment of the rank-1 array
!!    argument ARRAY specified by the segment handle SEGID of the layout.  Note
!!    very carefully the semantics of this function call.  The pointer result
!!    becomes aliased to (a portion of) ARRAY, and ARRAY must either be a
!!    pointer itself or have the TARGET attribute (Important!)  The function
!!    reference may appear directly in an expression or as the actual argument
!!    for an intent-in dummy argument, but in order to pass the target of the
!!    result to an intent-out or intent-inout dummy argument, a pointer must be
!!    associated with the function result and the pointer passed instead.
!!    ARRAY is of either intrinsic real or double precision types, and the
!!    result has the same type.
!!
!!  SIZE() returns the total size of the layout.  The array argument to
!!    DATA_PTR is expected to have this size.
!!
!!  GLOBAL_INDEX(SEGID, INDEX) returns the index into the layout that
!!    corresponds to the index INDEX of the segment SEGID of the layout.
!!    If DATA => LAYOUT%DATA_PTR(ARRAY, SEGID) then DATA(INDEX) and
!!    ARRAY(LAYOUT%GLOBAL_INDEX(SEGID, INDEX)) reference the same storage
!!    location.
!!

#include "f90_assert.fpp"

module data_layout_type

  implicit none
  private

  type :: seg_desc
    integer :: lb = 0, ub = 0, size = 0
  end type

  type :: seg_node
    type(seg_desc) :: seg
    type(seg_node), pointer :: next => null()
  end type

  type, public :: data_layout
    private
    integer :: nseg = 0
    integer :: size = 0
    type(seg_desc), allocatable :: seg(:)
    type(seg_node), pointer :: list => null()
  contains
    procedure :: alloc_segment
    procedure :: alloc_complete
    procedure :: layout_size
    procedure :: global_index
    procedure, private :: data_ptr_sp
    procedure, private :: data_ptr_dp
    generic :: data_ptr => data_ptr_sp, data_ptr_dp
    procedure, private :: segment_view_sp
    procedure, private :: segment_view_dp
    generic :: segment_view => segment_view_sp, segment_view_dp
    procedure, private :: segment_copy_sp
    procedure, private :: segment_copy_dp
    generic :: segment_copy => segment_copy_sp, segment_copy_dp
    final :: dealloc
  end type data_layout

contains

  integer function alloc_segment (this, size) result (segid)

    class(data_layout), intent(inout) :: this
    integer, intent(in) :: size

    type(seg_node), pointer :: new

    ASSERT(size >= 0)
    ASSERT(.not.allocated(this%seg))

    this%nseg  = this%nseg + 1
    segid      = this%nseg

    !! Create a new segment descriptor.
    allocate(new)
    new%seg%size = size
    new%seg%lb = this%size + 1
    new%seg%ub = this%size + size
    this%size  = new%seg%ub

    !! Prepend it to the list.
    new%next  => this%list
    this%list => new

  end function alloc_segment


  subroutine alloc_complete (this)

    class(data_layout), intent(inout) :: this

    integer :: id
    type(seg_node), pointer :: next

    ASSERT(.not.allocated(this%seg))

    !! Convert the linked list into a static array for quick lookup.
    !! The segment IDs were assigned in sequence and we pop things
    !! off the list in LIFO fashion.  (Perhaps it would be safer to
    !! have actually stored the ID in the seg_desc rather than just
    !! 'know' what it is implicitly.)
    allocate(this%seg(this%nseg))
    do id = this%nseg, 1, -1
      ASSERT(associated(this%list))
      this%seg(id) = this%list%seg
      next => this%list%next
      deallocate(this%list)
      this%list => next
    end do
    ASSERT(.not.associated(this%list))

  end subroutine alloc_complete


  elemental subroutine dealloc (this)
    type(data_layout), intent(inout) :: this
    type(seg_node), pointer :: next
    do while (associated(this%list))
      next => this%list%next
      deallocate(this%list)
      this%list => next
    end do
  end subroutine dealloc


  pure integer function layout_size (this)
    class(data_layout), intent(in) :: this
    layout_size = this%size
  end function layout_size


  integer function global_index (this, segid, index)
    class(data_layout), intent(in) :: this
    integer, intent(in) :: segid, index
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(index > 0 .and. index <= this%seg(segid)%size)
    global_index = index + this%seg(segid)%lb - 1
  end function global_index


  subroutine segment_view_sp (this, array, segid, view)
    class(data_layout), intent(in) :: this
    real, intent(in), target :: array(:)
    integer, intent(in) :: segid
    real, pointer :: view(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(allocated(this%seg))
    view => array(this%seg(segid)%lb:this%seg(segid)%ub)
  end subroutine segment_view_sp


  subroutine segment_view_dp (this, array, segid, view)
    class(data_layout), intent(in) :: this
    double precision, intent(in), target :: array(:)
    integer, intent(in) :: segid
    double precision, pointer :: view(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(allocated(this%seg))
    view => array(this%seg(segid)%lb:this%seg(segid)%ub)
  end subroutine segment_view_dp


  subroutine segment_copy_sp (this, array, segid, copy)
    class(data_layout), intent(in) :: this
    real, intent(in) :: array(:)
    integer, intent(in) :: segid
    real, intent(out) :: copy(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(allocated(this%seg))
    ASSERT(size(copy) >= this%seg(segid)%size)
    copy(:this%seg(segid)%size) = array(this%seg(segid)%lb:this%seg(segid)%ub)
  end subroutine segment_copy_sp


  subroutine segment_copy_dp (this, array, segid, copy)
    class(data_layout), intent(in) :: this
    double precision, intent(in) :: array(:)
    integer, intent(in) :: segid
    double precision, intent(out) :: copy(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(allocated(this%seg))
    ASSERT(size(copy) >= this%seg(segid)%size)
    copy(:this%seg(segid)%size) = array(this%seg(segid)%lb:this%seg(segid)%ub)
  end subroutine segment_copy_dp


  function data_ptr_sp (this, array, segid) result (ptr)
    class(data_layout), intent(in) :: this
    real, intent(in), target :: array(:)
    integer, intent(in) :: segid
    real, pointer :: ptr(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(allocated(this%seg))
    ptr => array(this%seg(segid)%lb:this%seg(segid)%ub)
  end function data_ptr_sp


  function data_ptr_dp (this, array, segid) result (ptr)
    class(data_layout), intent(in) :: this
    double precision, target, intent(in) :: array(:)
    integer, intent(in) :: segid
    double precision, pointer :: ptr(:)
    ASSERT(size(array) == this%size)
    ASSERT(segid > 0 .and. segid <= this%nseg)
    ASSERT(allocated(this%seg))
    ptr => array(this%seg(segid)%lb:this%seg(segid)%ub)
  end function data_ptr_dp

end module data_layout_type
