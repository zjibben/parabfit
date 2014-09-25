!!
!! BNDRY_FACE_FUNC_TYPE
!!
!! This module defines an implementation of the base class BNDRY_FUNC that
!! describes a function on the boundary faces of a mesh of type UNSTR_MESH.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for Fortran 2008, February 2014
!!
!! PROGRAMMING INTERFACE
!!
!! An instance of the derived type BNDRY_FACE_FUNC is intended to
!! encapsulate the data associated with a particular boundary condition.
!! Once initialized
!!
!!  INIT(MESH) initializes the object to begin receiving the specification of
!!    the face boundary data.  MESH is the UNSTR_MESH object over which the
!!    boundary data will be applied.  The object maintains an internal
!!    reference to the MESH, which
!!
!!  ADD(F, SETIDS, STAT, ERRMSG) assigns the allocatable SCALAR_FUNC class
!!    object F to compute the boundary data for the mesh faces
!!    with the mesh faces belonging to the face sets with IDs specified by
!!    the rank-1 array SETIDS.  The function F will be used to compute the
!!    boundary data, and its EVAL method should expect [t, x, y, z] as its
!!    vector argument.  This method must be called after INIT, and can be
!!    called multiple times, or even not at all.

#include "f90_assert.fpp"

module bndry_face_func_type

  use kinds, only: r8
  use bndry_func_class
  use unstr_mesh_type
  use scalar_func_class
  use scalar_func_containers
  implicit none
  private

  type, extends(bndry_func), public :: bndry_face_func
    private
    type(unstr_mesh), pointer :: mesh => null() ! reference only - do not own
    logical :: evaluated = .false.
    real(r8) :: tlast = -huge(1.0_r8)
    integer :: ngroup = -1
    integer, allocatable :: xgroup(:)
    type(scalar_func_box), allocatable :: farray(:)
    integer, allocatable :: hint(:)
    !! temporaries used during construction
    integer, allocatable :: tag(:)
    type(scalar_func_list) :: flist
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
  end type bndry_face_func

  !! Optimization hint values.
  integer, parameter :: BD_DATA_HINT_NONE    = 0
  integer, parameter :: BD_DATA_HINT_CONST   = 1
  integer, parameter :: BD_DATA_HINT_T_INDEP = 2
  integer, parameter :: BD_DATA_HINT_X_INDEP = 3

contains

  subroutine init (this, mesh)
    class(bndry_face_func), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    this%mesh => mesh
    this%ngroup = 0
    allocate(this%tag(size(mesh%face_set_mask)))
    this%tag = 0
  end subroutine init

  subroutine add (this, f, setids, stat, errmsg)
    class(bndry_face_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: f
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
    ASSERT(this%ngroup >= 0 .and. allocated(this%tag))
    call set_tag_array (this, setids, stat, errmsg)
    call this%flist%append (f)
  end subroutine add

  subroutine add_complete (this)

    use const_scalar_func_type

    class(bndry_face_func), intent(inout) :: this

    integer :: n, j

    !! Verify that THIS is in the correct state.
    ASSERT(this%ngroup >= 0 .and. allocated(this%tag))

    ASSERT(minval(this%tag) >= 0 .and. maxval(this%tag) <= this%ngroup)

    n = count(this%tag > 0)
    allocate(this%index(n), this%value(n), this%xgroup(this%ngroup+1))

    !! Prepare XGROUP: faces of group N will be FACES(XGROUP(N):XGROUP(N+1)-1).
    this%xgroup(1) = 1
    do n = 1, this%ngroup
      this%xgroup(n+1) = this%xgroup(n) + count(this%tag == n)
    end do

    !! Fill the FACES array; XGROUP(N) stores the next free location for condition N.
    do j = 1, size(this%tag)
      n = this%tag(j)
      if (n == 0) cycle
      this%index(this%xgroup(n)) = j
      this%xgroup(n) = 1 + this%xgroup(n)
    end do

    !! Restore XGROUP; XGROUP(N) is now the start of condition N+1 instead of N
    do n = this%ngroup, 1, -1
      this%xgroup(n+1) = this%xgroup(n)
    end do
    this%xgroup(1) = 1

    !! Convert the function list into the final function array.
    call scalar_func_list_to_box_array (this%flist, this%farray)

    deallocate(this%tag)

    !! For now we don't expose optimization hinting to the user,
    !! but we can determine directly which functions are constant.
    allocate(this%hint(this%ngroup))
    do n = 1, this%ngroup
      select type (f => this%farray(n)%f)
      type is (const_scalar_func)
        this%hint(n) = BD_DATA_HINT_CONST
      class default
        this%hint(n) = BD_DATA_HINT_NONE
      end select
    end do

  end subroutine add_complete


  subroutine compute (this, t)

    class(bndry_face_func), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: n, j
    real(r8) :: args(0:size(this%mesh%x,dim=1))

    !! Verify that THIS is in the correct state.
    ASSERT(allocated(this%index) .and. .not.allocated(this%tag))

    if (this%evaluated .and. t == this%tlast) return  ! values already set for this T

    args(0) = t
    do n = 1, this%ngroup
      associate(faces => this%index(this%xgroup(n):this%xgroup(n+1)-1), &
                value => this%value(this%xgroup(n):this%xgroup(n+1)-1))
        select case (this%hint(n))
        case (BD_DATA_HINT_CONST)
          if (.not.this%evaluated) value = this%farray(n)%f%eval(args)
        case (BD_DATA_HINT_X_INDEP)
          value = this%farray(n)%f%eval(args)
        case (BD_DATA_HINT_T_INDEP)
          if (.not.this%evaluated) then
            do j = 1, size(faces)
              args(1:) = sum(this%mesh%x(:,this%mesh%fnode(:,faces(j))),dim=2) / size(this%mesh%fnode,dim=1)
              value(j) = this%farray(n)%f%eval(args)
            end do
          end if
        case default
          do j = 1, size(faces)
            args(1:) = sum(this%mesh%x(:,this%mesh%fnode(:,faces(j))),dim=2) / size(this%mesh%fnode,dim=1)
            value(j) = this%farray(n)%f%eval(args)
          end do
        end select
      end associate
    end do

    this%tlast = t
    this%evaluated = .true.

  end subroutine compute

  !!
  !! This auxillary subroutine tags the faces identified by the list
  !! of face set IDs, checking for error conditions in the process.
  !!

  subroutine set_tag_array (this, setids, stat, errmsg)

    use string_utilities, only: i_to_c

    class(bndry_face_func), intent(inout) :: this
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: i, j
    logical :: mask(this%mesh%nface)
    integer(kind(this%mesh%face_set_mask)) :: bitmask

    !! Create the bitmask corresponding to SETIDS.
    bitmask = 0
    do i = 1, size(setids)
      do j = size(this%mesh%face_set_ID), 1, -1
        if (setids(i) == this%mesh%face_set_ID(j)) exit
      end do
      if (j == 0) then
        stat = 1
        errmsg = 'unknown face set ID: ' // i_to_c(setids(i))
        return
      end if
      bitmask = ibset(bitmask, j)
    end do

    !! Identify the faces specified by SETIDS.
    mask = (iand(bitmask, this%mesh%face_set_mask) /= 0)

    !! Check that these faces don't overlap those from preceding calls.
    if (any(mask .and. this%tag /= 0)) then
      stat = 1
      errmsg = 'face already associated with a function'
      return
    end if

    !! Verify that these faces are boundary faces.
    if (any(mask .and. .not.btest(this%mesh%face_set_mask,0))) then
      stat = 1
      errmsg = 'face not a boundary face'
      return
    end if

    !! Set the tag array.
    this%ngroup = 1 + this%ngroup
    !this%tag = merge(this%ngroup, this%tag, mask)
    where (mask) this%tag = this%ngroup

    stat = 0
    errmsg = ''

  end subroutine set_tag_array

end module bndry_face_func_type
