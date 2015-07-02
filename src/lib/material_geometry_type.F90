!!
!! MATERIAL_GEOMETRY_TYPE
!!
!! This module defines a class that encapsulates an exact geometry for multiple materials,
!! particularly used for initialization of the Vof scalar.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

#include "f90_assert.fpp"

module material_geometry_type

  use kinds, only: r8
  use logging_services
  implicit none
  private
  
  ! define shape types
  type :: object_t
     integer :: matl_id
   contains
     procedure :: location_is_inside
  end type object_t

  type, extends(object_t) :: plane_t
     real(r8), dimension(:), allocatable :: normal
     real(r8) :: plane_const
  end type plane_t
  
  type, extends(object_t) :: sphere_t
     real(r8), dimension(:), allocatable :: center
     real(r8) :: radius
  end type sphere_t
  
  type :: object_ptr
     class(object_t), pointer :: o
   contains
     final :: delete_object_ptr
  end type object_ptr

  ! region type
  type, abstract, public :: base_region
   contains
     procedure(index_at), deferred :: index_at
  end type base_region
  
  abstract interface
     integer function index_at (this, x)
       use kinds, only: r8
       import base_region
       class(base_region), intent(in) :: this
       real(r8), intent(in) :: x(3)
     end function index_at
  end interface

  ! object_geometry type
  type, extends(base_region), public :: object_geometry
     private
     type(object_ptr), dimension(:), allocatable :: object
   contains
     !procedure :: init
     procedure :: index_at => object_index_at
     procedure :: object_at
  end type object_geometry
  
  ! material_geometry type
  type, extends(object_geometry), public :: material_geometry
   contains
     procedure :: init
     procedure :: index_at => material_at
  end type material_geometry
  
contains
  
  subroutine delete_object_ptr (this)
    type(object_ptr), intent(inout) :: this
    if (associated(this%o)) deallocate(this%o)
  end subroutine delete_object_ptr
  
  subroutine init (this, params, user_matl_id)
    use parameter_list_type
    
    class(material_geometry), intent(inout) :: this
    type(parameter_list), intent(in)    :: params
    integer, dimension(:), intent(in)   :: user_matl_id

    integer :: i,id,stat
    type(parameter_list_iterator)       :: param
    type(parameter_list), pointer       :: property
    character(:), allocatable           :: context,errmsg
    class(object_t), pointer            :: obj_ptr
    
    ! loop through every item in params
    ! each item describes a region of the domain as a shape
    ! and the material present in that region

    ! allocate the object array
    allocate(this%object(params%count()))
    
    ! generate an iterator
    param = parameter_list_iterator(params)

    i = 1
    do while(.not.param%at_end())
       context = 'processing ' // param%name() // ': '
       
       ! allocate the object as the user-given shape
       select case (param%name())
       case('sphere')
          allocate(sphere_t :: this%object(i)%o)
       case('plane')
          allocate(plane_t :: this%object(i)%o)
       case('all')
          ! 'all' should always be the last item in the list
          allocate(object_t :: this%object(i)%o)
       case default
          call LS_fatal (context//'unrecognized shape "'//param%name()//'"')
       end select
       
       ! read in details on the shape
       property => param%sublist()
       call property%get ('material-id', id, stat=stat, errmsg=errmsg)
       if (stat /= 0) call LS_fatal (context//errmsg)
       this%object(i)%o%matl_id = index_of(id, user_matl_id)
       
       select type (obj_ptr => this%object(i)%o)
       type is (sphere_t)
          call property%get ('center'     , obj_ptr%center , stat=stat, errmsg=errmsg)
          if (stat /= 0) call LS_fatal (context//errmsg)
          call property%get ('radius'     , obj_ptr%radius , stat=stat, errmsg=errmsg)
          if (stat /= 0) call LS_fatal (context//errmsg)
          ! should do some check here to ensure center is of length 3
          ASSERT(size(obj_ptr%center)==3)
       type is (plane_t)
          call property%get ('normal'     , obj_ptr%normal , stat=stat, errmsg=errmsg)
          if (stat /= 0) call LS_fatal (context//errmsg)
          call property%get ('plane-const'     , obj_ptr%plane_const , stat=stat, errmsg=errmsg)
          if (stat /= 0) call LS_fatal (context//errmsg)
          ! should do some check here to ensure center is of length 3
          ASSERT(size(obj_ptr%normal)==3)
          ! make sure the normal is normalized
          obj_ptr%normal = obj_ptr%normal / sum(obj_ptr%normal**2)
       type is (object_t)
          ! already done
       class default
          call LS_fatal (context//'unrecognized type')
       end select
       
       i = i+1
       call param%next()
    end do
    
  end subroutine init
  
  logical function location_is_inside(this, x)
    ! returns whether or not the given location is inside the described shape
    class(object_t), intent(in) :: this
    real(r8), intent(in) :: x(3)

    select type (this)
    type is (plane_t)
       location_is_inside = (sum(x*this%normal) - this%plane_const <= 0)
    type is (sphere_t)
       location_is_inside = (sum((x-this%center)**2) <= this%radius**2)
    class default ! all
       location_is_inside = .true.
       !call LS_fatal('location_is_inside: unrecognized shape')
    end select
    
  end function location_is_inside
  
  integer function material_at(this, x)
    ! determines and returns the id of the material at point x
    class(material_geometry), intent(in) :: this
    real(r8), dimension(3), intent(in) :: x
    class(object_t), pointer :: object

    object => this%object_at(x)
    material_at = object%matl_id
    
  end function material_at

  integer function object_index_at (this, x)
    ! identifies and returns the index of object at point x
    class(object_geometry), intent(in) :: this
    real(r8), dimension(3), intent(in) :: x

    integer :: i
    
    ! loop through every shape defined in the initial geometry
    do i = 1,size(this%object)
       ! check if we are inside the shape
       ! these shapes are ordered such that if materials overlap,
       ! the first one hit is the actual material at this location
       ! the last 'shape' covers the remainder of the domain, so
       ! there should always be a material discovered
       if (this%object(i)%o%location_is_inside (x)) then
          object_index_at = i
          return
       end if
    end do
    call LS_fatal('material_at: did not find any material at this location')
    
  end function object_index_at
  
  function object_at(this, x)
    ! identifies and returns the object at point x
    class(object_geometry), intent(in) :: this
    real(r8), dimension(3), intent(in) :: x
    class(object_t), pointer :: object_at

    integer :: i
    
    ! loop through every shape defined in the initial geometry
    do i = 1,size(this%object)
       ! check if we are inside the shape
       ! these shapes are ordered such that if materials overlap,
       ! the first one hit is the actual material at this location
       ! the last 'shape' covers the remainder of the domain, so
       ! there should always be a material discovered
       if (this%object(i)%o%location_is_inside (x)) then
          object_at => this%object(i)%o
          return
       end if
    end do
    call LS_fatal('material_at: did not find any material at this location')
    
  end function object_at

  integer function index_of(n, array)
    ! returns the index of the first entry of value n in array
    ! returns -1 if not found
    integer, intent(in) :: n
    integer, dimension(:), intent(in) :: array

    integer :: i

    do i = 1,size(array)
       if (array(i)==n) then
          index_of = i
          return
       end if
    end do

    ! value not found
    index_of = -1
    !call LS_fatal ('value not found in array')    
    
  end function index_of
  
end module material_geometry_type
