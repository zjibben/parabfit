!!
!! MATERIAL_GEOMETRY_TYPE
!!
!! This module defines a class that encapsulates an exact geometry for multiple materials,
!! particularly used for initialization of the Vof scalar.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

!! TODO: Make a way to initialize this type without a parameter list type, which really
!!       relies on JSON input files. I want to set up a material geometry list in the code
!!       as a unit test.

#include "f90_assert.fpp"

module material_geometry_type

  use kinds, only: r8
  use region_class
  use scalar_func_class
  use logging_services
  implicit none
  private

  ! region type
  type, abstract, public :: base_region
  contains
    procedure(index_at), deferred :: index_at
    procedure(smallest_signed_distance), deferred :: signed_distance
  end type base_region

  abstract interface
    integer function index_at (this, x)
      use kinds, only: r8
      import base_region
      class(base_region), intent(in) :: this
      real(r8), intent(in) :: x(:)
    end function index_at

    real(r8) function smallest_signed_distance(this, x)
      import :: base_region, r8
      class(base_region), intent(in) :: this
      real(r8), intent(in) :: x(:)
    end function smallest_signed_distance
  end interface

  ! material_geometry type
  type, extends(base_region), public :: material_geometry
    private
    class(scalar_func), allocatable :: matl_index
    type(region_box), allocatable :: region(:)
  contains
    generic :: init => init_params, init_matl_index
    procedure, private :: init_params
    procedure, private :: init_matl_index
    procedure :: index_at => material_at
    procedure :: signed_distance
  end type material_geometry

contains

  subroutine init_params (this, params, user_matl_id)

    use parameter_list_type
    use region_factories
    use scalar_func_containers
    use scalar_func_factories, only: alloc_const_scalar_func, alloc_piecewise_scalar_func
    use array_utils, only: index_of

    class(material_geometry), intent(inout) :: this
    type(parameter_list),     intent(in)    :: params
    integer,                  intent(in)    :: user_matl_id(:)

    integer :: i,nregions,matl_id,matl_index,stat
    type(parameter_list_iterator)       :: param_iter
    type(parameter_list), pointer       :: mparams,rparams
    character(:), allocatable           :: context,errmsg
    type(scalar_func_box), allocatable :: subfunc(:)

    ! loop through every item in params
    ! each item describes a region of the domain as a shape
    ! and the material present in that region

    ! allocate the object array
    nregions = params%count()
    allocate(subfunc(nregions), this%region(nregions))

    ! generate an iterator
    param_iter = parameter_list_iterator(params)

    i = 1
    do while(.not.param_iter%at_end())
      context = 'processing ' // param_iter%name() // ': '

      mparams => param_iter%sublist()

      ! get the user material id defined for this region
      call mparams%get ('material-id', matl_id, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)

      ! convert from user material id to vof index
      matl_index = index_of(matl_id, user_matl_id)

      ! place as constant in subfunction
      call alloc_const_scalar_func (subfunc(i)%f, real(matl_index, r8))

      ! get the user-defined region
      if (.not.mparams%is_sublist('region')) &
          call LS_fatal (context//'missing "region" sublist parameter')
      rparams => mparams%sublist('region')
      call alloc_region (this%region(i)%r, rparams)

      i = i+1
      call param_iter%next()
    end do

    ! build piecewise-constant function giving the material index over given regions
    call alloc_piecewise_scalar_func (this%matl_index, subfunc, this%region)

  end subroutine init_params

  subroutine init_matl_index (this, matl_index, rgn)

    class(material_geometry), intent(out) :: this
    class(scalar_func), intent(in) :: matl_index
    type(region_box), intent(in) :: rgn(:)

    ! this%matl_index = matl_index
    allocate(this%matl_index, source=matl_index)
    this%region = rgn

  end subroutine init_matl_index

  ! determines and returns the id of the material at point x
  integer function material_at(this, x)
    class(material_geometry), intent(in) :: this
    real(r8),                 intent(in) :: x(:)
    material_at = nint(this%matl_index%eval(x))
  end function material_at

  ! signed distance to interface
  ! currently only working for 2-phases
  real(r8) function signed_distance(this, x)
    class(material_geometry), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = this%region(1)%r%signed_distance(x)
  end function signed_distance

end module material_geometry_type
