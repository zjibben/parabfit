!!
!! CYLINDER_REGION_TYPE
!!
!! A concrete implementation of the abstract base class REGION.
!! This implementation defines a region inside a given cylinder.
!!
!! Zechariah Jibben <zjibben@lanl.gov>
!! October 2016
!!

#include "f90_assert.fpp"

module cylinder_region_type

  use kinds, only: r8
  use consts, only: ndim
  use region_class
  implicit none
  private

  type, extends(region), public :: cylinder_region
    private
    real(r8), allocatable :: center(:), axis(:)
    real(r8)              :: radius, halfheight
  contains
    procedure :: location_is_inside
  end type cylinder_region

  interface cylinder_region
    procedure cylinder_region_value
  end interface cylinder_region

contains

  !! constructor for CYLINDER_REGION objects
  function cylinder_region_value (xc, axis, radius, halfheight) result(r)

    use array_utils, only: normalize
    
    real(r8), intent(in) :: xc(:), axis(:), radius, halfheight
    type(cylinder_region) :: r

    ASSERT(size(xc)==ndim)
    ASSERT(size(axis)==ndim)

    r%center = xc
    r%axis = normalize(axis)
    r%radius = radius
    r%halfheight = halfheight
    
  end function cylinder_region_value

  logical function location_is_inside (this, x)

    use array_utils, only: magnitude2

    class(cylinder_region), intent(in) :: this
    real(r8), intent(in) :: x(:)

    real(r8) :: xt(ndim), d, r

    ASSERT(size(x)==ndim)

    ! get distance from cylinder origin both along and orthogonal to the axis
    xt = x - this%center
    d = dot_product(xt,this%axis)
    r = magnitude2(xt - d*this%axis)

    location_is_inside = r <= this%radius**2 .and. abs(d) <= this%halfheight

  end function location_is_inside

end module cylinder_region_type
