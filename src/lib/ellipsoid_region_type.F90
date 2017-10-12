!!
!! ELLIPSOID_REGION_TYPE
!!
!! A concrete implementation of the abstract base class REGION.
!! This implementation defines a region inside a given ellipsoid.
!!
!! Zechariah Jibben <zjibben@lanl.gov>
!! September 2017
!!

#include "f90_assert.fpp"

module ellipsoid_region_type

  use kinds, only: r8
  use consts, only: ndim
  use region_class
  implicit none
  private

  type, extends(region), public :: ellipsoid_region
    private
    real(r8), allocatable :: center(:), axes(:)
  contains
    procedure :: location_is_inside
  end type ellipsoid_region

  interface ellipsoid_region
    procedure ellipsoid_region_value
  end interface ellipsoid_region

contains

  !! constructor for ELLIPSOID_REGION objects
  function ellipsoid_region_value (xc, axes) result(r)

    use array_utils, only: normalize

    real(r8), intent(in) :: xc(:), axes(:)
    type(ellipsoid_region) :: r

    ASSERT(size(xc)==ndim)
    ASSERT(size(axes)==ndim)

    r%center = xc
    r%axes = axes

  end function ellipsoid_region_value

  pure logical function location_is_inside (this, x)
    class(ellipsoid_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    location_is_inside = norm2((x - this%center) / this%axes) <= 1
  end function location_is_inside

end module ellipsoid_region_type
