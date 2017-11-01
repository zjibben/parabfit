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
    procedure :: signed_distance
  end type ellipsoid_region

  interface ellipsoid_region
    procedure ellipsoid_region_value
  end interface ellipsoid_region

contains

  !! constructor for ELLIPSOID_REGION objects
  function ellipsoid_region_value (xc, axes) result(r)

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
    location_is_inside = norm2((x - this%center) / this%axes) < 1
  end function location_is_inside

  ! the ellipsoid signed distance function is approximate
  ! given by Ivey and Moin "Accurate interface normal and curvature...", 2015
  ! good enough for second order initialization
  real(r8) function signed_distance(this, x)
    class(ellipsoid_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = &
        -sqrt(this%axes(1)*this%axes(2)*this%axes(3)) * (1 - norm2((x - this%center) / this%axes))
  end function signed_distance

end module ellipsoid_region_type
