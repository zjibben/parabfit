!!
!! SINUSOID_REGION_TYPE
!!
!! A concrete implementation of the abstract base class REGION.
!! This implementation defines a region inside a given sinusoid.
!!
!! Zechariah Jibben <zjibben@lanl.gov>
!! September 2017
!!

#include "f90_assert.fpp"

module sinusoid_region_type

  use kinds, only: r8
  use consts, only: ndim
  use region_class
  implicit none
  private

  type, extends(region), public :: sinusoid_region
    private
    real(r8), allocatable :: c(:)
  contains
    procedure :: location_is_inside
    procedure :: signed_distance
  end type sinusoid_region

  interface sinusoid_region
    procedure sinusoid_region_value
  end interface sinusoid_region

  real(r8), parameter :: pi = 4*atan(1.0_r8)

contains

  !! constructor for SINUSOID_REGION objects
  function sinusoid_region_value (coeffs) result(r)
    real(r8), intent(in) :: coeffs(:)
    type(sinusoid_region) :: r
    ASSERT(size(coeffs)==7)
    r%c = coeffs
  end function sinusoid_region_value

  pure logical function location_is_inside (this, x)
    class(sinusoid_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    location_is_inside = &
        x(3) <= this%c(1) + this%c(2)*cos(this%c(3)*pi*(x(1) - this%c(4))) + &
        this%c(5)*cos(this%c(6)*pi*(x(2)-this%c(7)))
  end function location_is_inside

  ! TODO: any real estimate for the signed distance?
  real(r8) function signed_distance(this, x)
    class(sinusoid_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = this%c(1) + this%c(2)*cos(this%c(3)*pi*(x(1) - this%c(4))) + &
        this%c(5)*cos(this%c(6)*pi*(x(2)-this%c(7))) - x(3)
  end function signed_distance

end module sinusoid_region_type
