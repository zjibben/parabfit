!!
!! FILL_REGION_TYPE
!!
!! A concrete implementation of the abstract base class REGION.
!! This implementation defines a region which fills the entire domain
!! (or what is left of it after allocating previous regions).
!!
!! Zechariah Jibben <zjibben@lanl.gov>
!! April 2016
!!

module fill_region_type

  use kinds, only: r8
  use region_class
  implicit none
  private

  type, extends(region), public :: fill_region
  contains
    procedure :: location_is_inside
    procedure :: signed_distance
  end type fill_region

contains

  logical function location_is_inside (this, x)
    class(fill_region), intent(in) :: this
    real(r8),           intent(in) :: x(:)
    location_is_inside = .true.
  end function location_is_inside

  real(r8) function signed_distance(this, x)
    class(fill_region), intent(in) :: this
    real(r8), intent(in) :: x(:)
    signed_distance = 0
  end function signed_distance

end module fill_region_type
