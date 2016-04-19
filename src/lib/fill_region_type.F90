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

  use region_class
  implicit none
  private

  type, extends(region), public :: fill_region
  contains
    procedure :: location_is_inside
  end type fill_region

contains

  logical function location_is_inside (this, x)
    use kinds, only: r8
    class(fill_region), intent(in) :: this
    real(r8),           intent(in) :: x(:)
    location_is_inside = .true.
  end function location_is_inside

end module fill_region_type
