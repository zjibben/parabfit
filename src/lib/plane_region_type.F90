!!
!! PLANE_REGION_TYPE
!!
!! A concrete implementation of the abstract base class REGION.
!! This implementation defines a region behind a given, infinitely extending, plane.
!!
!! Zechariah Jibben <zjibben@lanl.gov>
!! April 2016
!!

module plane_region_type

  use kinds, only: r8
  use region_class
  implicit none
  private

  type, extends(region), public :: plane_region
    private
    real(r8), allocatable :: normal(:)
    real(r8)              :: plane_const
  contains
    procedure :: location_is_inside
  end type plane_region

  interface plane_region
    procedure plane_region_value
  end interface plane_region

contains

  !! constructor for PLANE_REGION objects
  function plane_region_value (n, p) result(r)
    
    real(r8), intent(in) :: n(:), p
    type(plane_region) :: r

    r%normal = n
    r%plane_const = p
    
  end function plane_region_value

  logical function location_is_inside (this, x)

    class(plane_region), intent(in) :: this
    real(r8),            intent(in) :: x(:)

    location_is_inside = dot_product(x,this%normal) - this%plane_const <= 0.0_r8

  end function location_is_inside

end module plane_region_type
