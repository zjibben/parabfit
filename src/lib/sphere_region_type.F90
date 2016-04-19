!!
!! SPHERE_REGION_TYPE
!!
!! A concrete implementation of the abstract base class REGION.
!! This implementation defines a region inside a given sphere.
!!
!! Zechariah Jibben <zjibben@lanl.gov>
!! April 2016
!!

module sphere_region_type

  use kinds, only: r8
  use region_class
  implicit none
  private

  type, extends(region), public :: sphere_region
    private
    real(r8), allocatable :: center(:)
    real(r8)              :: radius
  contains
    procedure :: location_is_inside
  end type sphere_region

  interface sphere_region
    procedure sphere_region_value
  end interface sphere_region

contains

  !! constructor for SPHERE_REGION objects
  function sphere_region_value (xc, radius) result(r)
    
    real(r8), intent(in) :: xc(:), radius
    type(sphere_region) :: r

    r%center = xc
    r%radius = radius
    
  end function sphere_region_value

  logical function location_is_inside (this, x)

    use array_utils, only: magnitude2

    class(sphere_region), intent(in) :: this
    real(r8),             intent(in) :: x(:)

    location_is_inside = magnitude2(x-this%center) <= this%radius**2

  end function location_is_inside

end module sphere_region_type
