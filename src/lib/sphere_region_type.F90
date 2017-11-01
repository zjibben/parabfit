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
    procedure :: signed_distance
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
    class(sphere_region), intent(in) :: this
    real(r8),             intent(in) :: x(:)
    location_is_inside = norm2(x-this%center) < this%radius
  end function location_is_inside

  real(r8) function signed_distance(this, x)
    class(sphere_region), intent(in) :: this
    real(r8),             intent(in) :: x(:)
    signed_distance = norm2(x-this%center) - this%radius
  end function signed_distance

end module sphere_region_type
