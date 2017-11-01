!!
!! HALFSPHERE_REGION_TYPE
!!
!! A concrete implementation of the abstract base class REGION.
!! This implementation defines a region inside a given halfsphere.
!!
!! Zechariah Jibben <zjibben@lanl.gov>
!! April 2016
!!

module halfsphere_region_type

  use kinds, only: r8
  use region_class
  use sphere_region_type
  use plane_region_type
  implicit none
  private

  type, extends(region), public :: halfsphere_region
    private
    type(sphere_region), allocatable :: sphere
    type(plane_region),  allocatable :: plane
  contains
    procedure :: location_is_inside
    procedure :: signed_distance
  end type halfsphere_region

  interface halfsphere_region
    procedure halfsphere_region_value
  end interface halfsphere_region

contains

  !! constructor for HALFSPHERE_REGION objects
  function halfsphere_region_value (xc, radius, n) result(r)

    real(r8), intent(in) :: xc(:), radius, n(:)
    type(halfsphere_region) :: r

    allocate(r%sphere, source=sphere_region (xc, radius))
    allocate(r%plane, source=plane_region (n, dot_product(xc,n)))

  end function halfsphere_region_value

  logical function location_is_inside (this, x)

    class(halfsphere_region), intent(in) :: this
    real(r8),                 intent(in) :: x(:)

    location_is_inside = this%sphere%location_is_inside (x) .and. this%plane%location_is_inside (x)

  end function location_is_inside

  real(r8) function signed_distance(this, x)

    use array_utils, only: minmag

    class(halfsphere_region), intent(in) :: this
    real(r8), intent(in) :: x(:)

    signed_distance = this%plane%signed_distance(x)

    if (signed_distance < 0) &
        signed_distance = minmag(signed_distance, this%sphere%signed_distance(x))

  end function signed_distance

end module halfsphere_region_type
