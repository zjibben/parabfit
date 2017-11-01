!!
!! REGION_CLASS
!!
!! This module defines the abstract base class REGION that provides
!! an interface to a general geometric region.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! April 2016
!!

module region_class

  use kinds, only: r8
  implicit none
  private

  type, abstract, public :: region
  contains
    procedure(location_is_inside_region), deferred :: location_is_inside
    procedure(signed_distance_to_surface), deferred :: signed_distance
  end type region

  abstract interface
    logical function location_is_inside_region (this, x)
      import :: region, r8
      class(region), intent(in) :: this
      real(r8), intent(in) :: x(:)
    end function location_is_inside_region

    real(r8) function signed_distance_to_surface(this, x)
      import :: region, r8
      class(region), intent(in) :: this
      real(r8), intent(in) :: x(:)
    end function signed_distance_to_surface
  end interface

  type, public :: region_box
    class(region), allocatable :: r
  end type region_box

end module region_class
