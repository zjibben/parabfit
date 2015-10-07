!!
!! plane_type
!!
!! This module defines a plane type, along with routines for
!! calculating intersection points and distance.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!

module plane_type
  use kinds, only: r8
  use logging_services
  use polygon_type
  implicit none
  private

  type, public :: plane
    real(r8) :: rho       ! interface plane constant
    real(r8) :: normal(3) ! interface normal
  contains
    procedure :: signed_distance
    procedure :: intersects
    procedure :: intersection_point
  end type plane

contains

  !
  ! calculates the signed distance from a plane
  !
  real(r8) function signed_distance (this,x)
    class(plane), intent(in) :: this
    real(r8),     intent(in) :: x(3)

    ! if the signed distances are opposite, the two points
    ! are on opposite sides of the plane
    signed_distance = sum(x*this%normal) - this%rho
  end function signed_distance

  !
  ! check if the plane lies between two points in space
  !
  logical function intersects (this,x)
    class(plane), intent(in) :: this
    real(r8),     intent(in) :: x(3,2) ! tuple of x positions

    intersects = sign(1.0_r8,this%signed_distance(x(:,1)))/=sign(1.0_r8,this%signed_distance(x(:,2)))
  end function intersects

  !
  ! return the point where the line between x1 & x2 intersect with the given plane
  !
  function intersection_point (this,x)
    class(plane), intent(in) :: this
    real(r8),     intent(in) :: x(3,2)
    real(r8)                 :: intersection_point(3)

    real(r8)                 :: dx(3)

    dx = x(:,2)-x(:,1)
    intersection_point = x(:,1) - this%signed_distance (x(:,1)) / sum(dx*this%normal) * dx
  end function intersection_point

end module plane_type
