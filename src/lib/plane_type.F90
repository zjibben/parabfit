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
  use kinds,  only: r8
  use logging_services
  use polygon_type
  implicit none
  private

  ! dot(n,x) - rho = 0
  type, public :: plane
    real(r8) :: rho       ! interface plane constant
    real(r8) :: normal(3) ! interface normal
  contains
    procedure :: signed_distance
    procedure :: intersects
    procedure :: intersection_point
  end type plane

contains

  ! calculates the signed distance from a plane
  real(r8) function signed_distance (this,x)
    use consts, only: alpha

    class(plane), intent(in) :: this
    real(r8),     intent(in) :: x(3)

    signed_distance = sum(x*this%normal) - this%rho

    ! set distance to zero if it lies within alpha of the plane
    signed_distance = merge(signed_distance, 0.0_r8, abs(signed_distance) > alpha)
  end function signed_distance

  ! check if the plane lies between two points in space
  logical function intersects (this,x)
    use array_utils, only: isZero

    class(plane), intent(in) :: this
    real(r8),     intent(in) :: x(3,2) ! tuple of x positions

    real(r8) :: d1,d2

    d1 = this%signed_distance(x(:,1))
    d2 = this%signed_distance(x(:,2))
    
    ! if the signed distances have opposite signs, the two points are on opposite sides of the plane
    intersects = sign(1.0_r8,d1)/=sign(1.0_r8,d2) .or. isZero (d1) .or. isZero (d2)
    ! what does sign return for sign(1.0,0.0)?
  end function intersects

  ! return the point where the line between x1 & x2 intersect with the given plane
  function intersection_point (this,x)
    use array_utils, only: isZero

    class(plane), intent(in) :: this
    real(r8),     intent(in) :: x(3,2)
    real(r8)                 :: intersection_point(3)

    real(r8)                 :: dx(3),d1,d2

    if (.not.this%intersects (x)) call LS_fatal('edge does not intersect plane')

    d1 = this%signed_distance(x(:,1))
    d2 = this%signed_distance(x(:,2))
    
    if (isZero (d1)) then
      intersection_point = x(:,1)
    else if (isZero (d2)) then
      intersection_point = x(:,2)
    else
      dx = x(:,2)-x(:,1)
      intersection_point = x(:,1) - (d1/sum(dx*this%normal)) * dx
    end if

  end function intersection_point

end module plane_type
