!!
!! PRESCRIBED_VELOCITY_FIELDS
!!
!! This module contains a function for calculating various prescribed velocity fields
!! This could potentially be combined with velocity field initialization in all cases,
!! and with unstr_mesh_func
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

module prescribed_velocity_fields
  use kinds, only: r8
  use logging_services
  implicit none
  private

  public :: prescribed_velocity

contains

  function prescribed_velocity (x, t, field)
    real(r8), intent(in) :: x(3), t
    integer , intent(in) :: field
    real(r8)             :: prescribed_velocity(3)

    real(r8) :: periodT
    real(r8), parameter :: PI = 4.0_r8 * atan(1.0_r8)

    prescribed_velocity = 0.0_r8
    select case (field)
    case (1) ! deforming sphere
      periodT = 3.0_r8
      prescribed_velocity(1) = 2*sin(  PI*x(1))**2 * sin(2*PI*x(2))    * sin(2*PI*x(3))    * cos(PI*t/periodT)
      prescribed_velocity(2) = - sin(2*PI*x(1))    * sin(  PI*x(2))**2 * sin(2*PI*x(3))    * cos(PI*t/periodT)
      prescribed_velocity(3) = - sin(2*PI*x(1))    * sin(2*PI*x(2))    * sin(  PI*x(3))**2 * cos(PI*t/periodT)
    case (2) ! advecting plane in x-direction
      prescribed_velocity(1) = 0.5_r8 / 3.0_r8
    case (3) ! advecting plane in y-direction
      prescribed_velocity(2) = 0.5_r8 / 3.0_r8
    case (4) ! constant in x-y-z direction
      prescribed_velocity = 1.0_r8 / sqrt(3.0_r8)
    case default
      call LS_fatal ('unrecognized prescribed velocity field case')
    end select

  end function prescribed_velocity

end module prescribed_velocity_fields
