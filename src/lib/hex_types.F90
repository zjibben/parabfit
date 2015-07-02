!!
!!
!!
!!
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

module hex_types
  use kinds, only: r8
  use material_geometry_type
  use logging_services
  implicit none
  private

  ! hex type to make divide and conquer algorithm simpler
  type, public :: base_hex
     real(r8) :: node(3,8), volume
   contains
     procedure         :: calc_volume
     ! procedure         :: cell_center
     ! procedure         :: face_centers
     ! procedure         :: edge_centers
     ! procedure         :: contains_interface
  end type base_hex

  type, extends(base_hex), public :: reconstruction_hex
     real(r8) :: rho       ! interface plane constant
     real(r8) :: normal(3) ! interface normal
     real(r8) :: int_area  ! area of the interface for materials
     real(r8) :: vof
   ! contains
   !   procedure              :: locate_plane
  end type reconstruction_hex
  
  ! integer, parameter, public :: face_node(4,6) = &
  !      (/ &
  !      (/4,8,7,3/), &
  !      (/5,1,2,6/), &
  !      (/5,8,4,1/), &
  !      (/6,2,3,7/), &
  !      (/3,2,1,4/), &
  !      (/7,8,5,6/)  &
  !      /)
  
contains

  real(r8) function calc_volume (this)
    ! calculates the volume of a hex
    use cell_geometry, only: eval_hex_volumes
    class(base_hex), intent(in) :: this
    real(r8) :: cvol_tmp(8)
    call eval_hex_volumes(this%node, calc_volume, cvol_tmp)
  end function calc_volume
  
end module hex_types
