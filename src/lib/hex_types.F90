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
    procedure :: calc_volume
    ! procedure         :: cell_center
    ! procedure         :: face_centers

    ! procedure         :: edge_centers
    ! procedure         :: contains_interface
  end type base_hex

  type, extends(base_hex), public :: cell_data
    real(r8) :: face_area(6), face_normal(3,6)
  contains
    procedure :: init => init_cell_data
    procedure :: calc_face_areas_and_normals
  end type cell_data
  
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
  
  subroutine init_cell_data (this, node, volume, face_area, face_normal)
    class(cell_data), intent(out) :: this
    real(r8),         intent(in) :: node(3,8)
    real(r8),         intent(in), optional :: volume, face_area(6), face_normal(3,6)

    integer :: f
    
    this%node   = node

    if (present(volume)) then
      this%volume = volume
    else
      this%volume = this%calc_volume ()
    end if

    if (present(face_area) .and. present(face_normal)) then
      this%face_area = face_area
      do f = 1,6
        this%face_normal(:,f) = face_normal(:,f) / sqrt(sum(face_normal(:,f)**2))
      end do
    else
      call this%calc_face_areas_and_normals ()
    end if
    
  end subroutine init_cell_data

  ! calculates the volume of a hex
  real(r8) function calc_volume (this)
    use cell_geometry, only: eval_hex_volumes
    class(base_hex), intent(in) :: this
    real(r8) :: cvol_tmp(8)
    call eval_hex_volumes(this%node, calc_volume, cvol_tmp)
  end function calc_volume

  subroutine calc_face_areas_and_normals (this)
    use cell_geometry, only: quad_face_normal, vector_length
    
    class(cell_data), intent(inout) :: this

    integer :: f, v(4)

    do f = 1,6
      select case (f)
      case(1)
        v = [4,8,7,3] ! Left face
      case(2)
        v = [5,1,2,6] ! Right face
      case(3)
        v = [5,8,4,1] ! Front face
      case(4)
        v = [6,2,3,7] ! Back face
      case(5)
        v = [3,2,1,4] ! Bottom face
      case(6)
        v = [7,8,5,6] ! Top face
      end select

      this%face_normal(:,f) = quad_face_normal(this%node(:,v))
      this%face_area(f) = vector_length(this%face_normal(:,f))
      this%face_normal(:,f) = this%face_normal(:,f) / sqrt(sum(this%face_normal(:,f)**2))
    end do
  end subroutine calc_face_areas_and_normals
  
end module hex_types
