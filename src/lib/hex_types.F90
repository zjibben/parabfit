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

  public :: calculate_outward_normal ! this should really go somewhere else (maybe with the mesh?)
  
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

  ! truchas face ordering and node numbering
  ! integer, parameter, public :: face_node(4,6) = &
  !      [ &
  !      [4,8,7,3], & ! y-
  !      [5,1,2,6], & ! y+
  !      [5,8,4,1], & ! x+
  !      [6,2,3,7], & ! x-
  !      [3,2,1,4], & ! z-
  !      [7,8,5,6]  & ! z+
  !      ]

  ! ! truchas face ordering and pececillo node numbering
  ! integer, parameter, public :: face_node(4,6) = &
  !      [ &
  !      [2,6,5,1], & ! y-
  !      [7,3,4,8], & ! y+
  !      [7,6,2,3], & ! x+
  !      [8,4,1,5], & ! x-
  !      [1,4,3,2], & ! z-
  !      [5,6,7,8]  & ! z+
  !      ]

  ! pececillo face ordering and node numbering
  integer, parameter, public :: face_node(4,6) = &
       [ &
       [7,3,4,8], & ! y+
       [2,6,5,1], & ! y-
       [8,4,1,5], & ! x-
       [7,6,2,3], & ! x+
       [1,4,3,2], & ! z-
       [5,6,7,8]  & ! z+
       ]
  
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
        this%face_normal(:,f) = face_normal(:,f) / sqrt(sum(face_normal(:,f)**2)) ! normalize the input
      end do
    else
      call this%calc_face_areas_and_normals ()
    end if

    ! ensure the normals are outward facing
    do f = 1,6
      this%face_normal(:,f) = calculate_outward_normal (this%face_normal(:,f), sum(this%node, dim=2)/8.0_r8, &
           sum(this%node(:,face_node(:,f)), dim=2)/4.0_r8)
    end do
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

    integer :: f

    do f = 1,6
      this%face_normal(:,f) = quad_face_normal(this%node(:,face_node(:,f)))
      this%face_area(f) = vector_length(this%face_normal(:,f))
      this%face_normal(:,f) = this%face_normal(:,f) / sqrt(sum(this%face_normal(:,f)**2))
    end do
    
  end subroutine calc_face_areas_and_normals

  function calculate_outward_normal (normal, cell_center, face_center) result(outward_normal)
    real(r8), intent(in) :: normal(:), cell_center(:), face_center(:)
    real(r8)             :: outward_normal(3)
    
    real(r8) :: outward_dir(3)
    
    outward_dir = face_center - cell_center

    if ( sum(normal*outward_dir)>0.0_r8 ) then
      outward_normal = normal
    else
      outward_normal = - normal
    end if
    
  end function calculate_outward_normal
  
end module hex_types
