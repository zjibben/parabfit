!!
!! hex_types
!!
!! note: this could be derived from the polyhedron type
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

#include "f90_assert.fpp"

module hex_types

  use kinds, only: r8
  use consts, only: ndim, nfc, nvc
  use material_geometry_type
  use logging_services
  use plane_type
  implicit none
  private

  ! hex type to make divide and conquer algorithm simpler
  type, public :: base_hex
    real(r8) :: node(ndim,nvc), volume
  contains
    procedure :: calc_volume
    ! procedure         :: cell_center
    ! procedure         :: face_centers

    ! procedure         :: edge_centers
    ! procedure         :: contains_interface
  end type base_hex

  type, extends(base_hex), public :: cell_data
    real(r8) :: face_area(nfc), face_normal(ndim,nfc)
  contains
    procedure :: init => init_cell_data
    procedure :: calc_face_areas_and_normals
  end type cell_data

  type, extends(base_hex), public :: reconstruction_hex
    type(plane) :: P
    real(r8)    :: int_area  ! area of the interface for materials
    real(r8)    :: vof
    ! contains
    !   procedure              :: locate_plane
  end type reconstruction_hex

  integer, parameter, public :: hex_f(4,6) = reshape([ & ! face vertices
       3,4,8,7, & ! y+
       1,2,6,5, & ! y-
       1,5,8,4, & ! x-
       2,3,7,6, & ! x+
       1,4,3,2, & ! z-
       5,6,7,8],& ! z+
       shape(hex_f))
  ! [ & ! face vertices
  !      [ 3,4,8,7 ], & ! y+
  !      [ 1,2,6,5 ], & ! y-
  !      [ 1,5,8,4 ], & ! x-
  !      [ 2,3,7,6 ], & ! x+
  !      [ 1,4,3,2 ], & ! z-
  !      [ 5,6,7,8 ]  & ! z+
  !      ]

  integer, parameter, public :: hex_e(2,12) = reshape([ &
       1,2, & ! edge vertices
       2,3, &
       3,4, &
       4,1, &
       1,5, &
       2,6, &
       3,7, &
       4,8, &
       5,6, &
       6,7, &
       7,8, &
       8,5],&
       shape(hex_e))

  ! cube vertex positions for unit testing
  real(r8), parameter, public :: cube_v(3,8) = reshape([ &
       0.0_r8, 0.0_r8, 0.0_r8, &
       1.0_r8, 0.0_r8, 0.0_r8, &
       1.0_r8, 1.0_r8, 0.0_r8, &
       0.0_r8, 1.0_r8, 0.0_r8, &
       0.0_r8, 0.0_r8, 1.0_r8, &
       1.0_r8, 0.0_r8, 1.0_r8, &
       1.0_r8, 1.0_r8, 1.0_r8, &
       0.0_r8, 1.0_r8, 1.0_r8],&
       shape(cube_v))

  integer, parameter, public :: tet_fv(3,4) = reshape([&
      1,3,2,&
      1,2,4,&
      1,4,3,&
      2,3,4], shape(tet_fv))

  integer, parameter, public :: tet_ev(2,6) = reshape([&
      1,2,&
      1,3,&
      1,4,&
      2,3,&
      2,4,&
      3,4], shape(tet_ev))

  integer, parameter, public :: tet_fe(3,4) = reshape([&
        2,4,1,&
        1,5,3,&
        3,6,2,&
        4,6,5], shape(tet_fe))

  integer, parameter, public :: tet_ef(2,6) = reshape([&
      1,2,&
      1,3,&
      3,2,&
      1,4,&
      2,4,&
      3,4], shape(tet_ef))

  integer, parameter, public :: tet_vf(3,4) = reshape([&
      1,2,3,&
      1,2,4,&
      1,3,4,&
      2,3,4], shape(tet_vf))

contains

  subroutine init_cell_data (this, node, volume, face_area, face_normal, cfpar)

    use consts, only: ndim,nvc,nfc
    use array_utils, only: normalize

    class(cell_data),   intent(out) :: this
    real(r8),           intent(in)  :: node(:,:)
    real(r8), optional, intent(in)  :: volume, face_area(:), face_normal(:,:)
    integer,  optional, intent(in)  :: cfpar

    integer :: f

    ASSERT(all(shape(node)==[ndim,nvc]))

    this%node = node

    if (present(volume)) then
      this%volume = volume
    else
      this%volume = this%calc_volume ()
    end if

    if (present(face_area) .and. present(face_normal) .and. present(cfpar)) then
      ASSERT(size(face_area)==nfc)
      ASSERT(all(shape(face_normal)==[ndim,nfc]))
      this%face_area = face_area
      do f = 1,nfc
        !this%face_normal(:,f) = face_normal(:,f) / sqrt(sum(face_normal(:,f)**2)) ! normalize the input
        this%face_normal(:,f) = normalize(face_normal(:,f))
        if (btest(cfpar,f)) this%face_normal(:,f) = -this%face_normal(:,f) ! ensure the normals are outward facing
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

    integer :: f

    do f = 1,6
      this%face_normal(:,f) = quad_face_normal(this%node(:,hex_f(:,f)))
      this%face_area(f) = vector_length(this%face_normal(:,f))
      this%face_normal(:,f) = this%face_normal(:,f) / sqrt(sum(this%face_normal(:,f)**2))

      ! ensure the normals are outward facing
      this%face_normal(:,f) = calculate_outward_normal (this%face_normal(:,f), sum(this%node, dim=2)/8.0_r8, &
           sum(this%node(:,hex_f(:,f)), dim=2)/4.0_r8)
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
      outward_normal = -normal
    end if

  end function calculate_outward_normal

end module hex_types
