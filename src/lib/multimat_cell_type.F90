!!
!! MULTIMAT_CELL_TYPE
!!
!! This module provides a cell type that
!! describes internal material geometries
!! as child arbitrary polyhedra
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! Last revised 4 Nov 2012.
!!

module multimat_cell_type
  use kinds, only: r8
  use polyhedron_type
  implicit none
  private
  
  ! a multimat_cell is a polyhedron itself, describing
  ! the cell geometry, and also contains an array
  ! of polyhedra each describing the geometry of a
  ! particular material
  type, extends(polyhedron), public :: multimat_cell
    integer                       :: nmat
    !integer,          allocatable :: mat_id(:)
    type(polyhedron), allocatable :: mat_poly(:)
  contains
    procedure :: partition
    procedure :: volumes_behind_plane
  end type multimat_cell
  
contains

  ! given a set of VoFs, normals, and an order,
  ! create child polyhedra for each material
  subroutine partition (this, vof, norm)
    use consts, only: cutvof
    !use plane_type
    use locate_plane_nd_module

    class(multimat_cell), intent(inout) :: this
    real(r8),             intent(in)    :: vof(:), norm(:,:)

    type(polyhedron) :: tmp(2),remainder
    integer          :: m

    allocate(this%mat_poly(size(vof)))

    call remainder%init (this)

    do m = 1,size(vof)
      if (vof(m) < cutvof) cycle

      ! reconstruct the plane from the remaining free space
      ! use the plane to generate the polyhedron for this material,
      ! and update the free-space polyhedron
      tmp = remainder%split (locate_plane_nd (remainder, norm(:,m), vof(m)*this%vol))
      remainder = tmp(1)
      this%mat_poly(m) = tmp(2)

      ! TODO: save the last face as the interface reconstruction for dumping purposes
    end do
    
  end subroutine partition

  ! given a plane, find the volumes of each
  ! material behind that plane (flux volumes)
  function volumes_behind_plane (this, P) result(vol)
    use plane_type

    class(multimat_cell), intent(in) :: this
    class(plane),         intent(in) :: P
    real(r8)                         :: vol(size(this%mat_poly))

    integer :: m

    do m = 1,size(this%mat_poly)
      vol(m) = this%mat_poly(m)%volume_behind_plane (P)
    end do
    
  end function volumes_behind_plane

end module multimat_cell_type
