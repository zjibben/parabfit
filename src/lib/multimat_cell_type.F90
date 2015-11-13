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
    integer          :: m,nm,nmat_in_cell

    allocate(this%mat_poly(size(vof)))

    call remainder%init (this)

    !write(*,*) 'here1.1'

    nmat_in_cell = count(cutvof < vof(:) .and. vof(:) < 1.0_r8-cutvof)
    nm = 0
    write(*,*) 'nmat',nmat_in_cell
    write(*,*) 'norm',norm(:,1)
    write(*,*) 'norm',norm(:,2)

    do m = 1,size(vof)
      if (vof(m) < cutvof .or. 1.0_r8-cutvof < vof(m)) cycle
      nm = nm+1 ! update the counter of how many materials we've seen thus far
      
      ! reconstruct the plane from the remaining free space
      ! use the plane to generate the polyhedron for this material,
      ! and update the free-space polyhedron

      !write(*,'(i3,es14.4,3f10.4)') m,vof(m),norm(:,m)

      if (nm==nmat_in_cell) then
        ! if this is the final material in the cell,
        ! it gets the entire remainder of the polyhedron
        this%mat_poly(m) = remainder
      else
        ! if this is not the final material in the cell, split the cell
        tmp = remainder%split (locate_plane_nd (remainder, norm(:,m), vof(m)*this%vol))
        remainder = tmp(1)
        this%mat_poly(m) = tmp(2)
      end if

      ! TODO: save the last face as the interface reconstruction for dumping purposes
    end do

    !write(*,*) 'here1.3'
    
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
