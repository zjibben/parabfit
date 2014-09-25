!!
!! DISTRIBUTED_MESH
!!
!! This module provides a parallel data structure that enscapsulates the
!! description of a distributed mesh, and a few general procedures that operate
!! on instances of the structure.  The data structure is currently designed to
!! handle the requirements of mimetic discretizations over tetrahedral and
!! hexahedral meshes.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Last revised 10 Apr 2007.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!!  The derived type DIST_MESH has the following public components.  Depending
!!  on the method used to instantiate the mesh object, not all components will
!!  necessarily be defined.
!!
!!  o NNODE, NEDGE, NFACE, NCELL are the number of nodes, edges faces, and
!!    cells in the mesh.
!!
!!  o CNODE(:,:) is the cell node array; CNODE(k,j) is the node number of local
!!    node k of cell j.  CNODE is dimensioned (4,NCELL) for a tetrahedral mesh
!!    and (8,NCELL) for a hexahedral mesh.
!!
!!  o CFACE(:,:) is the cell face array; CFACE(k,j) is the face number of local
!!    face k of cell j.  CFACE is dimensioned (4,NCELL) for a tetrahedral mesh
!!    and (6,NCELL) for a hexahedral mesh.
!!
!!  o CFPAR(:) is the cell face parity bitmask; BTEST(CFPAR(j),k) returns true
!!    if the orientation of face CFACE(k,j) is inward with respect to cell j.
!!
!!  o CEDGE(:,:) is the cell edge array; CEDGE(k,j) is the edge number of local
!!    edge k of cell j.  CEDGE is dimensioned (6,NCELL) for a tetrahedral mesh
!!    and (12,NCELL) for a hexahedral mesh.
!!
!!  o CEPAR(:) is the cell edge parity bitmask; BTEST(CEPAR(j),k) returns true
!!    if the orientation of edge CEDGE(k,j) is opposite to the local orientation
!!    of that edge with respect to cell j.
!!
!!  o FNODE(:,:) is the face node array; FNODE(k,j) is the node number of local
!!    node k of face j.  FNODE is dimensioned (3,NFACE) for a tetrahedral mesh,
!!    which has triangular faces, and (4,NFACE) for a hexahedral mesh, which
!!    has quad faces.
!!
!!  o FEDGE(:,:) is the face edge array; FEDGE(k,j) is the edge number of local
!!    edge k of face j.  FEDGE is dimensioned (3,NFACE) for a tetrahedral mesh,
!!    which has triangular faces, and (4,NFACE) for a hexahedral mesh, which
!!    has quad faces.
!!
!!  o FEPAR(:) is the face edge parity bitmask; BTEST(FEPAR(j),k) returns true
!!    if the orientation of edge FEDGE(k,j) is opposite to the local orientation
!!    of that edge with respect to face j.
!!
!!  o ENODE(:,:) is the edge node array; ENODE(k,j) is the node number of local
!!    node k of edge j.  ENODE is dimensioned (2,NEDGE).
!!
!!  o XNODE(:) is the external node numbering array; XNODE(j) is the external
!!    number (as defined in the mesh file, e.g.) of local node j.
!!
!!  o XCELL(:) is the external cell numbering array; XCELL(j) is the external
!!    number (as defined in the mesh file, e.g.) of local cell j.
!!
!!  o CBLOCK(:) is the cell block ID array; CBLOCK(j) is the block ID assigned
!!    to cell j.  BLOCK_ID(:) is the list of unique IDs in no particular order.
!!    IDs must be positive integers.  These are legacy arrays tied to Cubit's
!!    element blocks; their use is deprecated in favor of the CELL_SET_* arrays.
!!
!!  o CELL_SET_ID(:) is the list of unique cell set IDs (positive integers), and
!!    CELL_SET_MASK(:) is the cell set bitmask array; BTEST(CELL_SET_MASK(j),k)
!!    returns true if cell j belongs to the cell set with ID CELL_SET_ID(k).
!!
!!  o FACE_SET_ID(:) is the list of unique face set IDs (positive integers), and
!!    FACE_SET_MASK(:) is the face set bitmask array; BTEST(FACE_SET_MASK(j),k)
!!    returns true if face j belongs to the face set with ID FACE_SET_ID(k).
!!    BTEST(FACE_SET_MASK(J),0) returns true if face j is a boundary face.
!!
!!  o X(:,:) is the array of node coordinates, dimensioned (3,NNODE).
!!
!!  o LENGTH(:) is the vector of edge lengths, dimensioned (NEDGE).
!!
!!  o AREA(:) is the vector of face areas, dimensioned (NFACE).
!!
!!  o VOLUME(:) is the vector of cell volumes, dimensioned (NCELL).  This may
!!    be a signed volume that depends on the intrinsic orientation of the cell
!!    given by the local node order of the cell.
!!
!!  o NORMAL(:,:) is the array of oriented face areas; NORMAL(:,j) is the
!!    oriented face area of face j.
!!
!!  o CORNER_VOLUME(:,:) is the array of cell corner volumes; CORNER_VOLUME(k,j)
!!    is the volume of the 'corner' of cell j adjacent to the local vertex k.
!!    This is only relevant to hexahedral meshes, where the 'corner' is the
!!    tetrahedron subtended by the vertex k and the 3 adjacent vertices of the
!!    cell.
!!
!!  o NNODE_ONP, NEDGE_ONP, NFACE_ONP are the number of nodes, edges, and
!!    faces that are owned by the underlying process (on-process).
!!
!!  o NODE_IP, EDGE_IP, FACE_IP, and CELL_IP are derived types that describe
!!    the partitioning of the node, edge, face and cell index sets, and they
!!    include information necessary to communication certain off-process
!!    node, edge, face and cell data between processes.
!!
!! UTILITY PROCEDURES
!!
!!  CALL GET_FACE_SET_IDS (MESH, FACES, SETIDS)
!!    TYPE(DIST_MESH), INTENT(IN) :: MESH
!!    INTEGER, INTENT(IN) :: FACES(:)
!!    INTEGER, POINTER :: SETIDS(:)
!!
!!    This procedure returns a list of face set IDs in the array SETIDS.
!!    A set ID is included in the list if and only if some face in the
!!    specified list of faces belongs to the set.  The procedure allocates
!!    the return array SETIDS.  This is a parallel procedure returning the
!!    same result on all processes; the list of (local) faces may be
!!    process-dependent, of course.
!!

#include "f90_assert.fpp"

module unstr_mesh_type

  use kinds, only: r8
  !use bitfield_type
  implicit none
  private

  type, public :: unstr_mesh
    character(:), allocatable :: mesh_type
    integer :: nnode=0, nedge=0, nface=0, ncell=0
    !! Primary indexing arrays which define the mesh topology.
    integer, allocatable :: cnode(:,:)  ! cell nodes
    !integer, allocatable :: cedge(:,:)  ! cell edges
    integer, allocatable :: cface(:,:)  ! cell faces

    integer, allocatable :: cfpar(:)    ! relative cell face orientation (bitfield)
    !integer, allocatable :: cepar(:)    ! relative cell edge orientation (bitfield)

    !! Secondary indexing arrays derivable from the primary indexing arrays.
    integer, allocatable :: fnode(:,:)  ! face nodes
    !integer, allocatable :: fedge(:,:)  ! face edges
    !integer, allocatable :: enode(:,:)  ! edge nodes

    integer, allocatable :: fcell(:,:)  ! face cells

    !integer, allocatable :: fepar(:)    ! relative face edge orientation (bitfield)

    !! Relationship to external numbering.
    integer, allocatable :: xnode(:)    ! external node number
    integer, allocatable :: xcell(:)    ! external cell number

    !! Cell block ID arrays.
    integer, allocatable :: block_id(:) ! user-assigned ID for each cell block.
    integer, allocatable :: cblock(:)   ! cell block index.

    !! Cell set arrays.
    !integer, allocatable :: cell_set_id(:)
    !integer, allocatable :: cell_set_mask(:)

    !! Face set arrays.
    integer, allocatable :: face_set_id(:)
    integer, allocatable :: face_set_mask(:)
    !type(bitfield), allocatable :: face_set_mask(:)

    real(r8), allocatable :: x(:,:)
    !real(r8), allocatable :: length(:)
    real(r8), allocatable :: area(:)
    real(r8), allocatable :: volume(:)

    real(r8), allocatable :: normal(:,:)
    real(r8), allocatable :: corner_volume(:,:)
  contains
    procedure :: get_face_set_IDs
    procedure :: dump
  end type unstr_mesh

contains

  subroutine get_face_set_IDs (this, faces, setids)

    class(unstr_mesh), intent(in) :: this
    integer, intent(in) :: faces(:)
    integer, allocatable, intent(out) :: setids(:)

    integer :: j, n
    integer(kind(this%face_set_mask)) :: bitmask
    !type(bitfield) :: bitmask

    bitmask = 0 !ZERO_BITFIELD
    do j = 1, size(faces)
      bitmask = ior(bitmask, this%face_set_mask(faces(j)))
    end do
    bitmask = ibclr(bitmask, pos=0) ! clear the boundary flag

    !! Create the list of involved side set IDS.
    n = 0 ! count first to allocate
    do j = 1, size(this%face_set_id)
      if (btest(bitmask,j)) n = n + 1
    end do
    allocate(setids(n))
    n = 0 ! now store the data
    do j = 1, size(this%face_set_id)
      if (btest(bitmask,j)) then
        n = n + 1
        setids(n) = this%face_set_id(j)
      end if
    end do

  end subroutine get_face_set_IDs

  subroutine dump (this)

    class(unstr_mesh), intent(in) :: this

    integer :: j, k

    write(*,'(a,i0)') 'nnode=', this%nnode
    write(*,'(a,i0)') 'nedge=', this%nedge
    write(*,'(a,i0)') 'nface=', this%nface
    write(*,'(a,i0)') 'ncell=', this%ncell

    if (allocated(this%x)) then
      write(*,'(/,a)') 'node coordinates:'
      do j = 1, this%nnode
        write(*,'("[",i0,"]",3es12.4)') j, this%x(:,j)
      end do
    end if

    if (allocated(this%cnode)) then
      write(*,'(/,a)') 'cell node connectivity:'
      do j = 1, this%ncell
        write(*,'("[",i0,"]",8(" ",i4,:))') j, this%cnode(:,j)
      end do
    end if

    if (allocated(this%cface)) then
      write(*,'(/,a)') 'cell face connectivity:'
      do j = 1, this%ncell
        write(*,'("[",i0,"]",sp,8(" ",i5,:))') j, &
          (merge(-this%cface(k,j), this%cface(k,j), btest(this%cfpar(j),pos=k)), &
            k=1, size(this%cface,1))
      end do
    end if

    if (allocated(this%fnode)) then
      write(*,'(/,a)') 'face node connectivity:'
      do j = 1, this%nface
        write(*,'("[",i0,"]",4(" ",i4,:))') j, this%fnode(:,j)
      end do
    end if

    if (allocated(this%block_id)) then
      write(*, '(/,a,(t16,16(" ",i0,:)))') 'cell block ids=', this%block_id
    end if

    if (allocated(this%cblock)) then
      write(*,'(/,a)') 'cell block membership:'
      write(*,'(("[",i0,"] ",i0))') (j, this%cblock(j), j = 1, this%ncell)
    end if

    if (allocated(this%xcell)) then
      write(*,'(/,a)') 'external cell number:'
      write(*,'(("[",i0,"] ",i0))') (j, this%xcell(j), j = 1, this%ncell)
    end if

    if (allocated(this%xnode)) then
      write(*,'(/,a)') 'external node number:'
      write(*,'(("[",i0,"] ",i0))') (j, this%xnode(j), j = 1, this%nnode)
    end if

    if (allocated(this%area) .and. allocated(this%normal)) then
      write(*,'(/,a)') 'face areas and directed areas:'
      do j = 1, this%nface
        write(*,'("[",i0,"]",es12.4," /",3es12.4)') j, this%area(j), this%normal(:,j)
      end do
    end if

    if (allocated(this%volume)) then
      if (allocated(this%corner_volume)) then
        write(*,'(/,a)') 'cell volumes and corner volumes:'
        do j = 1, this%ncell
          write(*,'("[",i0,"]",es12.4," /",8es12.4)') j, this%volume(j), this%corner_volume(:,j)
        end do
      else
        write(*,'(/,a)') 'cell volumes:'
        do j = 1, this%ncell
          write(*,'("[",i0,"]",es12.4)') j, this%volume(j)
        end do
      end if
    end if

  end subroutine dump

end module unstr_mesh_type
