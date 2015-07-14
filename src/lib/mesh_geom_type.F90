!!
!! MESH_GEOM
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! July 2014
!!
!!
!! cneighbor(f,n) is the cell id of the neighbor to cell n opposite local face id f
!! fneighbor(f,n) is the local face id of the face belonging to the neighbor to cell n opposite local face id f
!! vcell(:,v) are the cell ids of cells containing node v
!!

module mesh_geom_type
  use kinds, only: r8
  use unstr_mesh_type
  implicit none
  private

  type, public :: mesh_geom
    private
    integer,  allocatable, public :: cneighbor(:,:), fneighbor(:,:), vcell(:,:)
    real(r8), allocatable, public :: length(:)
  contains
    procedure :: init
  end type mesh_geom
  
contains

  subroutine init (this, mesh)
    class(mesh_geom), intent(out) :: this
    type(unstr_mesh), intent(in) :: mesh

    allocate(this%cneighbor(6,mesh%ncell), this%fneighbor(6,mesh%ncell), &
         this%vcell(8,mesh%nnode), this%length(mesh%nedge))
    
    call neighbor_ids (this%cneighbor, this%fneighbor, mesh)
    this%vcell = cells_neighboring_vertices (mesh)
    !this%length = edge_lengths (mesh%enode, mesh%nedge)
    
  end subroutine init

  ! function edge_lengths (enode, nedge)
  !   real(r8), intent(in)  :: enode(:,:)
  !   integer,  intent(in)  :: nedge
  !   real(r8)              :: edge_lengths(nedge)

  !   integer :: e

  !   do e = 1,nedge
  !     edge_lengths(e) = sqrt(sum( (enode(1,e)-enode(2,e))**2 ))
  !   end do

  ! end function edge_lengths
  
  function cells_neighboring_vertices (mesh)
    type(unstr_mesh), intent(in) :: mesh
    integer :: cells_neighboring_vertices(8,mesh%nnode)

    integer :: i,n,nid,j(mesh%nnode)

    j = 1
    cells_neighboring_vertices = -1 

    do i = 1,mesh%ncell ! loop through all cells
      do n = 1,8 ! loop through every node on that cell
        nid = mesh%cnode(n,i)
        cells_neighboring_vertices(j(nid),nid) = i
        j(nid) = j(nid) + 1
      end do
    end do

  end function cells_neighboring_vertices

  ! there should be a better way of doing this
  ! generate an array that, given a cell and local face id,
  ! returns the neighboring cell and its local face id
  subroutine neighbor_ids (cneighbor, fneighbor, mesh)
    use unstr_mesh_type
    
    type(unstr_mesh), intent(in) :: mesh
    integer, intent(inout) :: cneighbor(6,mesh%ncell), fneighbor(6,mesh%ncell)

    integer :: cell,face,i,f,il,fid, fcell(2,mesh%nface),flid(2,mesh%nface)

    call cells_connected_to_faces(fcell,flid,mesh)

    do i = 1,mesh%ncell
      do f = 1,6
        fid = mesh%cface(f,i)

        ! get neighbor cell id local to face
        if (fcell(1,fid)==i) then
          il = 2
        else
          il = 1
        end if

        fneighbor(f,i) = flid(il,fid)
        cneighbor(f,i) = fcell(il,fid)
      end do
    end do

  end subroutine neighbor_ids

  ! get the cell ids and local face ids for each face
  subroutine cells_connected_to_faces (fcell,flid,mesh)
    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: mesh
    integer, intent(out) :: fcell(2,mesh%nface),flid(2,mesh%nface)

    integer :: i,f,fid,j(mesh%nface)

    j = 1; fcell = -1; flid = -1
    do i = 1,mesh%ncell
      do f = 1,6
        fid = mesh%cface(f,i)

        fcell(j(fid),fid) = i
        flid(j(fid),fid)  = f

        j(fid) = j(fid)+1
      end do
    end do

  end subroutine cells_connected_to_faces
  
end module mesh_geom_type
