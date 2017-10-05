!!
!! MESH_GEOM
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! July 2015
!!
!! cneighbor(f,n) is the cell id of the neighbor to cell n opposite local face id f
!! caneighbor(1:nneighbor(n),n) is the set of all node neighbors to cell n
!! fneighbor(f,n) is the local face id of the face belonging to the neighbor to cell n
!!                opposite local face id f
!! fcell(n,f) is the cell id of local cell n connected to face f
!! flid(n,f) is the local face id of the global face f on the local cell n
!! vcell(:,v) are the cell ids of cells containing node v
!! xc(:,n) is the centroid of cell n
!! fc(:,f) is the centroid of face f
!! dx(:,f,i) is the vector from the centroid of its neighbor through face f to the centroid of cell i
!!           on a boundary, this is the vector from that face's centroid
!! dx_sclr(f,i) is the magnitude of dx(:,f,i)
!!

module mesh_geom_type

  use kinds,  only: r8
  use consts, only: ndim,nfc,nvf,nvc
  use unstr_mesh_type
  use set_type
  implicit none
  private

  type, public :: mesh_geom
    private
    integer,  allocatable, public :: cneighbor(:,:), fneighbor(:,:), vcell(:,:), fcell(:,:), &
        flid(:,:) !, caneighbor(:,:), nneighbor(:) !, nboundary_faces(:)
    real(r8), allocatable, public :: length(:), outnorm(:,:,:), xc(:,:), fc(:,:)
    type(set_integer), allocatable, public :: caneighbor(:)
  contains
    procedure :: init
  end type mesh_geom

contains

  subroutine init (this, mesh)

    use logging_services

    class(mesh_geom), intent(out) :: this
    type(unstr_mesh), intent(in)  :: mesh

    integer :: i,f

    ! set global variables
    select case (mesh%mesh_type)
    case ('TET')
      nfc = 4
      nvf = 3
      nvc = 4
    case ('HEX')
      nfc = 6
      nvf = 4
      nvc = 8
    case default
      call LS_fatal ('mesh_geom: unknown mesh type')
    end select

    allocate(this%cneighbor(nfc,mesh%ncell), this%fneighbor(nfc,mesh%ncell), &
        this%fcell(2,mesh%nface), this%flid(2,mesh%nface), this%length(mesh%nedge), &
        this%outnorm(ndim,nfc,mesh%ncell), this%xc(ndim,mesh%ncell), &
        this%fc(ndim,mesh%nface), this%caneighbor(mesh%ncell))

    call neighbor_ids (this%cneighbor, this%fneighbor, this%fcell, this%flid, mesh)
    this%vcell = cells_neighboring_vertices (mesh)
    call all_neighbor_ids (this%caneighbor, this%vcell, mesh)
    !this%length = edge_lengths (mesh%enode, mesh%nedge)

    !$omp parallel do default(private) shared(this,mesh,nfc,nvc)
    do i = 1,mesh%ncell
      ! calculate the outward normal
      do f = 1,nfc
        this%outnorm(:,f,i) = mesh%normal(:,mesh%cface(f,i)) / norm2(mesh%normal(:,mesh%cface(f,i)))

        ! check the orientation of the face with respect to this cell
        if (btest(mesh%cfpar(i),f)) this%outnorm(:,f,i) = -this%outnorm(:,f,i)
      end do

      ! calculate cell centroids
      ! note in truchas they do something far more complicated,
      ! probably taking the integral of x over the domain of the hex divided by the volume?
      this%xc(:,i) = sum(mesh%x(:,mesh%cnode(:,i)),dim=2) / nvc

      ! ! calculate the number of boundary faces touching this cell
      ! ! this is used to ensure values associated with the center of a cell touching a boundary
      ! ! are not updated until after all the boundary faces are updated
      ! this%nboundary_faces(i) = count(this%cneighbor(:,i) > 0)

    end do
    !$omp end parallel do

    ! calculate face centroids
    ! again, truchas does something different here
    !$omp parallel do
    do f = 1,mesh%nface
      this%fc(:,f) = sum(mesh%x(:,mesh%fnode(:,f)),dim=2) / nvf
    end do
    !$omp end parallel do

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
    integer, allocatable :: cells_neighboring_vertices(:,:)

    integer :: i,n,nid,j(mesh%nnode)

    ! get maximum number of cells attached to a node
    j = 0
    do i = 1,mesh%ncell
      do n = 1,nvc
        nid = mesh%cnode(n,i)
        j(nid) = j(nid) + 1
      end do
    end do

    allocate(cells_neighboring_vertices(maxval(j),mesh%nnode))
    cells_neighboring_vertices = -1

    j = 1
    do i = 1,mesh%ncell
      do n = 1,nvc
        nid = mesh%cnode(n,i)
        cells_neighboring_vertices(j(nid),nid) = i

        j(nid) = j(nid) + 1
      end do
    end do

  end function cells_neighboring_vertices

  ! generate an array that, given a cell and local face id,
  ! returns the neighboring cell and its local face id
  subroutine neighbor_ids (cneighbor, fneighbor, fcell, flid, mesh)
    use unstr_mesh_type

    integer,          intent(inout) :: cneighbor(:,:), fneighbor(:,:), fcell(:,:), flid(:,:)
    type(unstr_mesh), intent(in)    :: mesh

    integer :: cell,i,f,il,fid

    call cells_connected_to_faces(fcell,flid,mesh)

    !$omp parallel do default(private) shared(fneighbor,cneighbor,mesh,fcell,flid,nfc)
    do i = 1,mesh%ncell
      do f = 1,nfc
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
    !$omp end parallel do

  end subroutine neighbor_ids

  ! get the cell ids and local face ids for each face
  subroutine cells_connected_to_faces (fcell,flid,mesh)
    use unstr_mesh_type

    type(unstr_mesh), intent(in)  :: mesh
    integer,          intent(out) :: fcell(2,mesh%nface),flid(2,mesh%nface)

    integer :: i,f,fid,j(mesh%nface)

    j = 1; fcell = -1; flid = -1

    ! at each face of each cell, add the cell and local face id to a list
    do i = 1,mesh%ncell
      do f = 1,nfc
        fid = mesh%cface(f,i)

        fcell(j(fid),fid) = i
        flid(j(fid),fid)  = f

        j(fid) = j(fid)+1
      end do
    end do

  end subroutine cells_connected_to_faces


  subroutine all_neighbor_ids (caneighbor, vcell, mesh)

    type(set_integer), intent(inout) :: caneighbor(:)
    integer, intent(in) :: vcell(:,:)
    type(unstr_mesh), intent(in) :: mesh

    integer :: i,v,vc

    do i = 1,mesh%ncell
      do v = 1,nvc
        vc = mesh%cnode(v,i)
        call caneighbor(i)%add(pack(vcell(:,vc), mask=vcell(:,vc)>0 .and. vcell(:,vc)/=i))
      end do
    end do

  end subroutine all_neighbor_ids

end module mesh_geom_type
