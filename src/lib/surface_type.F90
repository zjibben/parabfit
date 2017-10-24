!!
!! surface_type
!!
!! This module defines a surface type, which is a collection of polygons,
!! along with routines for writing surfaces to file in various formats.
!!
!! WARNING: This type is currently not thread safe, and should not be used
!!          to store interface reconstructions in threaded environments.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!

module surface_type

  use kinds, only: r8
  use polygon_type
  implicit none
  private

  ! a surface is defined as a collection of polygons
  type, public :: surface
    private
    type(polygon_box), public, allocatable :: element(:)
    integer, public, allocatable :: cell_id(:)
  contains
    !procedure, private :: append_polyhedron
    procedure, private :: append_polygon
    procedure, private :: append_surf
    generic   :: append => append_polygon, append_surf !, append_polyhedron
    !procedure :: write_gmv
    procedure :: write_ply
    procedure :: purge
    procedure :: local_patch
    procedure :: local_centroid
  end type surface

contains

  ! append element to the end of array (increasing array size by 1)
  subroutine append_polygon (this, new_element, cell_id)
    class(surface), intent(inout) :: this
    class(polygon_box), intent(in)    :: new_element
    integer,        intent(in)    :: cell_id

    type(polygon_box), allocatable :: tmp(:)
    integer,       allocatable :: tmp2(:)
    integer                    :: N

    !if (new_element%nVerts < 3) return ! TODO: throw error here?
    if (new_element%n_elements < 1) return

    if (allocated(this%element)) then
      N = size(this%element)
      allocate(tmp(N+1), tmp2(N+1))

      tmp(:N) = this%element
      tmp(N+1) = new_element
      call move_alloc(tmp, this%element)

      tmp2(:N) = this%cell_id
      tmp2(N+1) = cell_id
      call move_alloc(tmp2, this%cell_id)
    else
      allocate(this%element(1), this%cell_id(1))
      this%element(1) = new_element
      this%cell_id(1) = cell_id
    end if

  end subroutine append_polygon

  ! ! append a polyhedron to the surface
  ! ! (may instead want to write a routine in the polyhedron type that returns a surface)
  ! subroutine append_polyhedron (this, poly)

  !   use polyhedron_type

  !   class(surface),    intent(inout) :: this
  !   class(polyhedron), intent(inout) :: poly

  !   integer       :: nV,N,f
  !   type(polygon) :: face
  !   type(polygon_box) :: tmp(size(this%element))

  !   if (allocated(this%element)) then
  !     N = size(this%element)
  !     tmp = this%element
  !     deallocate(this%element)
  !     allocate(this%element(N+1))
  !     this%element(1:N) = tmp
  !   else
  !     N = 0
  !     allocate(this%element(1))
  !   end if

  !   this%element(N+1)%n_elements = poly%nFaces
  !   allocate(this%element(N+1)%elements(poly%nFaces))
  !   do f = 1,poly%nFaces
  !     nV = count(poly%face_vid(:,f) /= 0) ! number of vertices on this face
  !     call face%init (poly%x(:,poly%face_vid(1:nV,f)), poly%face_normal(:,f))

  !     this%element(N+1)%elements(f) = face
  !   end do

  ! end subroutine append_polyhedron

  ! append elements of another surface to
  subroutine append_surf (this, surf)
    class(surface), intent(inout) :: this
    class(surface), intent(in)    :: surf

    type(polygon_box)                 :: tmp(size(this%element))
    integer                       :: N,Nsurf

    Nsurf = size(surf%element)
    if (Nsurf==0 .or. .not.allocated(surf%element)) return

    if (allocated(this%element)) then
      N = size(this%element)
      tmp = this%element
      deallocate(this%element)
      allocate(this%element(N+Nsurf))
      this%element(1:N) = tmp
      this%element(N+1:N+Nsurf) = surf%element
    else
      allocate(this%element(Nsurf))
      this%element = surf%element
    end if

  end subroutine append_surf

  subroutine purge (this)
    class(surface), intent(inout) :: this
    if (allocated(this%element)) deallocate(this%element)
    if (allocated(this%cell_id)) deallocate(this%cell_id)
  end subroutine purge

  ! subroutine write_gmv (this)
  !   use gmvwrite_fortran_binding
  !   use array_utils, only: xrange

  !   class(surface), intent(in) :: this

  !   character(16)              :: filename
  !   integer                    :: e,j,Nelements,Nverts
  !   integer,  allocatable      :: element_vid(:,:)
  !   real(r8), allocatable      :: x(:,:),fakedata(:)

  !   ! generate contiguous array of nodes, and array of vertex ids
  !   Nelements = size(this%element)
  !   Nverts    = sum(this%element(:)%nVerts)

  !   allocate(x(3,Nverts), element_vid(maxval(this%element(:)%nVerts),size(this%element)), &
  !        fakedata(Nelements))

  !   ! To plot the surface elements, gmv needs data there.
  !   ! Or, at least, paraview won't open the file without it. Paraview won't display
  !   ! surfaces anyways, though, so it could be paraview that is wrong.
  !   fakedata = 0.0_r8

  !   j = 1
  !   do e = 1,Nelements
  !     x(:,j:j+this%element(e)%nVerts-1) = this%element(e)%x
  !     element_vid(1:this%element(e)%nVerts,e) = xrange (j,j+this%element(e)%nVerts-1)
  !     j = j + this%element(e)%nVerts
  !   end do

  !   ! dump the gmv file
  !   write(filename,'(a,i4.4)') 'planes.gmv'

  !   call gmvwrite_openfile_ir_ascii (trim(filename)//C_NULL_CHAR, 4, 8)

  !   ! write plane nodes
  !   !call gmv_write_unstr_mesh (this%mesh)
  !   call gmvwrite_node_data (Nverts, x(1,:), x(2,:), x(3,:))
  !   call gmvwrite_cell_header (0) ! GMV requires the cell header, even if there are no cells
  !   call gmvwrite_nodeids (xrange (1,Nverts))

  !   ! write the plane descriptions
  !   call gmvwrite_surface_header (size(this%element))
  !   do e = 1,Nelements
  !     call gmvwrite_surface_data (this%element(e)%nVerts, element_vid(1:this%element(e)%nVerts,e))
  !   end do

  !   call gmvwrite_surfvars_header ()
  !   call gmvwrite_surfvars_name_data ('surf'//C_NULL_CHAR, fakedata)
  !   call gmvwrite_surfvars_endsvar ()

  !   call gmvwrite_closefile ()

  !   deallocate (x, element_vid, fakedata)

  ! end subroutine write_gmv

  ! write the polygons using the Stanford PLY format
  subroutine write_ply (this, fname)
    class(surface), intent(inout) :: this
    character(*),   intent(in)    :: fname

    integer                       :: e,j,i,Nelements,Nverts, k, fh

    if (.not.allocated(this%element)) return

    Nelements = sum(this%element(:)%n_elements)
    Nverts = 0
    do i = 1,size(this%element)
      Nverts = Nverts + sum(this%element(i)%elements(:)%nVerts)
    end do

    ! open file
    open(newunit=fh, file=trim(fname))

    ! write PLY header
    write(fh,'(a)')    'ply'
    write(fh,'(a)')    'format ascii 1.0'
    write(fh,'(a,i9)') 'element vertex ',Nverts
    write(fh,'(a)')    'property float32 x'
    write(fh,'(a)')    'property float32 y'
    write(fh,'(a)')    'property float32 z'
    write(fh,'(a,i9)') 'element face ',Nelements
    write(fh,'(a)')    'property list uint8 int32 vertex_index'
    write(fh,'(a)')    'end_header'

    ! write polygon data
    ! vertices
    do e = 1,size(this%element)
      do k = 1,this%element(e)%n_elements
        do j = 1,this%element(e)%elements(k)%nVerts
          write(fh,'(3f20.10)') this%element(e)%elements(k)%x(:,j)
        end do
      end do
    end do

    ! faces
    j = 0
    do e = 1,size(this%element)
      do k = 1,this%element(e)%n_elements
        write(fh,'(i9)',advance='no') this%element(e)%elements(k)%nVerts
        do i = j,j+this%element(e)%elements(k)%nVerts-1
          write(fh,'(i9)',advance='no') i
        end do
        write(fh,*)
        j = j+this%element(e)%elements(k)%nVerts
      end do
    end do

    ! clean up
    close (fh)

  end subroutine write_ply

  ! return the polygons immediately surrounding the given cell
  !
  ! note 1: This could never be executed, in the case we ask for curvature in a cell
  !         neighboring the interface, but it needs to be. I'm not sure how to deal
  !         with this yet, I'll have to look at how the height function method
  !         handles it. I assume you would want to look at the primary normal direction
  !         of the interface nearby, then grab the curvature in the cell directly
  !         "below" this one. That could catastrophically break down when the curvature
  !         is on the order of the mesh spacing, though. This might also add some
  !         complexity by requiring curvature computation in cells containing the
  !         interface first, then looking up the values in a table later.
  function local_patch (this, cell_id, gmesh, vof, verbose)

    use mesh_geom_type
    use set_type

    class(surface), intent(in) :: this
    integer, intent(in) :: cell_id
    type(mesh_geom), intent(in) :: gmesh
    real(r8), intent(in) :: vof(:)
    type(polygon_box), allocatable :: local_patch(:)
    logical, intent(in), optional :: verbose

    integer :: e, npolygons, i, j
    integer, allocatable :: neighbor(:), polygon_id(:)
    type(set_integer) :: ngbr
    logical :: verboseh

    ! get the neighboring cell ids
    !neighbor = pack(gmesh%cneighbor(:,cell_id), mask=gmesh%cneighbor(:,cell_id)>0) ! face neighbors
    neighbor = gmesh%caneighbor(cell_id)%elements ! node neighbors

    ! ! stretch out to neighbors of neighbors
    ! call ngbr%add (neighbor)
    ! do e = 1,size(neighbor)
    !   call ngbr%add (gmesh%caneighbor(neighbor(e))%elements)
    ! end do
    ! neighbor = ngbr%elements

    allocate(polygon_id(size(neighbor)+1))

    ! get polygons from these cells
    ! the order doesn't matter as long as the polygon associated with cell_id is first
    npolygons = 1
    do e = 1,size(this%cell_id)
      if (any(this%cell_id(e) == neighbor) .and. this%cell_id(e) /= cell_id) then
        !if (vof(this%cell_id(e)) < 1e-3_r8 .or. vof(this%cell_id(e)) > 1-1e-3_r8) cycle ! WARN?
        npolygons = npolygons + 1
        polygon_id(npolygons) = e
      else if (this%cell_id(e) == cell_id) then
        ! WARNING: see note 1
        polygon_id(1) = e
      end if
    end do

    ! return the local patch
    local_patch = this%element(polygon_id(:npolygons))

    verboseh = .false.
    if (present(verbose)) verboseh = verbose
    if (verboseh) then
      print *, 'RECONSTRUCTION POLYGON IDs'
      do e = 1,npolygons
        print *, e, this%cell_id(polygon_id(e)), local_patch(e)%n_elements
      end do
    end if

  end function local_patch

  ! return the centroid of the polygon box in this cell
  function local_centroid(this, cell_id) result(xc)

    use logging_services

    class(surface), intent(in) :: this
    integer, intent(in) :: cell_id
    real(r8) :: xc(3)

    integer :: e, i
    real(r8) :: area, a

    ! TODO: currently naive search. could do a binary search.
    do e = 1,size(this%cell_id)
      if (this%cell_id(e) == cell_id) then
        xc = 0; area = 0
        do i = 1,this%element(e)%n_elements
          a = this%element(e)%elements(i)%area2()
          xc = xc + this%element(e)%elements(i)%centroid2() * a
          area = area + a
        end do
        xc = xc / area
        exit
      end if
    end do

    if (e > size(this%cell_id)) &
        call LS_fatal ("could not find polygon in this cell")

  end function local_centroid

end module surface_type
