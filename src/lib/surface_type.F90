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

  public :: surface_unit_test

  ! a surface is defined as a collection of polygons
  type, public :: surface
    private
    type(polygon), public, allocatable :: element(:)
  contains
    procedure, private :: append_polyhedron
    procedure, private :: append_polygon
    procedure, private :: append_surf
    generic   :: append => append_polygon, append_surf, append_polyhedron
    procedure :: write_gmv
    procedure :: write_ply
    procedure :: purge
  end type surface

contains

  subroutine surface_unit_test ()
    use polyhedron_type
    use plane_type
    use array_utils, only: mag
    use hex_types,   only: hex_f,hex_e,cube_v

    type(surface)    :: surf
    type(polyhedron) :: cube
    type(plane)      :: P
    type(polygon)    :: element
    integer          :: f,nV

    ! generate a cube and use its faces as the definition of a surface, then print to file
    call cube%init (cube_v, hex_f, hex_e)

    do f = 1,cube%nfaces
      nV = count(cube%face_vid(:,f) /= 0) ! number of vertices on this face
      call element%init (cube%x(:,cube%face_vid(1:nV,f)))

      call surf%append (element)
    end do

    call surf%write_ply ('surf.ply')

    call surf%purge ()
    P%normal = [4.0_r8, 1.0_r8, 1.0_r8]
    P%normal = P%normal / mag (P%normal)
    P%rho    = 0.5_r8 / sqrt(3.0_r8)
    call surf%append (cube%intersection_verts (P))
    call surf%write_ply ('surf.ply')

  end subroutine surface_unit_test

  ! append element to the end of array (increasing array size by 1)
  subroutine append_polygon (this, new_element)
    class(surface), intent(inout) :: this
    class(polygon), intent(in)    :: new_element

    type(polygon), allocatable    :: tmp(:)
    integer                       :: N

    if (new_element%nVerts < 3) return

    if (allocated(this%element)) then
      N = size(this%element)
      tmp = this%element
      deallocate(this%element)
      allocate(this%element(N+1))
      this%element(1:N) = tmp
      this%element(N+1) = new_element
    else
      allocate(this%element(1))
      this%element(1) = new_element
    end if

  end subroutine append_polygon

  ! append a polyhedron to the surface
  ! (may instead want to write a routine in the polyhedron type that returns a surface)
  subroutine append_polyhedron (this, poly)

    use polyhedron_type

    class(surface),    intent(inout) :: this
    class(polyhedron), intent(inout) :: poly

    integer       :: nV,N,f
    type(polygon) :: face,tmp(size(this%element))

    if (allocated(this%element)) then
      N = size(this%element)
      tmp = this%element
      deallocate(this%element)
      allocate(this%element(N+poly%nFaces))
      this%element(1:N) = tmp
    else
      N = 0
      allocate(this%element(poly%nFaces))
    end if

    do f = N+1,N+poly%nFaces
      nV = count(poly%face_vid(:,f) /= 0) ! number of vertices on this face
      call face%init (poly%x(:,poly%face_vid(1:nV,f)), poly%face_normal(:,f))

      this%element(f) = face
    end do
    
  end subroutine append_polyhedron

  ! append elements of another surface to
  subroutine append_surf (this, surf)
    class(surface), intent(inout) :: this
    class(surface), intent(in)    :: surf

    type(polygon)                 :: tmp(size(this%element))
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
  end subroutine purge

  subroutine write_gmv (this)
    use gmvwrite_fortran_binding
    use array_utils, only: xrange

    class(surface), intent(in) :: this

    character(16)              :: filename
    integer                    :: e,j,Nelements,Nverts
    integer,  allocatable      :: element_vid(:,:)
    real(r8), allocatable      :: x(:,:),fakedata(:)

    ! generate contiguous array of nodes, and array of vertex ids
    Nelements = size(this%element)
    Nverts    = sum(this%element(:)%nVerts)

    allocate(x(3,Nverts), element_vid(maxval(this%element(:)%nVerts),size(this%element)), &
         fakedata(Nelements))

    ! To plot the surface elements, gmv needs data there.
    ! Or, at least, paraview won't open the file without it. Paraview won't display
    ! surfaces anyways, though, so it could be paraview that is wrong.
    fakedata = 0.0_r8 

    j = 1
    do e = 1,Nelements
      x(:,j:j+this%element(e)%nVerts-1) = this%element(e)%x
      element_vid(1:this%element(e)%nVerts,e) = xrange (j,j+this%element(e)%nVerts-1)
      j = j + this%element(e)%nVerts
    end do

    ! dump the gmv file
    write(filename,'(a,i4.4)') 'planes.gmv'

    call gmvwrite_openfile_ir_ascii (trim(filename)//C_NULL_CHAR, 4, 8)

    ! write plane nodes
    !call gmv_write_unstr_mesh (this%mesh)
    call gmvwrite_node_data (Nverts, x(1,:), x(2,:), x(3,:))
    call gmvwrite_cell_header (0) ! GMV requires the cell header, even if there are no cells
    call gmvwrite_nodeids (xrange (1,Nverts))

    ! write the plane descriptions
    call gmvwrite_surface_header (size(this%element))
    do e = 1,Nelements
      call gmvwrite_surface_data (this%element(e)%nVerts, element_vid(1:this%element(e)%nVerts,e))
    end do

    call gmvwrite_surfvars_header ()
    call gmvwrite_surfvars_name_data ('surf'//C_NULL_CHAR, fakedata)
    call gmvwrite_surfvars_endsvar ()

    call gmvwrite_closefile ()

    deallocate (x, element_vid, fakedata)

  end subroutine write_gmv

  ! write the polygons using the Stanford PLY format
  subroutine write_ply (this, fname)
    class(surface), intent(inout) :: this
    character(*),   intent(in)    :: fname

    integer                       :: e,j,i,Nelements,Nverts

    if (.not.allocated(this%element)) return

    Nelements = size(this%element)
    Nverts    = sum(this%element(:)%nVerts)

    ! open file
    open(99, file=trim(fname))

    ! write PLY header
    write(99,'(a)')    'ply'
    write(99,'(a)')    'format ascii 1.0'
    write(99,'(a,i9)') 'element vertex ',Nverts
    write(99,'(a)')    'property float32 x'
    write(99,'(a)')    'property float32 y'
    write(99,'(a)')    'property float32 z'
    write(99,'(a,i9)') 'element face ',Nelements
    write(99,'(a)')    'property list uint8 int32 vertex_index'
    write(99,'(a)')    'end_header'

    ! write polygon data
    ! vertices
    do e = 1,Nelements
      do j = 1,this%element(e)%nVerts
        write(99,'(3f20.10)') this%element(e)%x(:,j)
      end do
    end do

    ! faces
    j = 0
    do e = 1,Nelements
      write(99,'(i9)',advance='no') this%element(e)%nVerts
      do i = j,j+this%element(e)%nVerts-1
        write(99,'(i9)',advance='no') i
      end do
      write(99,*)
      j = j+this%element(e)%nVerts
    end do

    ! clean up
    close (99)

  end subroutine write_ply

end module surface_type
