!!
!! VOF_INIT
!!
!! This module provides routines for initializing the
!! VOF scalar based on shapes provided by the user
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! revised Nov 2015
!!

#include "f90_assert.fpp"

module vof_init

  use kinds, only: r8
  use material_geometry_type
  use logging_services
  use hex_types, only: hex_f, hex_e
  use polyhedron_type
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan ! DEBUGGING
  implicit none
  private

  ! hex type for the divide and conquer algorithm
  ! public to expose for testing
  type, extends(polyhedron), public :: dnc_cell
    private
    real(r8), allocatable :: cellc(:), facec(:,:), edgec(:,:)
    integer, allocatable :: mfacec(:), medgec(:), matl_at_node(:)
    integer :: mcellc, ntets
  contains
    procedure :: subcell
    procedure :: set_division_variables
    procedure :: cell_center
    procedure :: face_centers
    procedure :: edge_centers
    procedure :: contains_interface
    procedure :: vof
    procedure, private :: vof_from_nodes
    procedure, private :: init_polyhedron_copy => dnc_cell_copy
    procedure, private :: dnc_cell_copy
    !generic :: assignment(=) => dnc_cell_copy
  end type dnc_cell

  public :: vof_initialize

  ! TODO: make these user-specified parameters
  integer , parameter :: cell_vof_recursion_limit = 8

contains

  subroutine dnc_cell_copy (this, poly)

    class(dnc_cell), intent(out) :: this
    class(polyhedron), intent(in) :: poly

    ! parent variables
    this%nVerts = poly%nVerts
    this%nEdges = poly%nEdges
    this%nFaces = poly%nFaces
    this%vol = poly%vol
    this%tesselated = poly%tesselated

    if (allocated(poly%x)) this%x = poly%x
    if (allocated(poly%edge_vid)) this%edge_vid = poly%edge_vid
    if (allocated(poly%face_vid)) this%face_vid = poly%face_vid
    if (allocated(poly%face_normal)) this%face_normal = poly%face_normal
    if (allocated(poly%vertex_faces)) this%vertex_faces = poly%vertex_faces
    if (allocated(poly%edge_faces)) this%edge_faces = poly%edge_faces
    if (allocated(poly%face_eid)) this%face_eid = poly%face_eid

    if (associated(poly%tet)) allocate(this%tet(size(poly%tet)), source=poly%tet)

    select type(poly)
    type is (dnc_cell)
      ! child variables
      if (allocated(poly%cellc)) this%cellc = poly%cellc
      if (allocated(poly%facec)) this%facec = poly%facec
      if (allocated(poly%edgec)) this%edgec = poly%edgec
      !this%weight = poly%weight
      this%mcellc = poly%mcellc
      this%ntets = poly%ntets
      if (allocated(poly%mfacec)) this%mfacec = poly%mfacec
      if (allocated(poly%medgec)) this%medgec = poly%medgec
      if (allocated(poly%matl_at_node)) this%matl_at_node = poly%matl_at_node
    end select

  end subroutine dnc_cell_copy

  ! this subroutine initializes the vof values in all cells
  ! from a user input function. It calls a divide and conquer algorithm
  ! to calculate the volume fractions given an interface
  subroutine vof_initialize (mesh, gmesh, plist, vof, matl_ids, nmat) !cell_matls)

    use unstr_mesh_type
    use mesh_geom_type
    use parameter_list_type
    use timer_tree_type

    type(unstr_mesh),     intent(in)    :: mesh
    type(mesh_geom), intent(in) :: gmesh
    type(parameter_list), intent(in)    :: plist
    integer,              intent(in)    :: nmat
    integer,              intent(in)    :: matl_ids(:)
    real(r8),             intent(inout) :: vof(:,:)

    integer                 :: i,v, ierr
    type(dnc_cell)          :: cell
    type(material_geometry) :: matl_init_geometry

    print '(a)', "initializing vof ... "
    call start_timer("vof init")

    ! first, determine the initial state provided by the user, and assign
    ! the function which will determine what points are inside the materials
    call matl_init_geometry%init (plist, matl_ids)

    ! next, loop though every cell and check if all vertices lie within a single material
    ! if so, set the Vof for that material to 1.0.
    ! if not, divide the hex cell into 8 subdomains defined by the centroid and centers of faces
    ! repeat the process recursively for each subdomain up to a given threshold

    !$omp parallel do default(private) shared(vof,mesh,gmesh,matl_init_geometry,nmat)
    do i = 1,mesh%ncell
      ! initialize dnc_cell
      call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, gmesh%outnorm(:,:,i), &
          mesh%volume(i), tesselate=.false.)
      allocate(cell%matl_at_node(cell%nVerts))
      do v = 1,cell%nVerts
        cell%matl_at_node(v) = matl_init_geometry%index_at(cell%x(:,v))
      end do
      if (any(cell%matl_at_node==0)) call LS_fatal("parent node materials unset")

      ! calculate the vof
      vof(:,i) = cell%vof (matl_init_geometry, nmat, 0)

      !if (all(vof(:,i) /= 1)) print *, vof(:,i)
      if (any(ieee_is_nan(vof(:,i))) .or. abs(sum(vof(:,i)) - 1) > 1e-12_r8) then
        !print *, this%contains_interface(matl_geometry) .and. depth < cell_vof_recursion_limit
        call LS_Fatal("invalid vof init")
      end if
    end do
    !$omp end parallel do

    call stop_timer("vof init")
    print '(a)', "done"

  end subroutine vof_initialize

  ! calculates the volume fractions of materials in a cell
  recursive function vof (this, matl_geometry, nmat, depth) result(hex_vof)

    class(dnc_cell),     intent(inout) :: this
    class(base_region), intent(in) :: matl_geometry
    integer,            intent(in) :: nmat,depth
    real(r8)                       :: hex_vof(nmat)

    real(r8) :: sum_vol
    integer :: i
    type(dnc_cell) :: subcell

    ! if the cell contains an interface (and therefore has at least two materials
    ! and we haven't yet hit our recursion limit, divide the hex and repeat
    if (this%contains_interface(matl_geometry) .and. depth < cell_vof_recursion_limit) then
      ! tally the vof from sub-cells
      call this%set_division_variables(matl_geometry)
      hex_vof = 0; sum_vol = 0
      do i = 1,this%ntets
        call subcell%dnc_cell_copy(this%subcell(matl_geometry, i))
        if (any(subcell%matl_at_node==0) .or. .not.allocated(subcell%x) &
            .or. subcell%volume()==0) then
          print *, depth
          print *, subcell%matl_at_node
          call LS_Fatal("invalid subcell")
        end if
        hex_vof = hex_vof + subcell%vof(matl_geometry,nmat,depth+1) * subcell%volume()
        sum_vol = sum_vol + subcell%volume()
      end do
      hex_vof = hex_vof / sum_vol
    else
      ! if we are past the recursion limit
      ! or the cell does not contain an interface, calculate
      ! the vof in this hex based on the materials at its nodes
      hex_vof = this%vof_from_nodes(matl_geometry, nmat)
    end if

    if (any(ieee_is_nan(hex_vof)) .or. all(hex_vof==0) .or. abs(sum(hex_vof) - 1) > 1e-10_r8) then
      print *, hex_vof, sum(hex_vof)
      print *, depth, sum_vol
      print *, this%volume()
      this%vol = 0
      print *, this%volume()
      subcell = this
      subcell%vol = 0
      call subcell%tesselate()
      print *, subcell%volume(), sum_vol * this%volume() / subcell%volume()
      print *, this%contains_interface(matl_geometry) .and. depth < cell_vof_recursion_limit
      print *, this%matl_at_node
      print *, hex_vof
      call LS_Fatal("invalid init vof")
    end if

  end function vof

  ! This function takes a hex, described by its vertex positions,
  ! and determines if it contains an interface for a given material.
  ! This is done by checking if every vertex lies within the material, or not.
  pure logical function contains_interface(this, matl_geometry)

    class(dnc_cell),     intent(in) :: this
    class(base_region), intent(in) :: matl_geometry

    contains_interface = any(this%matl_at_node /= this%matl_at_node(1))

  end function contains_interface

  ! the material at each node contributes 1/8th of the vof in the given hex
  function vof_from_nodes (this, matl_geometry, nmat)

    class(dnc_cell),     intent(in) :: this
    class(base_region), intent(in) :: matl_geometry
    integer,            intent(in) :: nmat
    real(r8)                       :: vof_from_nodes(nmat)

    integer :: m

    !if (any(this%matl_at_node==0)) call LS_fatal("node materials unset")

    vof_from_nodes = 0
    do m = 1,nmat
      vof_from_nodes(m) = count(this%matl_at_node == m)
    end do
    vof_from_nodes = vof_from_nodes / this%nVerts

  end function vof_from_nodes

  ! returns the center of the given hex
  pure function cell_center (this)
    class(dnc_cell), intent(in) :: this
    real(r8) :: cell_center(3)
    cell_center = sum(this%x, dim=2) / this%nVerts
  end function cell_center

  ! returns an array of face centers
  pure function face_centers(this)

    class(dnc_cell), intent(in) :: this
    real(r8) :: face_centers(3,this%nFaces)

    integer :: f, nV

    do f = 1,this%nFaces
      nV = count(this%face_vid(:,f) > 0)
      face_centers(:,f) = sum(this%x(:,this%face_vid(:nV,f)), dim=2) / nV
    end do

  end function face_centers

  ! returns an array of edge centers
  pure function edge_centers(this)

    class(dnc_cell), intent(in) :: this
    real(r8) :: edge_centers(3,this%nEdges)

    integer :: e

    do e = 1,this%nEdges
      edge_centers(:,e) = sum(this%x(:,this%edge_vid(:,e)), dim=2) / 2
    end do

  end function edge_centers

  subroutine set_division_variables (this, matl_geometry)

    class(dnc_cell), intent(inout) :: this
    class(base_region), intent(in) :: matl_geometry

    integer :: j

    ! find center points
    this%cellc = this%cell_center()
    this%facec = this%face_centers()
    this%edgec = this%edge_centers()

    ! get material ids at all points above
    allocate(this%mfacec(this%nFaces), this%medgec(this%nEdges))
    this%mcellc = matl_geometry%index_at(this%cellc)
    do j = 1,this%nFaces
      this%mfacec(j) = matl_geometry%index_at(this%facec(:,j))
    end do
    do j = 1,this%nEdges
      this%medgec(j) = matl_geometry%index_at(this%edgec(:,j))
    end do

    this%ntets = 2*count(this%face_vid > 0)

  end subroutine set_division_variables

  type(dnc_cell) function subcell (this, matl_geometry, i)

    class(dnc_cell), intent(inout) :: this
    class(base_region), intent(in) :: matl_geometry
    integer, intent(in) :: i

    real(r8) :: xtmp(3,8)
    integer :: matl_at_node(8), ierr

    integer :: nV, f, v, feid, e, edge_side
    !integer :: t, ff, ee

    ! get face, face_edge_id, and edge_side of this tet
    ! assume all faces have the same number of nodes/edges
    ! i = (f-1)*2*nV + (feid-1)*2 + edge_side
    nV = count(this%face_vid(:,1) > 0)

    ! t = 0
    ! do ff = 1,this%nFaces
    !   do ee = 1,nV
    !     t = t+1

    !     f = (t-1) / (2*nV) + 1
    !     feid = (t - 1 - (f-1)*2*nV) / 2 + 1
    !     edge_side = t - (f-1)*2*nV - (feid-1)*2
    !     print *, ff, ee, 1
    !     print *, f, feid, edge_side
    !     print *, f==ff .and. ee==feid .and. edge_side==1
    !     print *


    !     t = t+1

    !     f = (t-1) / (2*nV) + 1
    !     feid = (t - 1 - (f-1)*2*nV) / 2 + 1
    !     edge_side = t - (f-1)*2*nV - (feid-1)*2
    !     print *, ff, ee, 2
    !     print *, f, feid, edge_side
    !     print *, f==ff .and. ee==feid .and. edge_side==2
    !     print *
    !   end do
    ! end do

    f = (i-1) / (2*nV) + 1
    feid = (i - 1 - (f-1)*2*nV) / 2 + 1
    edge_side = i - (f-1)*2*nV - (feid-1)*2
    ! print *, i, f, feid, edge_side
    ! stop
    e = this%face_eid(feid,f)

    xtmp(:,1) = this%cellc
    xtmp(:,2) = this%facec(:,f)
    if (edge_side==1) then
      xtmp(:,3) = this%x(:,this%face_vid(feid,f))
      xtmp(:,4) = this%edgec(:,e)
    else
      xtmp(:,3) = this%edgec(:,e)
      xtmp(:,4) = this%x(:,this%face_vid(modulo(feid,nV)+1,f))
    end if

    matl_at_node(1) = this%mcellc
    matl_at_node(2) = this%mfacec(f)
    if (edge_side==1) then
      matl_at_node(3) = this%matl_at_node(this%face_vid(feid,f))
      matl_at_node(4) = this%medgec(e)
    else
      matl_at_node(3) = this%medgec(e)
      matl_at_node(4) = this%matl_at_node(this%face_vid(modulo(feid,nV)+1,f))
    end if

    call subcell%init (ierr, xtmp(:,:4))
    subcell%matl_at_node = matl_at_node(:4)

    ! print *, 'here0', allocated(subcell%x), subcell%tesselated

    ! print *, i, f, feid, edge_side, e

    ! print *, xtmp(:,:4)
    ! call subcell%print_data()
    ! print *, 'here1', allocated(subcell%x), subcell%tesselated

    ! ! set subhex node positions
    ! select case(i)
    ! case(1)
    !   xtmp(:,1) = this%x(:,1)
    !   xtmp(:,2) = this%edgec(:,1)
    !   xtmp(:,3) = this%facec(:,5)
    !   xtmp(:,4) = this%edgec(:,4)
    !   xtmp(:,5) = this%edgec(:,5)
    !   xtmp(:,6) = this%facec(:,2)
    !   xtmp(:,7) = this%cellc(:)
    !   xtmp(:,8) = this%facec(:,3)

    !   matl_at_node(1) = this%matl_at_node(1)
    !   matl_at_node(2) = this%medgec(1)
    !   matl_at_node(3) = this%mfacec(5)
    !   matl_at_node(4) = this%medgec(4)
    !   matl_at_node(5) = this%medgec(5)
    !   matl_at_node(6) = this%mfacec(2)
    !   matl_at_node(7) = this%mcellc
    !   matl_at_node(8) = this%mfacec(3)

    ! case(2)
    !   xtmp(:,1) = this%edgec(:,1)
    !   xtmp(:,2) = this%x(:,2)
    !   xtmp(:,3) = this%edgec(:,2)
    !   xtmp(:,4) = this%facec(:,5)
    !   xtmp(:,5) = this%facec(:,2)
    !   xtmp(:,6) = this%edgec(:,6)
    !   xtmp(:,7) = this%facec(:,4)
    !   xtmp(:,8) = this%cellc(:)

    !   matl_at_node(1) = this%medgec(1)
    !   matl_at_node(2) = this%matl_at_node(2)
    !   matl_at_node(3) = this%medgec(2)
    !   matl_at_node(4) = this%mfacec(5)
    !   matl_at_node(5) = this%mfacec(2)
    !   matl_at_node(6) = this%medgec(6)
    !   matl_at_node(7) = this%mfacec(4)
    !   matl_at_node(8) = this%mcellc

    ! case(3)
    !   xtmp(:,1) = this%facec(:,5)
    !   xtmp(:,2) = this%edgec(:,2)
    !   xtmp(:,3) = this%x(:,3)
    !   xtmp(:,4) = this%edgec(:,3)
    !   xtmp(:,5) = this%cellc(:)
    !   xtmp(:,6) = this%facec(:,4)
    !   xtmp(:,7) = this%edgec(:,7)
    !   xtmp(:,8) = this%facec(:,1)

    !   matl_at_node(1) = this%mfacec(5)
    !   matl_at_node(2) = this%medgec(2)
    !   matl_at_node(3) = this%matl_at_node(3)
    !   matl_at_node(4) = this%medgec(3)
    !   matl_at_node(5) = this%mcellc
    !   matl_at_node(6) = this%mfacec(4)
    !   matl_at_node(7) = this%medgec(7)
    !   matl_at_node(8) = this%mfacec(1)

    ! case(4)
    !   xtmp(:,1) = this%edgec(:,4)
    !   xtmp(:,2) = this%facec(:,5)
    !   xtmp(:,3) = this%edgec(:,3)
    !   xtmp(:,4) = this%x(:,4)
    !   xtmp(:,5) = this%facec(:,3)
    !   xtmp(:,6) = this%cellc(:)
    !   xtmp(:,7) = this%facec(:,1)
    !   xtmp(:,8) = this%edgec(:,8)

    !   matl_at_node(1) = this%medgec(4)
    !   matl_at_node(2) = this%mfacec(5)
    !   matl_at_node(3) = this%medgec(3)
    !   matl_at_node(4) = this%matl_at_node(4)
    !   matl_at_node(5) = this%mfacec(3)
    !   matl_at_node(6) = this%mcellc
    !   matl_at_node(7) = this%mfacec(1)
    !   matl_at_node(8) = this%medgec(8)

    ! case(5)
    !   xtmp(:,1) = this%edgec(:,5)
    !   xtmp(:,2) = this%facec(:,2)
    !   xtmp(:,3) = this%cellc(:)
    !   xtmp(:,4) = this%facec(:,3)
    !   xtmp(:,5) = this%x(:,5)
    !   xtmp(:,6) = this%edgec(:,9)
    !   xtmp(:,7) = this%facec(:,6)
    !   xtmp(:,8) = this%edgec(:,12)

    !   matl_at_node(1) = this%medgec(5)
    !   matl_at_node(2) = this%mfacec(2)
    !   matl_at_node(3) = this%mcellc
    !   matl_at_node(4) = this%mfacec(3)
    !   matl_at_node(5) = this%matl_at_node(5)
    !   matl_at_node(6) = this%medgec(9)
    !   matl_at_node(7) = this%mfacec(6)
    !   matl_at_node(8) = this%medgec(12)

    ! case(6)
    !   xtmp(:,1) = this%facec(:,2)
    !   xtmp(:,2) = this%edgec(:,6)
    !   xtmp(:,3) = this%facec(:,4)
    !   xtmp(:,4) = this%cellc(:)
    !   xtmp(:,5) = this%edgec(:,9)
    !   xtmp(:,6) = this%x(:,6)
    !   xtmp(:,7) = this%edgec(:,10)
    !   xtmp(:,8) = this%facec(:,6)

    !   matl_at_node(1) = this%mfacec(2)
    !   matl_at_node(2) = this%medgec(6)
    !   matl_at_node(3) = this%mfacec(4)
    !   matl_at_node(4) = this%mcellc
    !   matl_at_node(5) = this%medgec(9)
    !   matl_at_node(6) = this%matl_at_node(6)
    !   matl_at_node(7) = this%medgec(10)
    !   matl_at_node(8) = this%mfacec(6)

    ! case(7)
    !   xtmp(:,1) = this%cellc(:)
    !   xtmp(:,2) = this%facec(:,4)
    !   xtmp(:,3) = this%edgec(:,7)
    !   xtmp(:,4) = this%facec(:,1)
    !   xtmp(:,5) = this%facec(:,6)
    !   xtmp(:,6) = this%edgec(:,10)
    !   xtmp(:,7) = this%x(:,7)
    !   xtmp(:,8) = this%edgec(:,11)

    !   matl_at_node(1) = this%mcellc
    !   matl_at_node(2) = this%mfacec(4)
    !   matl_at_node(3) = this%medgec(7)
    !   matl_at_node(4) = this%mfacec(1)
    !   matl_at_node(5) = this%mfacec(6)
    !   matl_at_node(6) = this%medgec(10)
    !   matl_at_node(7) = this%matl_at_node(7)
    !   matl_at_node(8) = this%medgec(11)

    ! case(8)
    !   xtmp(:,1) = this%facec(:,3)
    !   xtmp(:,2) = this%cellc(:)
    !   xtmp(:,3) = this%facec(:,1)
    !   xtmp(:,4) = this%edgec(:,8)
    !   xtmp(:,5) = this%edgec(:,12)
    !   xtmp(:,6) = this%facec(:,6)
    !   xtmp(:,7) = this%edgec(:,11)
    !   xtmp(:,8) = this%x(:,8)

    !   matl_at_node(1) = this%mfacec(3)
    !   matl_at_node(2) = this%mcellc
    !   matl_at_node(3) = this%mfacec(1)
    !   matl_at_node(4) = this%medgec(8)
    !   matl_at_node(5) = this%medgec(12)
    !   matl_at_node(6) = this%mfacec(6)
    !   matl_at_node(7) = this%medgec(11)
    !   matl_at_node(8) = this%matl_at_node(8)
    ! case default
    !   call LS_fatal("requested invalid subhex")
    ! end select

    ! call subcell%init (ierr, xtmp, hex_f, hex_e, tesselate=.false.)
    ! subcell%matl_at_node = matl_at_node

    ! for nonorthogonal meshes, sub-cell volumes are
    ! not equal, so the weight is the ratio of volumes
    !subcell%weight = subcell%volume() / this%volume()
    !subcell%weight = 1.0_r8 / 8.0_r8

    !if (subcell%weight==0) call LS_fatal("zero weight")
    if (any(subcell%matl_at_node==0)) call LS_fatal("node materials unset")
    if (subcell%volume()==0) then
      call this%print_data()
      call LS_fatal("zero subcell volume")
    end if
    !print *, subcell%matl_at_node

  end function subcell

end module vof_init
