!!
!! VOF_INIT
!!
!! This module provides routines for initializing the
!! VOF scalar based on shapes provided by the user
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! revised Nov 2015
!!

module vof_init

  use kinds, only: r8
  use material_geometry_type
  use logging_services
  use hex_types, only: hex_f, hex_e
  use polyhedron_type
  implicit none
  private

  ! hex type for the divide and conquer algorithm
  ! public to expose for testing
  type, extends(polyhedron), public :: dnc_hex
    private
    real(r8) :: cellc(3), facec(3,6), edgec(3,12), weight
    integer :: mcellc, mfacec(6), medgec(12), matl_at_node(8)
  contains
    procedure :: subcell
    procedure :: set_division_variables
    procedure :: cell_center
    procedure :: face_centers
    procedure :: edge_centers
    procedure :: contains_interface
    procedure :: vof
    procedure, private :: vof_from_nodes
  end type dnc_hex

  public :: vof_initialize

  ! TODO: make these user-specified parameters
  integer , parameter :: cell_vof_recursion_limit = 8

contains

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
    type(dnc_hex)           :: hex
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
      ! initialize dnc_hex
      call hex%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, gmesh%outnorm(:,:,i), &
          mesh%volume(i))
      do v = 1,8
        hex%matl_at_node(v) = matl_init_geometry%index_at(hex%x(:,v))
      end do

      ! calculate the vof
      vof(:,i) = hex%vof (matl_init_geometry, nmat, 0)
    end do
    !$omp end parallel do

    call stop_timer("vof init")
    print '(a)', "done"

  end subroutine vof_initialize

  ! calculates the volume fractions of materials in a cell
  recursive function vof (this, matl_geometry, nmat, depth) result(hex_vof)

    class(dnc_hex),     intent(inout) :: this
    class(base_region), intent(in) :: matl_geometry
    integer,            intent(in) :: nmat,depth
    real(r8)                       :: hex_vof(nmat)

    integer       :: i
    type(dnc_hex) :: subcell

    ! if the cell contains an interface (and therefore has at least two materials
    ! and we haven't yet hit our recursion limit, divide the hex and repeat
    if (this%contains_interface(matl_geometry) .and. depth < cell_vof_recursion_limit) then
      ! tally the vof from sub-cells
      call this%set_division_variables(matl_geometry)
      hex_vof = 0
      do i = 1,8
        subcell = this%subcell(matl_geometry, i)
        hex_vof = hex_vof + subcell%vof(matl_geometry,nmat,depth+1) * subcell%weight
      end do
    else
      ! if we are past the recursion limit
      ! or the cell does not contain an interface, calculate
      ! the vof in this hex based on the materials at its nodes
      hex_vof = this%vof_from_nodes(matl_geometry, nmat)
    end if

  end function vof

  ! This function takes a hex, described by its vertex positions,
  ! and determines if it contains an interface for a given material.
  ! This is done by checking if every vertex lies within the material, or not.
  pure logical function contains_interface(this, matl_geometry)

    class(dnc_hex),     intent(in) :: this
    class(base_region), intent(in) :: matl_geometry

    contains_interface = any(this%matl_at_node /= this%matl_at_node(1))

  end function contains_interface

  ! the material at each node contributes 1/8th of the vof in the given hex
  pure function vof_from_nodes (this, matl_geometry, nmat)

    class(dnc_hex),     intent(in) :: this
    class(base_region), intent(in) :: matl_geometry
    integer,            intent(in) :: nmat
    real(r8)                       :: vof_from_nodes(nmat)

    integer :: m

    do m = 1,nmat
      vof_from_nodes(m) = count(this%matl_at_node == m) / 8.0_r8
    end do

  end function vof_from_nodes

  ! returns the center of the given hex
  pure function cell_center (hex)
    class(dnc_hex), intent(in) :: hex
    real(r8) :: cell_center(3)
    cell_center = sum(hex%x, dim=2) / 8
  end function cell_center

  ! returns an array of face centers
  pure function face_centers(this)

    class(dnc_hex), intent(in) :: this
    real(r8) :: face_centers(3,this%nFaces)

    integer :: f

    do f = 1,this%nFaces
      face_centers(:,f) = sum(this%x(:,this%face_vid(:,f)), dim=2) / 4
    end do

  end function face_centers

  ! returns an array of edge centers
  pure function edge_centers(this)

    class(dnc_hex), intent(in) :: this
    real(r8) :: edge_centers(3,this%nEdges)

    integer :: e

    do e = 1,this%nEdges
      edge_centers(:,e) = sum(this%x(:,this%edge_vid(:,e)), dim=2) / 2
    end do

  end function edge_centers

  subroutine set_division_variables (this, matl_geometry)

    class(dnc_hex), intent(inout) :: this
    class(base_region), intent(in) :: matl_geometry

    integer :: j

    ! find center points
    this%cellc = this%cell_center()
    this%facec = this%face_centers()
    this%edgec = this%edge_centers()

    ! get material ids at all points above
    this%mcellc = matl_geometry%index_at(this%cellc)
    do j = 1,6
      this%mfacec(j) = matl_geometry%index_at(this%facec(:,j))
    end do
    do j = 1,12
      this%medgec(j) = matl_geometry%index_at(this%edgec(:,j))
    end do

  end subroutine set_division_variables

  type(dnc_hex) function subcell (this, matl_geometry, i)

    class(dnc_hex), intent(inout) :: this
    class(base_region), intent(in) :: matl_geometry
    integer, intent(in) :: i

    real(r8) :: xtmp(3,8)
    integer :: matl_at_node(8), ierr

    ! set subhex node positions
    select case(i)
    case(1)
      xtmp(:,1) = this%x(:,1)
      xtmp(:,2) = this%edgec(:,1)
      xtmp(:,3) = this%facec(:,5)
      xtmp(:,4) = this%edgec(:,4)
      xtmp(:,5) = this%edgec(:,5)
      xtmp(:,6) = this%facec(:,2)
      xtmp(:,7) = this%cellc(:)
      xtmp(:,8) = this%facec(:,3)

      matl_at_node(1) = this%matl_at_node(1)
      matl_at_node(2) = this%medgec(1)
      matl_at_node(3) = this%mfacec(5)
      matl_at_node(4) = this%medgec(4)
      matl_at_node(5) = this%medgec(5)
      matl_at_node(6) = this%mfacec(2)
      matl_at_node(7) = this%mcellc
      matl_at_node(8) = this%mfacec(3)

    case(2)
      xtmp(:,1) = this%edgec(:,1)
      xtmp(:,2) = this%x(:,2)
      xtmp(:,3) = this%edgec(:,2)
      xtmp(:,4) = this%facec(:,5)
      xtmp(:,5) = this%facec(:,2)
      xtmp(:,6) = this%edgec(:,6)
      xtmp(:,7) = this%facec(:,4)
      xtmp(:,8) = this%cellc(:)

      matl_at_node(1) = this%medgec(1)
      matl_at_node(2) = this%matl_at_node(2)
      matl_at_node(3) = this%medgec(2)
      matl_at_node(4) = this%mfacec(5)
      matl_at_node(5) = this%mfacec(2)
      matl_at_node(6) = this%medgec(6)
      matl_at_node(7) = this%mfacec(4)
      matl_at_node(8) = this%mcellc

    case(3)
      xtmp(:,1) = this%facec(:,5)
      xtmp(:,2) = this%edgec(:,2)
      xtmp(:,3) = this%x(:,3)
      xtmp(:,4) = this%edgec(:,3)
      xtmp(:,5) = this%cellc(:)
      xtmp(:,6) = this%facec(:,4)
      xtmp(:,7) = this%edgec(:,7)
      xtmp(:,8) = this%facec(:,1)

      matl_at_node(1) = this%mfacec(5)
      matl_at_node(2) = this%medgec(2)
      matl_at_node(3) = this%matl_at_node(3)
      matl_at_node(4) = this%medgec(3)
      matl_at_node(5) = this%mcellc
      matl_at_node(6) = this%mfacec(4)
      matl_at_node(7) = this%medgec(7)
      matl_at_node(8) = this%mfacec(1)

    case(4)
      xtmp(:,1) = this%edgec(:,4)
      xtmp(:,2) = this%facec(:,5)
      xtmp(:,3) = this%edgec(:,3)
      xtmp(:,4) = this%x(:,4)
      xtmp(:,5) = this%facec(:,3)
      xtmp(:,6) = this%cellc(:)
      xtmp(:,7) = this%facec(:,1)
      xtmp(:,8) = this%edgec(:,8)

      matl_at_node(1) = this%medgec(4)
      matl_at_node(2) = this%mfacec(5)
      matl_at_node(3) = this%medgec(3)
      matl_at_node(4) = this%matl_at_node(4)
      matl_at_node(5) = this%mfacec(3)
      matl_at_node(6) = this%mcellc
      matl_at_node(7) = this%mfacec(1)
      matl_at_node(8) = this%medgec(8)

    case(5)
      xtmp(:,1) = this%edgec(:,5)
      xtmp(:,2) = this%facec(:,2)
      xtmp(:,3) = this%cellc(:)
      xtmp(:,4) = this%facec(:,3)
      xtmp(:,5) = this%x(:,5)
      xtmp(:,6) = this%edgec(:,9)
      xtmp(:,7) = this%facec(:,6)
      xtmp(:,8) = this%edgec(:,12)

      matl_at_node(1) = this%medgec(5)
      matl_at_node(2) = this%mfacec(2)
      matl_at_node(3) = this%mcellc
      matl_at_node(4) = this%mfacec(3)
      matl_at_node(5) = this%matl_at_node(5)
      matl_at_node(6) = this%medgec(9)
      matl_at_node(7) = this%mfacec(6)
      matl_at_node(8) = this%medgec(12)

    case(6)
      xtmp(:,1) = this%facec(:,2)
      xtmp(:,2) = this%edgec(:,6)
      xtmp(:,3) = this%facec(:,4)
      xtmp(:,4) = this%cellc(:)
      xtmp(:,5) = this%edgec(:,9)
      xtmp(:,6) = this%x(:,6)
      xtmp(:,7) = this%edgec(:,10)
      xtmp(:,8) = this%facec(:,6)

      matl_at_node(1) = this%mfacec(2)
      matl_at_node(2) = this%medgec(6)
      matl_at_node(3) = this%mfacec(4)
      matl_at_node(4) = this%mcellc
      matl_at_node(5) = this%medgec(9)
      matl_at_node(6) = this%matl_at_node(6)
      matl_at_node(7) = this%medgec(10)
      matl_at_node(8) = this%mfacec(6)

    case(7)
      xtmp(:,1) = this%cellc(:)
      xtmp(:,2) = this%facec(:,4)
      xtmp(:,3) = this%edgec(:,7)
      xtmp(:,4) = this%facec(:,1)
      xtmp(:,5) = this%facec(:,6)
      xtmp(:,6) = this%edgec(:,10)
      xtmp(:,7) = this%x(:,7)
      xtmp(:,8) = this%edgec(:,11)

      matl_at_node(1) = this%mcellc
      matl_at_node(2) = this%mfacec(4)
      matl_at_node(3) = this%medgec(7)
      matl_at_node(4) = this%mfacec(1)
      matl_at_node(5) = this%mfacec(6)
      matl_at_node(6) = this%medgec(10)
      matl_at_node(7) = this%matl_at_node(7)
      matl_at_node(8) = this%medgec(11)

    case(8)
      xtmp(:,1) = this%facec(:,3)
      xtmp(:,2) = this%cellc(:)
      xtmp(:,3) = this%facec(:,1)
      xtmp(:,4) = this%edgec(:,8)
      xtmp(:,5) = this%edgec(:,12)
      xtmp(:,6) = this%facec(:,6)
      xtmp(:,7) = this%edgec(:,11)
      xtmp(:,8) = this%x(:,8)

      matl_at_node(1) = this%mfacec(3)
      matl_at_node(2) = this%mcellc
      matl_at_node(3) = this%mfacec(1)
      matl_at_node(4) = this%medgec(8)
      matl_at_node(5) = this%medgec(12)
      matl_at_node(6) = this%mfacec(6)
      matl_at_node(7) = this%medgec(11)
      matl_at_node(8) = this%matl_at_node(8)
    end select

    call subcell%init (ierr, xtmp, hex_f, hex_e)
    subcell%matl_at_node = matl_at_node

    ! for nonorthogonal meshes, sub-cell volumes are
    ! not equal, so the weight is the ratio of volumes
    subcell%weight = subcell%volume() / this%volume()
    !subcell%weight = 1.0_r8 / 8.0_r8

    if (subcell%weight==0) then
      call LS_fatal("zero weight")
    end if

  end function subcell

  end module vof_init
