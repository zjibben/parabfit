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
  use polyhedron_type
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan ! DEBUGGING
  implicit none
  private

  ! hex type for the divide and conquer algorithm
  ! public to expose for testing
  type, public :: dnc_cell
    private
    type(polyhedron) :: geom
    real(r8), allocatable :: cellc(:), facec(:,:), edgec(:,:)
    integer, allocatable :: mfacec(:), medgec(:), matl_at_node(:)
    integer :: mcellc, ntets
  contains
    procedure, private :: dnc_cell_init
    procedure, private :: init_from_polyhedron
    generic :: init => dnc_cell_init, init_from_polyhedron
    procedure :: get_subcell
    procedure :: set_division_variables
    procedure :: cell_center
    procedure :: face_centers
    procedure :: edge_centers
    procedure :: contains_interface
    procedure :: vof
    procedure, private :: vof_from_nodes
  end type dnc_cell

  public :: vof_initialize

  integer, parameter :: cell_vof_recursion_limit = 5

contains

  subroutine dnc_cell_init(this, ierr, matl_init_geometry, i, mesh, gmesh, tesselate)

    use unstr_mesh_type
    use mesh_geom_type

    class(dnc_cell), intent(out) :: this
    integer, intent(out) :: ierr
    type(material_geometry), intent(in) :: matl_init_geometry
    integer, intent(in) :: i
    class(unstr_mesh), intent(in) :: mesh
    class(mesh_geom), intent(in) :: gmesh
    logical, optional, intent(in) :: tesselate

    integer :: v

    call this%geom%init(ierr, i, mesh, gmesh, tesselate)
    !call this%geom%init(ierr, x, face_v, edge_v, face_normal, vol, tesselate)

    allocate(this%matl_at_node(this%geom%parent%nVerts))
    do v = 1,this%geom%parent%nVerts
      this%matl_at_node(v) = matl_init_geometry%index_at(this%geom%parent%x(:,v))
    end do
    if (any(this%matl_at_node==0)) call LS_fatal("parent node materials unset")

  end subroutine dnc_cell_init

  subroutine init_from_polyhedron(this, ierr, matl_init_geometry, x, face_v, edge_v, face_normal, &
      vol, tesselate)

    use unstr_mesh_type
    use mesh_geom_type

    class(dnc_cell), intent(out) :: this
    integer, intent(out) :: ierr
    type(material_geometry), intent(in) :: matl_init_geometry
    real(r8),           intent(in)  :: x(:,:)
    integer,            intent(in)  :: face_v(:,:), edge_v(:,:)
    real(r8), optional, intent(in)  :: face_normal(:,:), vol
    logical, optional, intent(in) :: tesselate

    integer :: v

    call this%geom%init(ierr, x, face_v, edge_v, face_normal, vol, tesselate)

    allocate(this%matl_at_node(this%geom%parent%nVerts))
    do v = 1,this%geom%parent%nVerts
      this%matl_at_node(v) = matl_init_geometry%index_at(this%geom%parent%x(:,v))
    end do
    if (any(this%matl_at_node==0)) call LS_fatal("parent node materials unset")

  end subroutine init_from_polyhedron

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

    integer :: i
    type(material_geometry) :: matl_init_geometry

    print '(a)', "initializing vof ... "
    call start_timer("vof init")

    call matl_init_geometry%init (plist, matl_ids)

    ! Dynamically schedule since only a few cells are mixed
    ! and they take the vast majority of the runtime.

    !$omp parallel do schedule(dynamic,100)
    do i = 1,mesh%ncell
      vof(:,i) = cell_vof(i, nmat, mesh, gmesh, matl_init_geometry)
    end do
    !$omp end parallel do

    call stop_timer("vof init")
    print '(a)', "done"

  end subroutine vof_initialize

  function cell_vof (i, nmat, mesh, gmesh, matl_init_geometry)

    use unstr_mesh_type
    use mesh_geom_type
    use parameter_list_type

    integer, intent(in) :: i, nmat
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    type(material_geometry), intent(in) :: matl_init_geometry
    real(r8) :: cell_vof(nmat)

    integer :: ierr
    type(dnc_cell) :: cell

    call cell%init (ierr, matl_init_geometry, i, mesh, gmesh, tesselate=.false.)
    cell_vof = cell%vof (matl_init_geometry, nmat, 0)

  end function cell_vof

  ! calculates the volume fractions of materials in a cell
  !
  !
  ! check if all vertices lie within a single material.
  ! if so, set the Vof for that material to 1.0.
  ! if not, divide the cell into tets and repeat recursively to a given threshold
  recursive function vof (this, matl_geometry, nmat, depth) result(hex_vof)

    use array_utils, only: isZero

    class(dnc_cell), intent(inout) :: this
    class(base_region), intent(in) :: matl_geometry
    integer,            intent(in) :: nmat,depth
    real(r8)                       :: hex_vof(nmat)

    real(r8) :: sum_vol
    integer :: i
    type(dnc_cell) :: subcell

    ! if the cell contains an interface (and therefore has at least two materials
    ! and we haven't yet hit our recursion limit, divide the hex and repeat
    if (depth < cell_vof_recursion_limit .and. this%contains_interface(matl_geometry)) then
      ! tally the vof from sub-cells
      call this%set_division_variables(subcell, matl_geometry)
      hex_vof = 0; sum_vol = 0
      do i = 1,this%ntets
        call this%get_subcell(subcell, matl_geometry, i)
        hex_vof = hex_vof + subcell%vof(matl_geometry,nmat,depth+1) * subcell%geom%volume()
        sum_vol = sum_vol + subcell%geom%volume()
      end do
      hex_vof = hex_vof / sum_vol
      ! if (.not.isZero(sum_vol - this%geom%volume())) then
      !   print *, sum_vol, this%geom%volume(), sum_vol - this%geom%volume()
      !   call LS_Fatal("volumes don't match")
      ! end if
    else
      ! if we are past the recursion limit
      ! or the cell does not contain an interface, calculate
      ! the vof in this hex based on the materials at its nodes
      hex_vof = this%vof_from_nodes(matl_geometry, nmat)
    end if

    ! if (any(ieee_is_nan(hex_vof)) .or. abs(sum(hex_vof) - 1) > 1e-10_r8) then
    !   print *, hex_vof, sum(hex_vof)
    !   print *, depth, sum_vol
    !   print *, this%geom%volume()
    !   this%geom%parent%vol = 0
    !   print *, this%geom%volume()
    !   subcell = this
    !   subcell%geom%parent%vol = 0
    !   call subcell%geom%tesselate()
    !   print *, subcell%geom%volume(), sum_vol * this%geom%volume() / subcell%geom%volume()
    !   print *, this%contains_interface(matl_geometry) .and. depth < cell_vof_recursion_limit
    !   print *, this%matl_at_node
    !   print *, hex_vof
    !   call LS_Fatal("invalid init vof")
    ! end if

  end function vof

  ! This function takes a hex, described by its vertex positions,
  ! and determines if it contains an interface for a given material.
  ! This is done by checking if every vertex lies within the material, or not.
  pure logical function contains_interface(this, matl_geometry)
    class(dnc_cell),     intent(in) :: this
    class(base_region), intent(in) :: matl_geometry
    contains_interface = any(this%matl_at_node /= this%matl_at_node(1))
  end function contains_interface

  ! calculate the volume fraction from nodes. used at the bottom depth
  !
  ! If we have >2 materials in this simulation, we need to calculate the vof
  ! from nodes only, using the average of the node materials. If we
  ! have 2 materials, we can reconstruct
  ! a plane from the signed distance function. Reconstructing the plane
  ! gives us much greater accuracy, which allows us to go to a lower depth
  ! and save significantly on computational cost.
  !
  ! TODO: We could switch this dynamically, so instead of looking at how many
  !       materials are in the entire simulation, we look at how many are in
  !       this cell. Then we can change the depth required for 2-phase cells,
  !       vs many-material cells.
  function vof_from_nodes (this, matl_geometry, nmat)

    use plane_type
    use array_utils, only: normalize
    use cell_geometry, only: cross_product

    class(dnc_cell), intent(inout) :: this
    class(base_region), intent(in) :: matl_geometry
    integer, intent(in) :: nmat
    real(r8) :: vof_from_nodes(nmat)

    integer :: m, ierr, i, e
    real(r8) :: x(3,3), s1, s2
    type(plane) :: P
    type(polyhedron) :: poly

    if (nmat > 2 .or. .not.this%contains_interface(matl_geometry)) then
      vof_from_nodes = 0
      do m = 1,nmat
        vof_from_nodes(m) = count(this%matl_at_node == m)
      end do
      vof_from_nodes = vof_from_nodes / this%geom%parent%nVerts
    else
      ! get 3 points intersecting the edges of this cell
      i = 0
      do e = 1,this%geom%parent%nEdges
        if (this%matl_at_node(this%geom%parent%edge_vid(1,e)) /= &
            this%matl_at_node(this%geom%parent%edge_vid(2,e))) then
          s1 = matl_geometry%signed_distance(this%geom%parent%x(:,this%geom%parent%edge_vid(1,e)))
          s2 = matl_geometry%signed_distance(this%geom%parent%x(:,this%geom%parent%edge_vid(2,e)))

          i = i + 1
          x(:,i) = this%geom%parent%x(:,this%geom%parent%edge_vid(1,e)) + &
              s1/(s1-s2) * (this%geom%parent%x(:,this%geom%parent%edge_vid(2,e)) - &
              this%geom%parent%x(:,this%geom%parent%edge_vid(1,e)))
          if (i == 3) exit
        end if
      end do
      ASSERT(e <= this%geom%parent%nEdges)
      if (e <= this%geom%parent%nEdges) then
        ! construct a plane from the 3 points
        P%normal = normalize(cross_product(x(:,2) - x(:,1), x(:,3) - x(:,1)))
        P%rho = dot_product(x(:,1), P%normal)

        ! find a node that is guaranteed not to intersect the interface
        ! flip the plane normal if it is outside the halfspace but inside the material
        do i = 1,this%geom%parent%nVerts
          s1 = P%signed_distance(this%geom%parent%x(:,i))

          if (abs(s1) > 0) then
            ! WARN: right now this is hardwired for material 2 being the "droplet"
            !       we will need to fix this
            if ((s1 > 0 .and. this%matl_at_node(i) == 2) .or. &
                (s1 < 0 .and. this%matl_at_node(i) /= 2)) then
              P%normal = -P%normal
              P%rho = -P%rho
            end if
            exit
          end if
        end do

        ! get the volume fraction behind the plane in this cell
        ! copy into new polyhedron structure because face normals
        ! aren't calculated at each level of the initialization
        ! maybe a better way of getting around this?
        ! WARN: the hardwiring for material 2 as the droplet is also here
        call poly%init(ierr, this%geom%parent%x)
        vof_from_nodes(2) = poly%volume_behind_plane(P, ierr) / this%geom%volume()
        vof_from_nodes(1) = 1 - vof_from_nodes(2)
      else
        ! interface doesn't actually intersect this subcell
        ! could happen if it just barely touches an edge
        vof_from_nodes = 0
        vof_from_nodes(this%matl_at_node(1)) = 1
      end if
    end if

  end function vof_from_nodes

  ! returns the center of the given hex
  pure function cell_center (this)
    class(dnc_cell), intent(in) :: this
    real(r8) :: cell_center(3)
    cell_center = sum(this%geom%parent%x, dim=2) / this%geom%parent%nVerts
  end function cell_center

  ! returns an array of face centers
  pure function face_centers(this)

    class(dnc_cell), intent(in) :: this
    real(r8) :: face_centers(3,this%geom%parent%nFaces)

    integer :: f, nV

    do f = 1,this%geom%parent%nFaces
      nV = count(this%geom%parent%face_vid(:,f) > 0)
      face_centers(:,f) = sum(this%geom%parent%x(:,this%geom%parent%face_vid(:nV,f)), dim=2) / nV
    end do

  end function face_centers

  ! returns an array of edge centers
  pure function edge_centers(this)

    class(dnc_cell), intent(in) :: this
    real(r8) :: edge_centers(3,this%geom%parent%nEdges)

    integer :: e

    do e = 1,this%geom%parent%nEdges
      edge_centers(:,e) = sum(this%geom%parent%x(:,this%geom%parent%edge_vid(:,e)), dim=2) / 2
    end do

  end function edge_centers

  subroutine set_division_variables (this, subcell, matl_geometry)

    use hex_types, only: tet_fv, tet_ev, tet_fe, tet_ef, tet_vf

    class(dnc_cell), intent(inout) :: this
    class(dnc_cell), intent(out) :: subcell
    class(base_region), intent(in) :: matl_geometry

    integer :: j

    ! find center points
    this%cellc = this%cell_center()
    this%facec = this%face_centers()
    this%edgec = this%edge_centers()

    ! get material ids at all points above
    if (allocated(this%mfacec)) deallocate(this%mfacec)
    if (allocated(this%medgec)) deallocate(this%medgec)
    allocate(this%mfacec(this%geom%parent%nFaces), this%medgec(this%geom%parent%nEdges))
    this%mcellc = matl_geometry%index_at(this%cellc)
    do j = 1,this%geom%parent%nFaces
      this%mfacec(j) = matl_geometry%index_at(this%facec(:,j))
    end do
    do j = 1,this%geom%parent%nEdges
      this%medgec(j) = matl_geometry%index_at(this%edgec(:,j))
    end do

    this%ntets = 2*count(this%geom%parent%face_vid > 0)

    ! initialize subcell structure (will be a tet)
    subcell%geom%parent%nVerts = 4
    subcell%geom%parent%nEdges = 6
    subcell%geom%parent%nFaces = 4
    subcell%geom%nchildren = 0
    subcell%geom%parent%edge_vid = tet_ev
    subcell%geom%parent%face_vid = tet_fv
    subcell%geom%parent%face_eid = tet_fe
    subcell%geom%parent%edge_faces = tet_ef
    subcell%geom%parent%vertex_faces = tet_vf

    allocate(subcell%geom%parent%x(3,4), subcell%matl_at_node(4))

  end subroutine set_division_variables

  subroutine get_subcell (this, subcell, matl_geometry, i)

    use cell_geometry, only: tet_volume

    class(dnc_cell), intent(in) :: this
    class(dnc_cell), intent(inout) :: subcell
    class(base_region), intent(in) :: matl_geometry
    integer, intent(in) :: i

    integer :: nV, f, v, feid, e, edge_side

    ! get face, face_edge_id, and edge_side of this tet
    ! assume all faces have the same number of nodes/edges
    ! i = (f-1)*2*nV + (feid-1)*2 + edge_side
    ! WARN: assumes all faces have the same number of vertices.
    !       good for bricks and tets, not prisms.
    nV = count(this%geom%parent%face_vid(:,1) > 0)

    f = (i-1) / (2*nV) + 1
    feid = (i - 1 - (f-1)*2*nV) / 2 + 1
    edge_side = i - (f-1)*2*nV - (feid-1)*2
    e = this%geom%parent%face_eid(feid,f)

    ! node geometry
    subcell%geom%parent%x(:,1) = this%cellc
    subcell%geom%parent%x(:,2) = this%facec(:,f)
    if (edge_side==1) then
      subcell%geom%parent%x(:,3) = this%geom%parent%x(:,this%geom%parent%face_vid(feid,f))
      subcell%geom%parent%x(:,4) = this%edgec(:,e)
    else
      subcell%geom%parent%x(:,3) = this%edgec(:,e)
      subcell%geom%parent%x(:,4) = &
          this%geom%parent%x(:,this%geom%parent%face_vid(modulo(feid,nV)+1,f))
    end if
    subcell%geom%parent%vol = tet_volume(subcell%geom%parent%x)

    ! materials at nodes
    subcell%matl_at_node(1) = this%mcellc
    subcell%matl_at_node(2) = this%mfacec(f)
    if (edge_side==1) then
      subcell%matl_at_node(3) = this%matl_at_node(this%geom%parent%face_vid(feid,f))
      subcell%matl_at_node(4) = this%medgec(e)
    else
      subcell%matl_at_node(3) = this%medgec(e)
      subcell%matl_at_node(4) = this%matl_at_node(this%geom%parent%face_vid(modulo(feid,nV)+1,f))
    end if

    !if (subcell%weight==0) call LS_fatal("zero weight")
    if (any(subcell%matl_at_node==0)) call LS_fatal("node materials unset")
    if (.not.allocated(subcell%geom%parent%x)) call LS_Fatal("invalid subcell")
    if (subcell%geom%volume()==0) then
      call this%geom%print_data()
      call LS_fatal("zero subcell volume")
    end if

  end subroutine get_subcell

end module vof_init
