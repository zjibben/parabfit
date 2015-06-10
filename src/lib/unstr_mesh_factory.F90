!!
!! UNSTR_MESH_FACTORY
!!
!! This module provides a generic procedure for instantiating an UNSTR_MESH
!! type object.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!! PROGRAMMING INTERFACE
!!
!!  The generic function NEW_UNSTR_MESH provides several different methods of
!!  instantiating an UNSTR_MESH object that is returned as the target of the
!!  pointer result of the function.
!!
!!  TYPE(UNSTR_MESH), POINTER :: MESH
!!
!!  MESH => NEW_UNSTR_MESH(PATH)
!!
!!    This form reads an Exodus/Genesis format mesh file.  The character
!!    argument PATH specifies the path of the file.  A relative path is
!!    interpreted relative to the execution directory.
!!
!!  MESH => NEW_UNSTR_MESH(XMIN, XMAX, NX [,EPS])
!!
!!    This form internally generates a uniform 3D hexahedral mesh. The meshed
!!    domain is a brick $(x_0,x_1)\times(y_0,y_1)\times(z_0,z_1)$ and is
!!    specified by giving the coordinates of two diagonally-opposed corners
!!    corners XMIN(:) = $(x_0,y_0,z_0)$ and XMAX(:) = $(x_1,y_1,z_1)$. The
!!    domain is diced up into equal-sized hexahedra by subdividing the domain
!!    into equal-sized intervals along each of the coordinate directions. The
!!    number of intervals in each direction is given by the integer array NX.
!!    If the optional scalar real argument EPS is given, the coordinates of
!!    the mesh nodes are perturbed by a random value in (-EPS,EPS).  Boundary
!!    nodes are not perturbed in directions normal to the boundary, so that
!!    the domain itself is not altered.  This is useful for testing purposes.
!!
!!  MESH => NEW_UNSTR_MESH(PARAMS)
!!
!!    This high-level form instantiates the mesh as specified by the parameter
!!    list PARAMS.  The parameter list should take one of the following forms
!!    that correspond to one of the preceding mesh instantiation methods.
!!
!!    Exodus/Genesis mesh: { "mesh-file": <string> }
!!
!!    Internally generated mesh:
!!    { "corner-min": <real-vector>, "corner-max": <real-vector>,
!!      "intervals":  <integer-vector>,
!!      "perturb-nodes": <real> // this is optional
!!    }
!!

#include "f90_assert.fpp"

module unstr_mesh_factory

  use unstr_mesh_type
  use logging_services
  !use hexahedral_mesh_support
  implicit none
  private

  public :: new_unstr_mesh

  interface new_unstr_mesh
    procedure new_unstr_mesh_exodus
    procedure new_unstr_mesh_regular
    procedure new_unstr_mesh_params
  end interface

contains

  !! Instantiate a new mesh object from an ExodusII mesh file
  function new_unstr_mesh_exodus (path) result (mesh)

    use exodus_file_type
    use permutations
    use unstr_mesh_tools, only: get_cell_neighbor_array, label_mesh_faces
    use cell_geometry

    character(*), intent(in) :: path
    type(unstr_mesh), pointer :: mesh

    type(exodus_file) :: file
    integer :: stat, ncface
    integer, allocatable :: cnhbr(:,:)
    character(:), allocatable :: errmsg

    !! Open the Exodus file and get some basic info about the mesh.
    call file%open (path, stat, errmsg)
    if (stat /= 0) call LS_fatal ('Unable to open Exodus file "' // path // '": ' // errmsg)

    allocate(mesh)

    !! Identify the type of mesh.
    mesh%mesh_type = mesh_type(file)
    if (mesh%mesh_type == 'UNKNOWN') &
        call LS_fatal ('mesh type not recognized; must be pure hex or pure tet')

    !! Read the mesh.
    call read_exodus_mesh (file, mesh, stat, errmsg)
    if (stat /= 0) call LS_fatal ('error reading Exodus mesh: ' // errmsg)


    !! At this point the mesh is completely defined mathematically (node
    !! coordinates and the cell-to-node connectivity information).  What
    !! remains is to derive additional data structures, like a labeling
    !! of mesh faces, good orderings of the mesh entities, and geometry info.

    !! Number the mesh faces
    select case (mesh%mesh_type)
    case ('TET')
      ncface = 4
    case ('HEX')
      ncface = 6
    case default
      INSIST(.false.) ! shouldn't be here
    end select

    !! Generate the array of cell face neighbors.  The topology of the mesh
    !! is examined in the process.  This cell-cell connectivity information
    !allocate(cnhbr, mold=mesh%cface)
    allocate(cnhbr(ncface,mesh%ncell))
    call get_cell_neighbor_array (mesh%cnode, cnhbr, stat)
    if (stat /= 0) call LS_fatal ('bad mesh topology detected')

    !! Generate a good ordering of the cells and reorder the existing
    !! cell-based mesh arrays.  The ordering will be based on the cell-cell
    !! adjacency graph; other choices are possible.
    allocate(mesh%xcell(mesh%ncell))
    call order_cells (cnhbr, mesh%xcell)
    call reorder (mesh%cnode,  mesh%xcell)
    call reorder (mesh%cblock, mesh%xcell)

    !! Generate a good ordering of the nodes and reorder the existing
    !! node-based arrays.  The ordering is tied to the cell ordering.
    allocate(mesh%xnode(mesh%nnode))
    call order_facets (mesh%cnode, mesh%xnode)
    call reorder (mesh%x, mesh%xnode)

    !! Generate a numbering of the mesh faces and the corresponding cell-to-face
    !! no need to generate a good ordering; the numbering process produces the
    !! same ordering as order_facets would.
    allocate(mesh%cface(ncface,mesh%ncell), mesh%cfpar(mesh%ncell))
    call label_mesh_faces (mesh%cnode, mesh%nface, mesh%cface, mesh%cfpar)

    !! Read side set data from the mesh file.
    !! we n
    call read_exodus_side_sets (file, mesh)
    call tag_boundary_faces (mesh%cface, mesh%face_set_mask)

    !! Initialize the mesh geometry components.
    call unstr_mesh_geom_init (mesh)

    call file%close

  contains

    !! Return the type of mesh described by the Exodus file.  Currently, only
    !! 3D meshes composed entirely of 4-node tetrahedral elements (return value
    !! 'TET') or entirely of 8-node hexahedral elements (return value 'HEX')
    !! are supported by the UNSTR_MESH data structure.  Anything else yields a
    !! return value 'UNKNOWN'.

    function mesh_type (file)
      use string_utilities, only: raise_case
      type(exodus_file), intent(in) :: file
      character(:), allocatable :: mesh_type, elem_type
      integer :: n
      mesh_type = 'UNKNOWN'
      select case (file%num_dim)
      case (3)
        elem_type = raise_case(file%elem_blk(1)%elem_type)
        select case (elem_type(1:3))
        case ('TET')
          !! Check that all blocks are 4-node tetrahedra
          do n = 1, file%num_elem_blk
            if (raise_case(file%elem_blk(n)%elem_type(1:3)) /= 'TET') return
            if (file%elem_blk(n)%num_nodes_per_elem /= 4) return
          end do
          mesh_type = 'TET'
        case ('HEX')
          !! Check that all blocks are 8-node hexahedra
          do n = 1, file%num_elem_blk
            if (raise_case(file%elem_blk(n)%elem_type(1:3)) /= 'HEX') return
            if (file%elem_blk(n)%num_nodes_per_elem /= 8) return
          end do
          mesh_type = 'HEX'
        end select
      end select
    end function mesh_type

    !! Reads the bare Exodus mesh, initializing the following components of the
    !! UNSTR_MESH data structure: NNODE, NCELL, X, CNODE, BLOCK_ID, and CBLOCK.
    !! It assumes the Exodus mesh is of homogeneous type, comprised entirely of
    !! elements of a single type.  Note: this is more general than is currently
    !! supported, but there is no loss in allowing the greater generality here.

    subroutine read_exodus_mesh (file, mesh, stat, errmsg)

      use exodus_file_type

      type(exodus_file), intent(in) :: file
      type(unstr_mesh), intent(inout) :: mesh
      integer, intent(out) :: stat
      character(:), allocatable, intent(out) :: errmsg

      integer :: n, nvert, offset

      mesh%ncell = file%num_elem
      mesh%nnode = file%num_node
      mesh%block_id = file%elem_blk%ID
      nvert = file%elem_blk(1)%num_nodes_per_elem ! assume homog mesh

      !! Read the connectivity data.
      allocate(mesh%cnode(nvert,mesh%ncell), mesh%cblock(mesh%ncell))
      offset = 0
      do n = 1, size(file%elem_blk)
        associate(connect => mesh%cnode(:,offset+1:offset+file%elem_blk(n)%num_elem), &
                  elem_blk_id => mesh%cblock(offset+1:offset+file%elem_blk(n)%num_elem))
          call file%get_connect (file%elem_blk(n)%ID, connect, stat, errmsg)
          if (stat /= 0) return
          elem_blk_id = file%elem_blk(n)%ID
        end associate
        offset = offset + file%elem_blk(n)%num_elem
      end do
      INSIST(offset == file%num_elem)

      !! Read the node coordinates.
      allocate(mesh%x(file%num_dim,mesh%nnode))
      call file%get_coord (mesh%x, stat, errmsg)

    end subroutine read_exodus_mesh

    !! Reads the side set data from the Exodus mesh file, translating it into
    !! the face set data stored in the mesh object.  Defines the FACE_SET_ID
    !! and FACE_SET_MASK mesh components.
    !! Requires that cface has been defined.

    subroutine read_exodus_side_sets (file, mesh)

      !use bitfield_type
      use permutations

      type(exodus_file), intent(in) :: file
      type(unstr_mesh), intent(inout) :: mesh

      integer :: j, k, n, error
      integer, allocatable :: elem_list(:), side_list(:)
      integer :: elem_map(mesh%ncell)

      !! The Exodus convention for numbering the sides of an element differs from
      !! that used by Truchas.  These arrays give the Exodus-to-Truchas mappings.
      integer, pointer :: side_map(:) => null()
      integer, target, save :: TET_SIDE_MAP(4), HEX_SIDE_MAP(6)
      data TET_SIDE_MAP/3,1,2,4/, HEX_SIDE_MAP/2,4,1,3,5,6/

      !! Link to the appropriate side mapping.
      select case (mesh%mesh_type)
      case ('TET')
        side_map => TET_SIDE_MAP
      case ('HEX')
        side_map => HEX_SIDE_MAP
      case default
        INSIST(.false.)
      end select

      !! The elements have been renumbered from their original Exodus
      !! numbering; generatate the old-to-new element number map.
      call invert_perm (mesh%xcell, elem_map)

      mesh%face_set_id = file%side_set%ID

      allocate(mesh%face_set_mask(mesh%nface))
      mesh%face_set_mask = 0 !ZERO_BITFIELD  ! unset all bits

      do k = 1, size(mesh%face_set_id)

        !! Read the side set data.
        n = file%side_set(k)%num_side_in_set
        allocate(elem_list(n), side_list(n))
        call file%get_side_set (file%side_set(k)%ID, elem_list, side_list, error)
        !error = ex_get_side_set(file%exoid, file%side_set(k)%ID, elem_list, side_list)
        INSIST(error == 0)  !TODO: replace with proper error handling

        !! Data sanity check (should we trust it?)
        INSIST(minval(elem_list) >= 1)
        INSIST(maxval(elem_list) <= size(mesh%cface,dim=2))
        INSIST(minval(side_list) >= 1)
        INSIST(maxval(side_list) <= size(mesh%cface,dim=1))

        !! Set bit K for all faces belonging to this side set.
        do j = 1, size(elem_list)
          n = mesh%cface(side_map(side_list(j)),elem_map(elem_list(j)))
          mesh%face_set_mask(n) = ibset(mesh%face_set_mask(n),pos=k)
        end do
        deallocate(elem_list, side_list)

      end do

    end subroutine read_exodus_side_sets

  end function new_unstr_mesh_exodus

  !! Generate a good ordering of the cells.  The cell-cell face
  !! adjacency info is given by CNHBR (0 values mark a boundary face),
  !! and the ordering is returned in the new-to-old mapping array
  !! PERM. 
  !!
  !! Correctness of this code and any mesh-partitioning code can be
  !! determined independently.  However, the performance of two will
  !! be determined by surface-to-volume ratios, and that ration will
  !! be determined by this code and the partitioning coded.

  subroutine order_cells (cnhbr, perm)
    integer, intent(in) :: cnhbr(:,:)
    integer, intent(out) :: perm(:)

    !! After some testing, RCM turns out to be poor compared to a
    !! simpler ordering.

    !call order_cells_RCM(cnhbr, perm)
    call order_cells_naive(cnhbr, perm)
  end subroutine order_cells
  
  subroutine order_cells_naive(cnhbr,perm)
    integer, intent(in) :: cnhbr(:,:)
    integer, intent(out) :: perm(:)
    integer :: j
    perm = [(j,j=1,size(perm),1)]
  end subroutine order_cells_naive

  subroutine order_cells_RCM (cnhbr, perm)
    use permutations
    integer, intent(in) :: cnhbr(:,:)
    integer, intent(out) :: perm(:)
    integer :: scratch(size(cnhbr,dim=1))
    logical :: ordered(size(cnhbr,dim=2))
    integer :: j,k,t,min_deg,min_idx,num_ordered,order_index
    ASSERT(size(perm) == size(cnhbr,dim=2))
    ! Nothing is ordered at the beginning
    num_ordered = 0
    ordered = .false.
    ! Find a vertex with minimal degree -- the first vertex
    min_idx = 1
    min_deg = vert_degree(min_idx)
    do j=2,size(cnhbr,dim=2)
       if (vert_degree(j) < min_deg) then
          min_idx = j
          min_deg = vert_degree(min_idx)
       end if
    end do
    ! Add the first vertex index
    call add_index(min_idx)
    ! Perform the RCM iteration
    order_index = 0
    do while (num_ordered .lt. size(cnhbr,dim=2))
       order_index = order_index+1
       ! Sort the neighbors of the currently ordered cell.
       scratch=cnhbr(:,perm(order_index))
       do j=1,size(scratch)-1
          do k=j+1,size(scratch)
             if (vert_degree(scratch(j)) > vert_degree(scratch(k))) then
                t = scratch(j)
                scratch(j) = scratch(k)
                scratch(k) = t
             end if
          end do
       end do ! Sorting scratch
       do j=1,size(scratch) 
          call add_index(scratch(j))
       end do ! Adding elements of scratch
    end do ! Ordering points
    ASSERT(is_perm(perm))
  contains
    function vert_degree(n) result(d)
      integer :: d
      integer, intent(in) :: n
      if (n >= 1 .and. n <= size(cnhbr,dim=2)) then      
         d = count(cnhbr(:,n).gt.0)
      else
         ! Handle the case that we are asked for the degree of a
         ! non-cell by returning a degree that is higher than the
         ! degree of any actual cell.
         d = 1 + size(cnhbr,dim=2) 
      end if
    end function vert_degree
    subroutine add_index(n)
      integer, intent(in) :: n
      ! Silently ignore requests to add indicies for non-cells
      if (n >= 1 .and. n <= size(cnhbr,dim=2)) then
         if (.not. ordered(n)) then
            num_ordered = num_ordered+1
            perm(num_ordered) = n
            ordered(n) = .true.
         end if
      end if
    end subroutine add_index
  end subroutine order_cells_RCM

  !! This procedure generates a new numbering of mesh facets (nodes, edges, or
  !! faces) that is well-ordered.  FACET is a cell-based array where FACET(:,j)
  !! are the facets of one type belonging to cell j.  The facets are renumbered
  !! according to the order they are encountered in the FACET array, and assumes
  !! the cells are already well-ordered.  The numbering permutation (new-to-old)
  !! is returned in the array PERM.  The values of the FACET array are mapped to
  !! the new values.

  subroutine order_facets (facet, perm)

    use permutations

    integer, intent(inout) :: facet(:,:)
    integer, intent(out)   :: perm(:)

    integer :: j, k, n

    ASSERT(minval(facet) == 1 .and. maxval(facet) == size(perm))

    !! Number the facets consecutively as they are encountered
    !! in the FACET array; PERM is the old-to-new numering.
    n = 0
    perm = 0
    do j = 1, size(facet,dim=2)
      do k = 1, size(facet,dim=1)
        if (perm(facet(k,j)) /= 0) cycle  ! already numbered
        n = n + 1
        perm(facet(k,j)) = n
      end do
    end do
    ASSERT(is_perm(perm))

    !! Map the FACET array values; PERM is old-to-new.
    do j = 1, size(facet,2)
      do k = 1, size(facet,1)
        facet(k,j) = perm(facet(k,j))
      end do
    end do

    !! Return the generated facet ordering (new-to-old).
    call invert_perm (perm)

  end subroutine order_facets

  !! This routine tags boundary faces in the FACE_SET_MASK array.  Bit 0 is
  !! set if the face is on the boundary; otherwise it is unset.  None of the
  !! other bits, which correspond to face sets, are modified.

  subroutine tag_boundary_faces (cface, face_set_mask)
    integer, intent(in) :: cface(:,:)
    integer, intent(inout) :: face_set_mask(:)
    integer :: k, j, tag(size(face_set_mask))
    tag = 0
    do j = 1, size(cface,dim=2)
      do k = 1, size(cface,dim=1)
        tag(cface(k,j)) = tag(cface(k,j)) + 1
      end do
    end do
    INSIST(all(tag >= 1 .and. tag <= 2))
    do j = 1, size(tag)
      if (tag(j) == 1) then
        face_set_mask(j) = ibset(face_set_mask(j),pos=0)
      else
        face_set_mask(j) = ibclr(face_set_mask(j),pos=0)
      end if
    end do
  end subroutine tag_boundary_faces

  !! This routine computes the geometrical quantities and secondary
  !! topological connectivity that is derivable from the primitive mesh.

  subroutine unstr_mesh_geom_init (mesh)

    use unstr_mesh_tools, only: get_face_node_array
    use cell_geometry

    type(unstr_mesh), intent(inout) :: mesh

    integer :: j

    select case (mesh%mesh_type)
    case ('TET')

      allocate(mesh%fnode(3,mesh%nface), mesh%fcell(2,mesh%nface))
      call get_face_node_array (mesh%cnode, mesh%cface, mesh%cfpar, mesh%fnode, mesh%fcell)
      ASSERT(minval(mesh%fcell(1,:)) > 0)
      ASSERT(all(mesh%fcell(2,:) == 0 .eqv. btest(mesh%face_set_mask,pos=0)))

      !! Directed face areas and face areas.
      allocate(mesh%normal(3,mesh%nface), mesh%area(mesh%nface))
      do j = 1, mesh%nface
        mesh%normal(:,j) = tri_face_normal(mesh%x(:,mesh%fnode(:,j)))
        mesh%area(j) = vector_length(mesh%normal(:,j))
      end do

      allocate(mesh%volume(mesh%ncell))
      do j = 1, mesh%ncell
        mesh%volume(j) = tet_volume(mesh%x(:,mesh%cnode(:,j)))
      end do

    case ('HEX')

      allocate(mesh%fnode(4,mesh%nface), mesh%fcell(2,mesh%nface))
      call get_face_node_array (mesh%cnode, mesh%cface, mesh%cfpar, mesh%fnode, mesh%fcell)
      ASSERT(minval(mesh%fcell(1,:)) > 0)
      ASSERT(all(mesh%fcell(2,:) == 0 .eqv. btest(mesh%face_set_mask,pos=0)))

      !! Directed face areas and face areas.
      allocate(mesh%normal(3,mesh%nface), mesh%area(mesh%nface))
      do j = 1, mesh%nface
        mesh%normal(:,j) = quad_face_normal(mesh%x(:,mesh%fnode(:,j)))
        mesh%area(j) = vector_length(mesh%normal(:,j))
      end do

      allocate(mesh%volume(mesh%ncell), mesh%corner_volume(8,mesh%ncell))
      do j = 1, mesh%ncell
        call eval_hex_volumes (mesh%x(:,mesh%cnode(:,j)), mesh%volume(j), mesh%corner_volume(:,j))
      end do

    case default
      INSIST(.false.) ! should not be here
    end select

  end subroutine unstr_mesh_geom_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function new_unstr_mesh_regular (xmin, xmax, nx, eps) result (mesh)

    use kinds, only: r8
    use unstr_mesh_tools, only: label_mesh_faces

    real(r8), intent(in) :: xmin(:), xmax(:)
    integer,  intent(in) :: nx(:)
    real(r8), intent(in), optional :: eps
    type(unstr_mesh), pointer :: mesh

    integer  :: i, j, k, n
    real(r8) :: x(3)

    ASSERT(size(xmin) == 3)
    ASSERT(size(xmax) == 3)
    ASSERT(size(nx) == 3)
    ASSERT(all(xmin < xmax))
    ASSERT(all(nx > 0))

    allocate(mesh)

    mesh%mesh_type = 'HEX'
    mesh%ncell = product(nx)
    mesh%nnode = product(nx+1)

    !! Generate the regular connectivity
    allocate(mesh%cnode(8,mesh%ncell))
    n = 0
    do k = 1, nx(3)
      do j = 1, nx(2)
        do i = 1, nx(1)
          n = cell_index(i,j,k)
          mesh%cnode(1,n) = node_index(i,j,k)
          mesh%cnode(2,n) = node_index(i+1,j,k)
          mesh%cnode(3,n) = node_index(i+1,j+1,k)
          mesh%cnode(4,n) = node_index(i,j+1,k)
          mesh%cnode(5,n) = node_index(i,j,k+1)
          mesh%cnode(6,n) = node_index(i+1,j,k+1)
          mesh%cnode(7,n) = node_index(i+1,j+1,k+1)
          mesh%cnode(8,n) = node_index(i,j+1,k+1)
        end do
      end do
    end do

    !! All cells go into a single cell block
    mesh%block_id = [1]
    allocate(mesh%cblock(mesh%ncell))
    mesh%cblock = 1

    !! Retain the cell numbering.
    allocate(mesh%xcell(mesh%ncell))
    do j = 1, mesh%ncell
      mesh%xcell(j) = j
    end do

    !! Generate the node coordinates.
    allocate(mesh%x(3,mesh%nnode))
    do k = 1, nx(3)+1
      x(3) = ((nx(3)-k+1)/real(nx(3),r8))*xmin(3) + ((k-1)/real(nx(3),r8))*xmax(3)
      do j = 1, nx(2)+1
        x(2) = ((nx(2)-j+1)/real(nx(2),r8))*xmin(2) + ((j-1)/real(nx(2),r8))*xmax(2)
        do i = 1, nx(1)+1
          x(1) = ((nx(1)-i+1)/real(nx(1),r8))*xmin(1) + ((i-1)/real(nx(1),r8))*xmax(1)
          mesh%x(:,node_index(i,j,k)) = x
        end do
      end do
    end do

    if (present(eps)) call randomize_coord (eps)

    !! Retain the node numbering.
    allocate(mesh%xnode(mesh%nnode))
    do j = 1, mesh%nnode
      mesh%xnode(j) = j
    end do

    !! Generate a numbering of the mesh faces and the cell-to-face map.
    !! Faces are numbered as encountered when interating through the cells.
    allocate(mesh%cface(6,mesh%ncell), mesh%cfpar(mesh%ncell))
    call label_mesh_faces (mesh%cnode, mesh%nface, mesh%cface, mesh%cfpar)

    !! Generate side set data: 6 sides of the brick, IDs 1 through 6.
    mesh%face_set_id = [1, 2, 3, 4, 5, 6]
    allocate(mesh%face_set_mask(mesh%nface))
    mesh%face_set_mask = 0 !ZERO_BITFIELD  ! unset all bits
    call tag_boundary_faces (mesh%cface, mesh%face_set_mask)
    !! Side sets 1 (x=xmin) and 2 (x=xmax)
    do k = 1, nx(3)
      do j = 1, nx(2)
        n = mesh%cface(3,cell_index(1,j,k))
        mesh%face_set_mask(n) = ibset(mesh%face_set_mask(n),pos=1)
        n = mesh%cface(4,cell_index(nx(1),j,k))
        mesh%face_set_mask(n) = ibset(mesh%face_set_mask(n),pos=2)
      end do
    end do
    !! Side sets 3 (y=ymin) and 4 (y=ymax)
    do k = 1, nx(3)
      do i = 1, nx(1)
        n = mesh%cface(2,cell_index(i,1,k))
        mesh%face_set_mask(n) = ibset(mesh%face_set_mask(n),pos=3)
        n = mesh%cface(1,cell_index(i,nx(2),k))
        mesh%face_set_mask(n) = ibset(mesh%face_set_mask(n),pos=4)
      end do
    end do
    !! Side sets 5 and 6 (z=zmin) and (z=zmax)
    do j = 1, nx(2)
      do i = 1, nx(1)
        n = mesh%cface(5,cell_index(i,j,1))
        mesh%face_set_mask(n) = ibset(mesh%face_set_mask(n),pos=5)
        n = mesh%cface(6,cell_index(i,j,nx(3)))
        mesh%face_set_mask(n) = ibset(mesh%face_set_mask(n),pos=6)
      end do
    end do

    !! Initialize the mesh geometry components.
    call unstr_mesh_geom_init (mesh)

  contains

    integer function node_index (i, j, k)
      integer, intent(in) :: i, j, k
      node_index = i + ((j-1) + (k-1)*(nx(2)+1))*(nx(1)+1)
    end function node_index

    integer function cell_index (i, j, k)
      integer, intent(in) :: i, j, k
      cell_index = i + ((j-1) + (k-1)*nx(2))*nx(1)
    end function cell_index

    subroutine randomize_coord (eps)
      real(r8), intent(in) :: eps
      integer :: i, j, k, n
      logical :: mask(3)
      real(r8) :: dx(3)
      ASSERT(eps >= 0.0_r8)
      if (eps == 0.0_r8) return
      do k = 1, nx(3)
        mask(3) = (k > 1 .and. k < nx(3))
        do j = 1, nx(2)
          mask(2) = (j > 1 .and. j < nx(2))
          do i = 1, nx(1)
            mask(1) = (i > 1 .and. i < nx(1))
            call random_number (dx) ! in [0,1)
            dx = eps*(2*dx - 1) ! in [-eps,eps)
            n = node_index(i,j,k)
            mesh%x(:,n) = mesh%x(:,n) + merge(dx, 0.0_r8, mask)
          end do
        end do
      end do
    end subroutine randomize_coord

  end function new_unstr_mesh_regular

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function new_unstr_mesh_params (params) result (mesh)

    use kinds, only: r8
    use parameter_list_type
    use logging_services

    type(parameter_list) :: params
    type(unstr_mesh), pointer :: mesh

    integer :: stat
    character(:), allocatable :: errmsg, mesh_file, context
    real(r8), allocatable :: xmin(:), xmax(:)
    integer,  allocatable :: nx(:)
    real(r8) :: eps

    context = 'processing ' // params%name() // ': '

    if (params%is_parameter('mesh-file')) then
      !! Read an Exodus/Genesis mesh file.
      call params%get('mesh-file', mesh_file, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      mesh => new_unstr_mesh_exodus(mesh_file)
    else
      !! Internally-generated regular mesh.
      call params%get('corner-min', xmin, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (size(xmin) /= 3) call LS_fatal (context//'3-vector value expected for "corner-min"')
      call params%get('corner-max', xmax, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (size(xmax) /= 3) call LS_fatal (context//'3-vector value expected for "corner-max"')
      if (any(xmax <= xmin)) call LS_fatal (context//'require "corner-min" < "corner-max"')
      call params%get('intervals', nx, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (size(nx) /= 3) call LS_fatal (context//'3-vector value expected for "intervals"')
      if (any(nx < 1)) call LS_fatal (context//'"intervals" must be > 0')
      call params%get('perturb-nodes', eps, default=0.0_r8, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (eps < 0.0_r8) call LS_fatal (context//'"perturb-nodes" must be >= 0.0')
      mesh => new_unstr_mesh_regular(xmin, xmax, nx, eps)
    end if

  end function new_unstr_mesh_params

end module unstr_mesh_factory
