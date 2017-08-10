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
  use hex_types
  implicit none
  private

  ! hex type for the divide and conquer algorithm
  ! public to expose for testing
  type, extends(base_hex), public :: dnc_hex
    private
    integer :: matl_at_node(8)
  contains
    procedure          :: divide
    procedure          :: cell_center
    procedure          :: face_centers
    procedure          :: edge_centers
    procedure          :: contains_interface
    procedure          :: vof
    procedure, private :: vof_from_nodes
  end type dnc_hex

  public :: vof_initialize

  ! TODO: make these user-specified parameters
  integer , parameter :: cell_vof_recursion_limit = 12

contains

  ! this subroutine initializes the vof values in all cells
  ! from a user input function. It calls a divide and conquer algorithm
  ! to calculate the volume fractions given an interface
  subroutine vof_initialize (mesh, plist, vof, matl_ids, nmat) !cell_matls)

    use unstr_mesh_type
    use parameter_list_type
    use timer_tree_type

    type(unstr_mesh),     intent(in)    :: mesh
    type(parameter_list), intent(in)    :: plist
    integer,              intent(in)    :: nmat
    integer,              intent(in)    :: matl_ids(:)
    real(r8),             intent(inout) :: vof(:,:)

    integer                 :: i,v
    type(dnc_hex)           :: hex
    type(material_geometry) :: matl_init_geometry

    call start_timer("vof init")

    ! first, determine the initial state provided by the user, and assign
    ! the function which will determine what points are inside the materials
    call matl_init_geometry%init (plist, matl_ids)

    ! next, loop though every cell and check if all vertices lie within a single material
    ! if so, set the Vof for that material to 1.0.
    ! if not, divide the hex cell into 8 subdomains defined by the centroid and centers of faces
    ! repeat the process recursively for each subdomain up to a given threshold

    !$omp parallel do default(private) shared(vof,mesh,matl_init_geometry,nmat)
    do i = 1,mesh%ncell
      ! initialize dnc_hex
      hex%node = mesh%x(:,mesh%cnode(:,i))
      do v = 1,8
        hex%matl_at_node(v) = matl_init_geometry%index_at(hex%node(:,v))
      end do

      ! calculate the vof
      vof(:,i) = hex%vof (matl_init_geometry, nmat, 0)
    end do
    !$omp end parallel do

    call stop_timer("vof init")

  end subroutine vof_initialize

  ! calculates the volume fractions of materials in a cell
  recursive function vof (this, matl_geometry, nmat, depth) result(hex_vof)

    class(dnc_hex),     intent(in) :: this
    class(base_region), intent(in) :: matl_geometry
    integer,            intent(in) :: nmat,depth
    real(r8)                       :: hex_vof(nmat)

    integer       :: i
    type(dnc_hex) :: subhex(8)
    real(r8)      :: this_volume

    !this_volume = this%calc_volume()

    ! if the cell contains an interface (and therefore has at least two materials
    ! and we haven't yet hit our recursion limit, divide the hex and repeat
    if (this%contains_interface(matl_geometry) .and. depth < cell_vof_recursion_limit) then
      ! divide into 8 smaller hexes
      subhex = this%divide(matl_geometry)

      ! tally the vof from subhexes
      ! WARN: For nonorthogonal meshes, the subhex volumes are not all equal,
      !       so we must calculate the ratio of volumes.
      hex_vof = 0
      do i = 1,8
        hex_vof = hex_vof + subhex(i)%vof(matl_geometry,nmat,depth+1) !* &
            !subhex(i)%calc_volume() / this_volume
      end do
      hex_vof = hex_vof / 8

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
    cell_center = sum(hex%node, dim=2) / 8
  end function cell_center

  ! returns an array of face centers
  pure function face_centers(this)

    use consts,    only: nfc
    use hex_types, only: hex_f

    class(dnc_hex), intent(in) :: this
    real(r8) :: face_centers(3,6)

    integer :: f

    do f = 1,nfc
      face_centers(:,f) = sum(this%node(:,hex_f(:,f)), dim=2) / 4
    end do

  end function face_centers

  ! returns an array of edge centers
  pure function edge_centers(this)

    use hex_types, only: hex_e

    class(dnc_hex), intent(in) :: this
    real(r8) :: edge_centers(3,12)

    integer :: e

    do e = 1,12
      !edge_centers(:,e) = sum(this%node(:,hex_e(:,e)), dim=2) / 2
      edge_centers(:,e) = (this%node(:,hex_e(1,e)) + this%node(:,hex_e(2,e))) / 2
    end do

  end function edge_centers

  ! this function takes a hex and returns an array of 8 hexes obtained from dividing the input
  function divide(this, matl_geometry) result(subhex)

    class(dnc_hex), intent(in) :: this
    class(base_region), intent(in) :: matl_geometry
    type(dnc_hex) :: subhex(8)

    real(r8) :: cellc(3), facec(3,6), edgec(3,12)
    integer :: i,j, mcellc, mfacec(6), medgec(12)

    ! do i = 1,8
    !   do j = 1,8
    !     subhex(i)%node(:,j) = 0.5_r8 * (this%node(:,i) + this%node(:,j))
    !   end do
    ! end do

    ! find center points
    cellc = this%cell_center()
    facec = this%face_centers()
    edgec = this%edge_centers()

    ! get material ids at all points above
    mcellc = matl_geometry%index_at(cellc)
    do i = 1,6
      mfacec(i) = matl_geometry%index_at(facec(:,i))
    end do
    do i = 1,12
      medgec(i) = matl_geometry%index_at(edgec(:,i))
    end do

    ! set subhex node positions
    subhex(1)%node(:,1) = this%node(:,1)
    subhex(1)%node(:,2) = edgec(:,1)
    subhex(1)%node(:,3) = facec(:,5)
    subhex(1)%node(:,4) = edgec(:,4)
    subhex(1)%node(:,5) = edgec(:,5)
    subhex(1)%node(:,6) = facec(:,2)
    subhex(1)%node(:,7) = cellc(:)
    subhex(1)%node(:,8) = facec(:,3)

    subhex(1)%matl_at_node(1) = this%matl_at_node(1)
    subhex(1)%matl_at_node(2) = medgec(1)
    subhex(1)%matl_at_node(3) = mfacec(5)
    subhex(1)%matl_at_node(4) = medgec(4)
    subhex(1)%matl_at_node(5) = medgec(5)
    subhex(1)%matl_at_node(6) = mfacec(2)
    subhex(1)%matl_at_node(7) = mcellc
    subhex(1)%matl_at_node(8) = mfacec(3)

    subhex(2)%node(:,1) = edgec(:,1)
    subhex(2)%node(:,2) = this%node(:,2)
    subhex(2)%node(:,3) = edgec(:,2)
    subhex(2)%node(:,4) = facec(:,5)
    subhex(2)%node(:,5) = facec(:,2)
    subhex(2)%node(:,6) = edgec(:,6)
    subhex(2)%node(:,7) = facec(:,4)
    subhex(2)%node(:,8) = cellc(:)

    subhex(2)%matl_at_node(1) = medgec(1)
    subhex(2)%matl_at_node(2) = this%matl_at_node(2)
    subhex(2)%matl_at_node(3) = medgec(2)
    subhex(2)%matl_at_node(4) = mfacec(5)
    subhex(2)%matl_at_node(5) = mfacec(2)
    subhex(2)%matl_at_node(6) = medgec(6)
    subhex(2)%matl_at_node(7) = mfacec(4)
    subhex(2)%matl_at_node(8) = mcellc

    subhex(3)%node(:,1) = facec(:,5)
    subhex(3)%node(:,2) = edgec(:,2)
    subhex(3)%node(:,3) = this%node(:,3)
    subhex(3)%node(:,4) = edgec(:,3)
    subhex(3)%node(:,5) = cellc(:)
    subhex(3)%node(:,6) = facec(:,4)
    subhex(3)%node(:,7) = edgec(:,7)
    subhex(3)%node(:,8) = facec(:,1)

    subhex(3)%matl_at_node(1) = mfacec(5)
    subhex(3)%matl_at_node(2) = medgec(2)
    subhex(3)%matl_at_node(3) = this%matl_at_node(3)
    subhex(3)%matl_at_node(4) = medgec(3)
    subhex(3)%matl_at_node(5) = mcellc
    subhex(3)%matl_at_node(6) = mfacec(4)
    subhex(3)%matl_at_node(7) = medgec(7)
    subhex(3)%matl_at_node(8) = mfacec(1)

    subhex(4)%node(:,1) = edgec(:,4)
    subhex(4)%node(:,2) = facec(:,5)
    subhex(4)%node(:,3) = edgec(:,3)
    subhex(4)%node(:,4) = this%node(:,4)
    subhex(4)%node(:,5) = facec(:,3)
    subhex(4)%node(:,6) = cellc(:)
    subhex(4)%node(:,7) = facec(:,1)
    subhex(4)%node(:,8) = edgec(:,8)

    subhex(4)%matl_at_node(1) = medgec(4)
    subhex(4)%matl_at_node(2) = mfacec(5)
    subhex(4)%matl_at_node(3) = medgec(3)
    subhex(4)%matl_at_node(4) = this%matl_at_node(4)
    subhex(4)%matl_at_node(5) = mfacec(3)
    subhex(4)%matl_at_node(6) = mcellc
    subhex(4)%matl_at_node(7) = mfacec(1)
    subhex(4)%matl_at_node(8) = medgec(8)

    subhex(5)%node(:,1) = edgec(:,5)
    subhex(5)%node(:,2) = facec(:,2)
    subhex(5)%node(:,3) = cellc(:)
    subhex(5)%node(:,4) = facec(:,3)
    subhex(5)%node(:,5) = this%node(:,5)
    subhex(5)%node(:,6) = edgec(:,9)
    subhex(5)%node(:,7) = facec(:,6)
    subhex(5)%node(:,8) = edgec(:,12)

    subhex(5)%matl_at_node(1) = medgec(5)
    subhex(5)%matl_at_node(2) = mfacec(2)
    subhex(5)%matl_at_node(3) = mcellc
    subhex(5)%matl_at_node(4) = mfacec(3)
    subhex(5)%matl_at_node(5) = this%matl_at_node(5)
    subhex(5)%matl_at_node(6) = medgec(9)
    subhex(5)%matl_at_node(7) = mfacec(6)
    subhex(5)%matl_at_node(8) = medgec(12)

    subhex(6)%node(:,1) = facec(:,2)
    subhex(6)%node(:,2) = edgec(:,6)
    subhex(6)%node(:,3) = facec(:,4)
    subhex(6)%node(:,4) = cellc(:)
    subhex(6)%node(:,5) = edgec(:,9)
    subhex(6)%node(:,6) = this%node(:,6)
    subhex(6)%node(:,7) = edgec(:,10)
    subhex(6)%node(:,8) = facec(:,6)

    subhex(6)%matl_at_node(1) = mfacec(2)
    subhex(6)%matl_at_node(2) = medgec(6)
    subhex(6)%matl_at_node(3) = mfacec(4)
    subhex(6)%matl_at_node(4) = mcellc
    subhex(6)%matl_at_node(5) = medgec(9)
    subhex(6)%matl_at_node(6) = this%matl_at_node(6)
    subhex(6)%matl_at_node(7) = medgec(10)
    subhex(6)%matl_at_node(8) = mfacec(6)

    subhex(7)%node(:,1) = cellc(:)
    subhex(7)%node(:,2) = facec(:,4)
    subhex(7)%node(:,3) = edgec(:,7)
    subhex(7)%node(:,4) = facec(:,1)
    subhex(7)%node(:,5) = facec(:,6)
    subhex(7)%node(:,6) = edgec(:,10)
    subhex(7)%node(:,7) = this%node(:,7)
    subhex(7)%node(:,8) = edgec(:,11)

    subhex(7)%matl_at_node(1) = mcellc
    subhex(7)%matl_at_node(2) = mfacec(4)
    subhex(7)%matl_at_node(3) = medgec(7)
    subhex(7)%matl_at_node(4) = mfacec(1)
    subhex(7)%matl_at_node(5) = mfacec(6)
    subhex(7)%matl_at_node(6) = medgec(10)
    subhex(7)%matl_at_node(7) = this%matl_at_node(7)
    subhex(7)%matl_at_node(8) = medgec(11)

    subhex(8)%node(:,1) = facec(:,3)
    subhex(8)%node(:,2) = cellc(:)
    subhex(8)%node(:,3) = facec(:,1)
    subhex(8)%node(:,4) = edgec(:,8)
    subhex(8)%node(:,5) = edgec(:,12)
    subhex(8)%node(:,6) = facec(:,6)
    subhex(8)%node(:,7) = edgec(:,11)
    subhex(8)%node(:,8) = this%node(:,8)

    subhex(8)%matl_at_node(1) = mfacec(3)
    subhex(8)%matl_at_node(2) = mcellc
    subhex(8)%matl_at_node(3) = mfacec(1)
    subhex(8)%matl_at_node(4) = medgec(8)
    subhex(8)%matl_at_node(5) = medgec(12)
    subhex(8)%matl_at_node(6) = mfacec(6)
    subhex(8)%matl_at_node(7) = medgec(11)
    subhex(8)%matl_at_node(8) = this%matl_at_node(8)

  end function divide

end module vof_init
