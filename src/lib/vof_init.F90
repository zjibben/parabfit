!!
!!
!!
!!
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

module vof_init
  use kinds, only: r8
  use material_geometry_type
  use logging_services
  use hex_types
  implicit none
  private

  ! hex type to make divide and conquer algorithm simpler
  type, extends(base_hex) :: dnc_hex
  contains
    procedure         :: divide
    procedure         :: cell_center
    procedure         :: face_centers
    procedure         :: edge_centers
    procedure         :: contains_interface
    procedure         :: vof
    procedure,private :: vof_from_nodes
  end type dnc_hex
  
  public :: vof_initialize

  integer , parameter :: cell_vof_recursion_limit  = 100     ! should be made user-specified parameter
  real(r8), parameter :: cell_vof_tolerance_factor = 1e-3_r8 ! should be made user-specified parameter
contains

  subroutine vof_initialize (mesh, plist, vof, matl_ids, nmat) !cell_matls)
    use unstr_mesh_type

    ! this subroutine initializes the vof values in all cells
    ! from a user input function. It calls a divide and conquer algorithm
    ! to calculate the volume fractions given an interface
    type(unstr_mesh), pointer, intent(in) :: mesh
    type(parameter_list), intent(in)      :: plist
    integer, intent(in) :: nmat
    integer, intent(in) :: matl_ids(:)
    real(r8), intent(inout) :: vof(:,:)

    ! local variables
    integer :: i
    real(r8) :: tolerance
    type(dnc_hex) :: hex
    type(material_geometry) :: matl_init_geometry

    ! first, determine the initial state provided by the user, and assign
    ! the function which will determine what points are inside the materials
    call matl_init_geometry%init(plist, matl_ids)

    tolerance = cell_vof_tolerance_factor * minval(mesh%volume(:))

    ! next, loop though every cell and check if all verteces lie within a single material
    ! if so, set the Vof for that material to 1.0.
    ! if not, divide the hexahedral cell into 8 subdomains defined by the centroid and centers of faces
    ! repeat the process recursively for each subdomain up to a given threshold

    !$omp parallel do default(private) shared(vof,mesh,matl_init_geometry,nmat,tolerance)
    do i = 1,mesh%ncell
      ! make a hex type out of the cell before calculating the vof
      ! this will make it easier to handle
      hex%node = mesh%x(:,mesh%cnode(:,i))

      ! calculate the vof
      vof(:,i) = hex%vof(matl_init_geometry, nmat, 0, tolerance)
    end do
    !$omp end parallel do

  end subroutine vof_initialize

  recursive function vof(this, matl_geometry, nmat, depth, tolerance) result(hex_vof)
    ! calculates the volume fractions of materials in a cell
    class(dnc_hex), intent(in) :: this
    class(base_region), intent(in) :: matl_geometry
    integer , intent(in) :: nmat
    integer , intent(in) :: depth
    real(r8), intent(in) :: tolerance
    real(r8) :: hex_vof(nmat)

    integer       :: i
    type(dnc_hex) :: subhex(8)
    real(r8)      :: this_volume

    this_volume = this%calc_volume()

    ! if the cell contains an interface (and therefore has at least two materials
    ! and we haven't yet hit our recursion limit, divide the hex and repeat
    if (this%contains_interface(matl_geometry) .and. depth < cell_vof_recursion_limit .and. this_volume > tolerance) then
      !if (this%depth <= cell_vof_recursion_limit) then
      ! divide into 8 smaller hexes
      subhex = this%divide()

      ! tally the vof from subhexes
      hex_vof = 0.0_r8
      do i = 1,8
        hex_vof = hex_vof + subhex(i)%vof(matl_geometry,nmat,depth+1,tolerance) * subhex(i)%calc_volume()/this_volume
      end do

    else
      ! if we are at (or somehow past) the recursion limit or the cell does not contain an interface,
      ! calculate the vof in this hex based on the materials at its nodes
      hex_vof = this%vof_from_nodes(matl_geometry, nmat)
    end if

  end function vof

  logical function contains_interface(this, matl_geometry)
    ! This function takes a hex, described by its vertex positions,
    ! and determines if it contains an interface for a given material.
    ! This is done by checking if every vertex lies within the material, or not.
    class(dnc_hex), intent(in) :: this
    class(base_region), intent(in) :: matl_geometry

    integer :: i,reference_id
    real(r8) :: facec(3,6),edgec(3,12)

    contains_interface = .false.
    reference_id = matl_geometry%index_at(this%node(:,1))

    ! if any node doesn't have the same material as the first node, there's an interface in this cell
    do i = 2,8
      contains_interface = contains_interface .or. reference_id /= matl_geometry%index_at(this%node(:,i))
    end do

  end function contains_interface

  function vof_from_nodes(this, matl_geometry, nmat)
    ! the material at each node contributes 1/8th of the vof in the given hex
    class(dnc_hex), intent(in) :: this
    class(base_region), intent(in) :: matl_geometry
    integer, intent(in) :: nmat
    real(r8) :: vof_from_nodes(nmat)

    integer :: n,m

    ! tally up the vof
    vof_from_nodes = 0.0_r8
    do n = 1,8
      m = matl_geometry%index_at(this%node(:,n))
      vof_from_nodes(m) = vof_from_nodes(m) + 0.125_r8
    end do

  end function vof_from_nodes

  function cell_center(hex)
    ! returns the center of the given hex
    class(dnc_hex), intent(in) :: hex
    real(r8) :: cell_center(3)
    cell_center = sum(hex%node(:,:), dim=2) / 8.0_r8
  end function cell_center

  function face_centers(this)
    ! returns an array of face centers
    class(dnc_hex), intent(in) :: this
    real(r8) :: face_centers(3,6)
    face_centers(:,1) = sum(this%node(:,(/1,2,5,6/)), dim=2) / 4.0_r8
    face_centers(:,2) = sum(this%node(:,(/2,3,6,7/)), dim=2) / 4.0_r8
    face_centers(:,3) = sum(this%node(:,(/3,4,7,8/)), dim=2) / 4.0_r8
    face_centers(:,4) = sum(this%node(:,(/4,1,5,8/)), dim=2) / 4.0_r8
    face_centers(:,5) = sum(this%node(:,(/1,2,3,4/)), dim=2) / 4.0_r8
    face_centers(:,6) = sum(this%node(:,(/5,6,7,8/)), dim=2) / 4.0_r8
  end function face_centers

  function edge_centers(this)
    ! returns an array of edge centers
    class(dnc_hex), intent(in) :: this
    real(r8) :: edge_centers(3,12)
    edge_centers(:, 1) = sum(this%node(:,(/1,2/)), dim=2) / 2.0_r8
    edge_centers(:, 2) = sum(this%node(:,(/2,3/)), dim=2) / 2.0_r8
    edge_centers(:, 3) = sum(this%node(:,(/3,4/)), dim=2) / 2.0_r8
    edge_centers(:, 4) = sum(this%node(:,(/4,1/)), dim=2) / 2.0_r8
    edge_centers(:, 5) = sum(this%node(:,(/1,5/)), dim=2) / 2.0_r8
    edge_centers(:, 6) = sum(this%node(:,(/2,6/)), dim=2) / 2.0_r8
    edge_centers(:, 7) = sum(this%node(:,(/3,7/)), dim=2) / 2.0_r8
    edge_centers(:, 8) = sum(this%node(:,(/4,8/)), dim=2) / 2.0_r8
    edge_centers(:, 9) = sum(this%node(:,(/5,6/)), dim=2) / 2.0_r8
    edge_centers(:,10) = sum(this%node(:,(/6,7/)), dim=2) / 2.0_r8
    edge_centers(:,11) = sum(this%node(:,(/7,8/)), dim=2) / 2.0_r8
    edge_centers(:,12) = sum(this%node(:,(/8,5/)), dim=2) / 2.0_r8
  end function edge_centers

  function divide(this) result(subhex)
    ! this function takes a hex and returns an array of 8 hexes obtained from dividing the input
    class(dnc_hex), intent(in) :: this
    type(dnc_hex)              :: subhex(8)

    real(r8) :: cellc(3), facec(3,6), edgec(3,12)

    ! find center points
    cellc = this%cell_center()
    facec = this%face_centers()
    edgec = this%edge_centers()

    ! set subhex node positions
    subhex(1)%node(:,1) = this%node(:,1)
    subhex(1)%node(:,2) = edgec(:,1)
    subhex(1)%node(:,3) = facec(:,5)
    subhex(1)%node(:,4) = edgec(:,4)
    subhex(1)%node(:,5) = edgec(:,5)
    subhex(1)%node(:,6) = facec(:,1)
    subhex(1)%node(:,7) = cellc(:)
    subhex(1)%node(:,8) = facec(:,4)

    subhex(2)%node(:,1) = edgec(:,1)
    subhex(2)%node(:,2) = this%node(:,2)
    subhex(2)%node(:,3) = edgec(:,2)
    subhex(2)%node(:,4) = facec(:,5)
    subhex(2)%node(:,5) = facec(:,1)
    subhex(2)%node(:,6) = edgec(:,6)
    subhex(2)%node(:,7) = facec(:,2)
    subhex(2)%node(:,8) = cellc(:)

    subhex(3)%node(:,1) = facec(:,5)
    subhex(3)%node(:,2) = edgec(:,2)
    subhex(3)%node(:,3) = this%node(:,3)
    subhex(3)%node(:,4) = edgec(:,3)
    subhex(3)%node(:,5) = cellc(:)
    subhex(3)%node(:,6) = facec(:,2)
    subhex(3)%node(:,7) = edgec(:,7)
    subhex(3)%node(:,8) = facec(:,3)

    subhex(4)%node(:,1) = edgec(:,4)
    subhex(4)%node(:,2) = facec(:,5)
    subhex(4)%node(:,3) = edgec(:,3)
    subhex(4)%node(:,4) = this%node(:,4)
    subhex(4)%node(:,5) = facec(:,4)
    subhex(4)%node(:,6) = cellc(:)
    subhex(4)%node(:,7) = facec(:,3)
    subhex(4)%node(:,8) = edgec(:,8)

    subhex(5)%node(:,1) = edgec(:,1)
    subhex(5)%node(:,2) = facec(:,1)
    subhex(5)%node(:,3) = cellc(:)
    subhex(5)%node(:,4) = facec(:,4)
    subhex(5)%node(:,5) = this%node(:,5)
    subhex(5)%node(:,6) = edgec(:,9)
    subhex(5)%node(:,7) = facec(:,6)
    subhex(5)%node(:,8) = edgec(:,12)

    subhex(6)%node(:,1) = facec(:,1)
    subhex(6)%node(:,2) = edgec(:,6)
    subhex(6)%node(:,3) = facec(:,2)
    subhex(6)%node(:,4) = cellc(:)
    subhex(6)%node(:,5) = edgec(:,9)
    subhex(6)%node(:,6) = this%node(:,6)
    subhex(6)%node(:,7) = edgec(:,10)
    subhex(6)%node(:,8) = facec(:,6)

    subhex(7)%node(:,1) = cellc(:)
    subhex(7)%node(:,2) = facec(:,2)
    subhex(7)%node(:,3) = edgec(:,7)
    subhex(7)%node(:,4) = facec(:,3)
    subhex(7)%node(:,5) = facec(:,6)
    subhex(7)%node(:,6) = edgec(:,10)
    subhex(7)%node(:,7) = this%node(:,7)
    subhex(7)%node(:,8) = edgec(:,11)

    subhex(8)%node(:,1) = facec(:,4)
    subhex(8)%node(:,2) = cellc(:)
    subhex(8)%node(:,3) = facec(:,3)
    subhex(8)%node(:,4) = edgec(:,8)
    subhex(8)%node(:,5) = edgec(:,12)
    subhex(8)%node(:,6) = facec(:,6)
    subhex(8)%node(:,7) = edgec(:,11)
    subhex(8)%node(:,8) = this%node(:,8)
  end function divide

end module vof_init
