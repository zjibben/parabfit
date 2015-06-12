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
  use vof_tools, only: cell_materials,distinct_matls,distinct_entries,index_of
  use material_geometry_type
  use logging_services
  implicit none
  private

  ! hex type to make divide and conquer algorithm simpler
  type, public :: hex_cell
     real(r8) :: node(3,8)
    integer  :: depth
   contains
     procedure         :: divide
     procedure         :: volume
     procedure         :: cell_center
     procedure         :: face_centers
     procedure         :: edge_centers
     procedure         :: contains_interface
     procedure         :: vof
     procedure,private :: vof_from_nodes
  end type hex_cell

  public :: vof_initialize

  integer, parameter :: cell_vof_recursion_limit = 4 ! should be made user-specified parameter
contains

  subroutine vof_initialize (mesh, plist, cell_matls)
    use unstr_mesh_type
    
    ! this subroutine initializes the vof values in all cells
    ! from a user input function. It calls a divide and conquer algorithm
    ! to calculate the volume fractions given an interface
    type(unstr_mesh), pointer, intent(in)             :: mesh
    type(parameter_list), intent(in)                  :: plist
    type(cell_materials), dimension(:), intent(inout) :: cell_matls

    ! local variables
    integer :: i
    type(hex_cell) :: hex
    type(material_geometry) :: matl_init_geometry
    
    ! first, determine the initial state provided by the user, and assign
    ! the function which will determine what points are inside the materials
    call matl_init_geometry%init(plist)
    
    ! next, loop though every cell and check if all verteces lie within a single material
    ! if so, set the Vof for that material to 1.0.
    ! if not, divide the hexahedral cell into 8 subdomains defined by the centroid and centers of faces
    ! repeat the process recursively for each subdomain up to a given threshold
    do i = 1,mesh%ncell
       ! make a hex type out of the cell before calculating the vof
       ! this will make it easier to handle
       hex%node = mesh%x(:,mesh%cnode(:,i))
       hex%depth = 0

       ! calculate the vof
       cell_matls(i) = hex%vof(matl_init_geometry)
    end do
    
  end subroutine vof_initialize

  recursive function vof(this, matl_geometry) result(hex_vof)
    ! calculates the volume fractions of materials in a cell
    class(hex_cell), intent(in) :: this
    type(material_geometry), intent(in) :: matl_geometry
    type(cell_materials)              :: hex_vof

    integer :: i,m,hm
    type(hex_cell)      , dimension(8) :: subhex
    type(cell_materials), dimension(8) :: subhex_vof

    ! if the cell contains an interface (and therefore has at least two materials
    ! and we haven't yet hit our recursion limit, divide the hex and repeat
    if (this%contains_interface(matl_geometry) .and. this%depth < cell_vof_recursion_limit) then
       !if (this%depth <= cell_vof_recursion_limit) then
       ! divide into 8 smaller hexes
       subhex = this%divide()
       
       ! calculate the vof in the subhexes
       do i = 1,8
          subhex_vof(i) = subhex(i)%vof(matl_geometry)
       end do
       
       ! get distinct materials and number of distinct materials from subhexes
       call distinct_matls(hex_vof%nmat, hex_vof%matl, subhex_vof)
       
       ! build up current hex's vof from information given by subhexes
       hex_vof%matl(:)%vof = 0.0_r8
       do i = 1,8
          do m = 1,subhex_vof(i)%nmat
             hm = index_of(subhex_vof(i)%matl(m)%id, hex_vof%matl(:)%id)
             hex_vof%matl(hm)%vof = hex_vof%matl(hm)%vof + 0.125_r8*subhex_vof(i)%matl(m)%vof
          end do
       end do
    else
       ! if we are at (or somehow past) the recursion limit or the cell does not contain an interface,
       ! calculate the vof in this hex based on the materials at its nodes
       hex_vof = this%vof_from_nodes(matl_geometry)
    end if
    
  end function vof
  
  logical function contains_interface(this, matl_geometry)
    ! This function takes a hex, described by its vertex positions,
    ! and determines if it contains an interface for a given material.
    ! This is done by checking if every vertex lies within the material, or not.
    class(hex_cell), intent(in) :: this
    type(material_geometry), intent(in) :: matl_geometry

    integer :: i,reference_id
    real(r8) :: facec(3,6),edgec(3,12)
    
    contains_interface = .false.
    reference_id = matl_geometry%material_at(this%node(:,1))

    ! if any node doesn't have the same material as the first node, there's an interface in this cell
    do i = 2,8
       contains_interface = contains_interface .or. reference_id /= matl_geometry%material_at(this%node(:,i))
    end do
    
  end function contains_interface
  
  type(cell_materials) function vof_from_nodes(this, matl_geometry)
    ! the material at each node contributes 1/8th of the vof in the given hex
    class(hex_cell), intent(in) :: this
    type(material_geometry), intent(in) :: matl_geometry
    
    integer :: n,m
    integer, dimension(8) :: matl_at_node
    
    ! get the material ids at each node
    do n = 1,8
       matl_at_node(n) = matl_geometry%material_at(this%node(:,n))
    end do
    
    ! count and store the distinct materials on the nodes
    call distinct_entries(vof_from_nodes%matl, matl_at_node)
    vof_from_nodes%nmat = size(vof_from_nodes%matl)
    
    ! sort?
    
    ! tally up the vof
    vof_from_nodes%matl(:)%vof = 0.0_r8
    do n = 1,8
       ! get the material at each node consider it to be taking up 1/8th of the volume
       m = index_of(matl_at_node(n), vof_from_nodes%matl(:)%id)
       vof_from_nodes%matl(m)%vof = vof_from_nodes%matl(m)%vof + 0.125_r8
    end do
    
  end function vof_from_nodes
  
  real(r8) function volume (this)
    ! calculates the volume of a hex
    use cell_geometry, only: eval_hex_volumes
    class(hex_cell), intent(in) :: this
    real(r8) :: cvol_tmp(8)
    call eval_hex_volumes(this%node, volume, cvol_tmp)
  end function volume
  
  function cell_center(hex)
    ! returns the center of the given hex
    class(hex_cell), intent(in) :: hex
    real(r8), dimension(3) :: cell_center
    cell_center = sum(hex%node(:,:), dim=2) / 8.0_r8
  end function cell_center

  function face_centers(this)
    ! returns an array of face centers
    class(hex_cell), intent(in) :: this
    real(r8), dimension(3,6) :: face_centers
    face_centers(:,1) = sum(this%node(:,(/1,2,5,6/)), dim=2) / 4.0_r8
    face_centers(:,2) = sum(this%node(:,(/2,3,6,7/)), dim=2) / 4.0_r8
    face_centers(:,3) = sum(this%node(:,(/3,4,7,8/)), dim=2) / 4.0_r8
    face_centers(:,4) = sum(this%node(:,(/4,1,5,8/)), dim=2) / 4.0_r8
    face_centers(:,5) = sum(this%node(:,(/1,2,3,4/)), dim=2) / 4.0_r8
    face_centers(:,6) = sum(this%node(:,(/5,6,7,8/)), dim=2) / 4.0_r8
  end function face_centers

  function edge_centers(this)
    ! returns an array of edge centers
    class(hex_cell), intent(in) :: this
    real(r8), dimension(3,12) :: edge_centers
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
    class(hex_cell), intent(in)  :: this
    type(hex_cell), dimension(8) :: subhex
    
    real(r8) :: cellc(3), facec(3,6), edgec(3,12)
    integer :: i
    
    ! find center points
    cellc = this%cell_center()
    facec = this%face_centers()
    edgec = this%edge_centers()
    
    ! set subhexes at 1 lower depth
    subhex(:)%depth = this%depth+1

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
