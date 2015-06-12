!!
!!
!!
!!
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

module vof_tools
  use kinds, only: r8
  use logging_services
  implicit none
  private
  
  ! the material type contains the material id and vof
  type :: material
     integer  :: id
     real(r8) :: vof
  end type material

  ! the cell_materials type contains the number of materials present in a given cell
  ! plus the material id and vof for each fluid in that cell
  type, public :: cell_materials
     integer                                   :: nmat
     type(material), dimension(:), allocatable :: matl
   contains
     final :: delete_cell_materials
  end type cell_materials

  public :: matl_get_vof
  public :: distinct_matls
  public :: distinct_entries
  public :: index_of
  public :: volume_of_matl

  interface distinct_matls
     procedure distinct_matls_material
     procedure distinct_matls_integer
  end interface distinct_matls

  interface distinct_entries
     procedure distinct_entries_material
     procedure distinct_entries_integer
  end interface distinct_entries
  
contains

  subroutine delete_cell_materials (this)
    type(cell_materials), intent(inout) :: this
    if (allocated(this%matl)) deallocate(this%matl)
  end subroutine delete_cell_materials
  
  subroutine matl_get_vof(vof, cell_matls, nmat, matl_id)
    ! this subroutine takes the ragged cell_matls array and returns
    ! an equivalent regular array
    real(r8), dimension(:,:), intent(inout) :: vof
    type(cell_materials), dimension(:), intent(in) :: cell_matls
    integer, dimension(:), intent(in) :: matl_id
    integer, intent(in) :: nmat

    integer :: i,m

    vof = 0.0_r8
    do i = 1,size(cell_matls)
       do m = 1,cell_matls(i)%nmat
          vof(i,index_of(cell_matls(i)%matl(m)%id, matl_id)) = cell_matls(i)%matl(m)%vof
       end do
    end do
    
  end subroutine matl_get_vof

  
  subroutine distinct_matls_integer(nmat, matls, cell_matls)
    integer,                            intent(out) :: nmat
    integer, dimension(:), allocatable, intent(out) :: matls
    type(cell_materials), dimension(:), intent(in)  :: cell_matls

    integer :: i,n,m, nmat_all
    integer, dimension(:), allocatable :: matl_ids

    ! first make a flat array of materials in each cell
    nmat_all = sum(cell_matls(:)%nmat)
    allocate(matl_ids(nmat_all))
    
    i = 1
    do n = 1,size(cell_matls)
       do m = 1,cell_matls(n)%nmat
          matl_ids(i) = cell_matls(n)%matl(m)%id
          i = i+1
       end do
    end do

    ! count distinct materials
    call distinct_entries(matls, matl_ids)
    nmat = size(matls) !ndistinct(matl_ids)

    ! clean up
    deallocate(matl_ids)
    
  end subroutine distinct_matls_integer

  subroutine distinct_matls_material(nmat, matls, cell_matls)
    integer,                            intent(out) :: nmat
    type(material), dimension(:), allocatable, intent(out) :: matls
    type(cell_materials), dimension(:), intent(in)  :: cell_matls

    integer :: i,n,m, nmat_all
    integer, dimension(:), allocatable :: matl_ids

    ! first make a flat array of materials in each cell
    nmat_all = sum(cell_matls(:)%nmat)
    allocate(matl_ids(nmat_all))
    
    i = 1
    do n = 1,size(cell_matls)
       do m = 1,cell_matls(n)%nmat
          matl_ids(i) = cell_matls(n)%matl(m)%id
          i = i+1
       end do
    end do

    ! count distinct materials
    call distinct_entries(matls, matl_ids)
    nmat = size(matls) !ndistinct(matl_ids)
    
    ! clean up
    deallocate(matl_ids)
    
  end subroutine distinct_matls_material

  integer function ndistinct(array)
    ! return the number of distinct values in an integer array
    integer, dimension(:), intent(in) :: array

    integer :: i
    
    ndistinct = 1
    do i = 2,size(array)
       if (.not.any(array(1:i-1) == array(i))) then
          ndistinct = ndistinct+1
       end if
    end do
    
  end function ndistinct
  
  subroutine distinct_entries_integer(distinct, array)
    ! return an array of the distinct entries in an integer array
    integer, dimension(:), allocatable, intent(out) :: distinct
    integer, dimension(:), intent(in) :: array

    integer :: i
    
    allocate(distinct(1))
    
    distinct(1) = array(1)
    do i = 2,size(array)
       if (.not.any(distinct == array(i))) call append(distinct, array(i))
    end do
    
  end subroutine distinct_entries_integer

  subroutine append(array, new_entry)
    integer, dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: new_entry

    if (allocated(array)) then
       call reallocate(array, size(array)+1)
    else
       allocate(array(1))
    end if
    array(size(array)) = new_entry
  end subroutine append
  
  subroutine reallocate(array, new_size)
    integer, dimension(:), allocatable, intent(inout) :: array
    integer, intent(in) :: new_size

    integer, dimension(:), allocatable :: tmp
    
    allocate(tmp(new_size))
    tmp(1:size(array)) = array
    !deallocate(array)
    call move_alloc(tmp, array)
  end subroutine reallocate

  
  subroutine distinct_entries_material(distinct, array)
    ! return an array of the distinct entries in an integer array
    type(material), dimension(:), allocatable, intent(out) :: distinct
    integer       , dimension(:), intent(in) :: array

    integer :: i,j
    
    allocate(distinct(ndistinct(array)))

    distinct(1)%id = array(1)
    j = 2
    do i = 2,size(array)
       if (.not.any(array(1:i-1) == array(i))) then
          distinct(j)%id = array(i)
          j = j+1
       end if
    end do
    
  end subroutine distinct_entries_material
  
  integer function index_of(n, array)
    ! returns the index of the first entry of value n in array
    ! returns -1 if not found
    integer, intent(in) :: n
    integer, dimension(:), intent(in) :: array

    integer :: i

    do i = 1,size(array)
       if (array(i)==n) then
          index_of = i
          return
       end if
    end do

    ! value not found
    index_of = -1
    !call LS_fatal ('value not found in array')    
    
  end function index_of

  real(r8) function volume_of_matl(id, cell_matls, cell_vol)
    ! returns the total volume of a material through the domain
    integer, intent(in) :: id
    type(cell_materials), dimension(:), intent(in) :: cell_matls
    real(r8), dimension(:), intent(in)             :: cell_vol

    integer :: i,m

    volume_of_matl = 0.0_r8
    do i = 1,size(cell_matls)
       m = index_of(id, cell_matls(i)%matl(:)%id)
       if (m>0) volume_of_matl = volume_of_matl + cell_matls(i)%matl(m)%vof * cell_vol(i)
    end do
    
  end function volume_of_matl
  
end module vof_tools
