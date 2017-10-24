!!
!! MESH_SUBSET_TYPE
!!
!! This module defines a type for managing subsets of the mesh.
!! It also provides a class for creating logical tests for defining subsets.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2017
!!

module mesh_subset_type

  implicit none
  private

  type, public :: mesh_subset
    integer :: ncell
    integer, allocatable :: cell_id(:)
  contains
    procedure :: init
  end type mesh_subset

  type, abstract, public :: subset_check
  contains
    procedure(is_inside_subset), deferred :: inside_subset
  end type subset_check

  abstract interface
    logical function is_inside_subset(this,i)
      import subset_check
      class(subset_check), intent(in) :: this
      integer, intent(in) :: i
    end function is_inside_subset
  end interface

contains

  subroutine init(this, check, mesh)

    use unstr_mesh_type

    class(mesh_subset), intent(out) :: this
    class(subset_check), intent(in) :: check
    type(unstr_mesh), intent(in) :: mesh

    integer :: i, j

    ! count the number of cells in this subset
    this%ncell = 0
    do i = 1,mesh%ncell
      if (check%inside_subset(i)) this%ncell = this%ncell + 1
    end do

    ! list the cells
    allocate(this%cell_id(this%ncell))
    j = 1
    do i = 1,mesh%ncell
      if (check%inside_subset(i)) then
        this%cell_id(j) = i
        j = j + 1
      end if
    end do

  end subroutine init

end module mesh_subset_type
