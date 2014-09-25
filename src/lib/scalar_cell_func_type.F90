module scalar_cell_func_type

  use kinds, only: r8
  use scalar_mesh_func_class
  implicit none
  private

  type, extends(scalar_mesh_func), public :: scalar_cell_func
    type(unstr_mesh), pointer :: mesh => null() ! reference only -- do not own
  contains
    procedure :: init
    procedure :: add
    procedure :: add_complete
    procedure :: compute
  end type scalar_cell_func

contains

  subroutine init (this, mesh)
    class(scalar_cell_func), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh

    this%mesh => mesh
  end subroutine init

  subroutine add (this, f, setids, stat, errmsg)
    class(scalar_cell_func), intent(inout) :: this
    class(scalar_func), allocatable, intent(inout) :: f
    integer, intent(in) :: setids(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg
  end subroutine add

end module scalar_cell_func_type
