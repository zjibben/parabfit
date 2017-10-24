module mixed_cell_subset_constructor

  use kinds, only: r8
  use unstr_mesh_type
  use mesh_subset_type
  implicit none
  private

  public :: mixed_cell_subset

  type, extends(subset_check) :: mixed_cell_check
    private
    real(r8) :: cutvof
    real(r8), pointer :: vof(:,:) ! reference only -- do not own
  contains
    procedure :: init => init_mixed_cell_check
    procedure :: inside_subset => is_mixed_cell
  end type mixed_cell_check

contains

  type(mesh_subset) function mixed_cell_subset(vof, mesh)

    use consts, only: cutvof

    real(r8), intent(in), target :: vof(:,:)
    type(unstr_mesh) :: mesh

    type(mixed_cell_check) :: check

    call check%init(vof, cutvof)
    call mixed_cell_subset%init(check, mesh)

  end function mixed_cell_subset

  subroutine init_mixed_cell_check(this, vof, cutvof)

    class(mixed_cell_check), intent(out) :: this
    real(r8), intent(in), target :: vof(:,:)
    real(r8), intent(in) :: cutvof

    this%vof => vof
    this%cutvof = cutvof

  end subroutine init_mixed_cell_check

  logical function is_mixed_cell(this, i)
    class(mixed_cell_check), intent(in) :: this
    integer, intent(in) :: i
    is_mixed_cell = any(this%vof(:,i) > this%cutvof .and. this%vof(:,i) < 1 - this%cutvof)
  end function is_mixed_cell

end module mixed_cell_subset_constructor
