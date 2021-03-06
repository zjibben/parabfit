!!
!! SCALAR_MESH_FUNC_CLASS
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! September 2014
!!

module scalar_mesh_func_class

  use kinds, only: r8
  implicit none
  private

  type, abstract, public :: scalar_mesh_func
    real(r8), allocatable, public :: value(:)
  contains
    procedure(compute), deferred :: compute
  end type scalar_mesh_func

  abstract interface
    subroutine compute (this, t)
      import r8, scalar_mesh_func
      class(scalar_mesh_func), intent(inout) :: this
      real(r8), intent(in) :: t
    end subroutine
  end interface

end module scalar_mesh_func_class
