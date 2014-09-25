!!
!! UNSTR_MESH_FUNC
!!
!! Provides a parameter list-driven procedure for defining a cell-based array.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Aug 2013.  Adapted Sep 2014.
!!
!! PROGRAMMING INTERFACE
!!
!!  The subroutine COMPUTE_MESH_FUNC defines the values of a cell-based array
!!  according to the specification given by a parameter list.  This is a quick
!!  and dirty implementation that needs refinement; e.g., inappropriate use of
!!  assertions for input error checking (FIXME).
!!
!!  The expected parameter list consists of one or more sublists.  Each sublist
!!  defines the value of the array over a region of the mesh.  The sublist
!!  regions must be disjoint and cover the entire mesh.  Each sublist consists
!!  of the parameter "cell-blocks" and either the parameter "constant" or a
!!  sublist "function".  The "cell-blocks" parameter is an array of cell block
!!  IDs that define the region.  To define a constant value for all cells in
!!  that region assign that value to the parameter "constant".  Alternatively,
!!  if the "function" sublist is defined, its value is used to define a scalar
!!  function of (x,y,z) whose value at a cell center is assigned to the array
!!  element corresponding to the cell.  For the range of functions that can be
!!  defined see scalar_func_factories.F90.
!!

#include "f90_assert.fpp"

module unstr_mesh_func

  use kinds, only: r8
  use unstr_mesh_type
  use scalar_func_class
  use scalar_func_factories
  use parameter_list_type
  implicit none
  private

  public :: compute_mesh_func

contains

  subroutine compute_mesh_func (mesh, params, f)

    type(unstr_mesh), intent(in) :: mesh
    type(parameter_list), intent(inout) :: params
    real(r8), intent(out) :: f(:)

    logical, allocatable :: mask(:)
    type(parameter_list_iterator) :: regit
    type(parameter_list), pointer :: sublist
    class(scalar_func), allocatable :: func
    real(r8) :: value, xc(size(mesh%x,dim=1))
    integer, allocatable :: blkids(:)
    integer :: n, j

    ASSERT(size(f) == mesh%ncell)

    allocate(mask(mesh%ncell))
    mask = .false.
    regit = parameter_list_iterator(params)!, sublists_only=.true.)
    do while (.not.regit%at_end())
      sublist => regit%sublist()
      call sublist%get ('cell-blocks', blkids)
      if (sublist%is_sublist('function')) then
        sublist => sublist%sublist('function')
        call alloc_scalar_func (func, sublist)
        do n = 1, size(blkids)
          INSIST(.not.any(mesh%cblock == blkids(n) .and. mask))
          do j = 1, mesh%ncell
            if (mesh%cblock(j) /= blkids(n)) cycle
            xc = sum(mesh%x(:,mesh%cnode(:,j)),dim=2) / size(mesh%cnode,dim=1)
            f(j) = func%eval(xc)
          end do
          mask = mask .or. (mesh%cblock == blkids(n))
        end do
      else
        call sublist%get ('constant', value)
        do n = 1, size(blkids)
          INSIST(.not.any(mesh%cblock == blkids(n) .and. mask))
          where (mesh%cblock == blkids(n)) f = value
          mask = mask .or. (mesh%cblock == blkids(n))
        end do
      end if
      call regit%next()
    end do
    INSIST(all(mask))

  end subroutine compute_mesh_func

end module unstr_mesh_func
