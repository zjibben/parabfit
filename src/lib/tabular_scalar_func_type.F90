!!
!! TABULAR_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class SCALAR_FUNC.  This
!! implementation defines a tabular function given by user-specified data
!! points with intervening linear interpolation.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!

#include "f90_assert.fpp"

module tabular_scalar_func_type

  use kinds, only: r8
  use scalar_func_class
  implicit none
  private

  type, extends(scalar_func), public :: tabular_scalar_func
    private
    real(r8), allocatable :: x(:), y(:)
  contains
    procedure :: eval
  end type tabular_scalar_func

  !! Defined constructor
  interface tabular_scalar_func
    procedure tabular_scalar_func_value
  end interface

contains

  !! Constructor for TABULAR_SCALAR_FUNC objects
  function tabular_scalar_func_value (x, y) result (f)
    real(r8), intent(in) :: x(:), y(:)
    type(tabular_scalar_func) :: f
    ASSERT(size(x) > 1)
    ASSERT(size(y) == size(x))
    f%x = x
    f%y = y
  end function tabular_scalar_func_value

  function eval (this, x) result (fx)
    class(tabular_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)  ! only x(1) is used
    real(r8) :: fx
    integer :: n, j, j1, j2
    n = size(this%x)
    if (x(1) <= this%x(1)) then
      fx = this%y(1)
    else if (x(1) >= this%x(n)) then
      fx = this%y(n)
    else
      !! Binary search to find the interval x(j1) < x <= x(j2), j2 = j1+1.
      j1 = 1; j2 = n
      do while (j2 - j1 > 1)
        j = (j1 + j2) / 2
        if (x(1) > this%x(j)) then
          j1 = j
        else
          j2 = j
        end if
      end do
      !! Linearly interpolate over the interval
      fx = ((this%x(j2)-x(1))*this%y(j1) + (x(1)-this%x(j1))*this%y(j2))/(this%x(j2)-this%x(j1))
    end if
  end function eval

end module tabular_scalar_func_type
