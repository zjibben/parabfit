!!
!! PIECEWISE_SCALAR_FUNC_TYPE
!!
!! A concrete implementation of the abstract base class SCALAR_FUNC.  This
!! implementation defines a piecewise function 
!!
!! Zechariah Jibben <zjibben@lanl.gov>
!! April 2016
!!

#include "f90_assert.fpp"

module piecewise_scalar_func_type

  use kinds, only: r8
  use scalar_func_class
  use scalar_func_containers
  use region_class
  implicit none
  private
  
  type, extends(scalar_func), public :: piecewise_scalar_func
    private
    type(scalar_func_box), allocatable :: f(:) ! subfunctions
    type(region_box),      allocatable :: r(:) ! regions
  contains
    procedure :: eval
  end type piecewise_scalar_func
  
  !! Defined constructor
  interface piecewise_scalar_func
    procedure piecewise_scalar_func_value
  end interface piecewise_scalar_func

contains

  !! Constructor for PIECEWISE_SCALAR_FUNC objects.
  function piecewise_scalar_func_value (subfunc, region) result (f)

    type(scalar_func_box), intent(in) :: subfunc(:)
    type(region_box),      intent(in) :: region(:)
    type(piecewise_scalar_func) :: f

    INSIST(size(subfunc)==size(region))

    f%f = subfunc
    f%r = region

  end function piecewise_scalar_func_value

  function eval (this, x) result (fx)

    use logging_services

    class(piecewise_scalar_func), intent(in) :: this
    real(r8), intent(in) :: x(:)
    real(r8) :: fx

    integer :: r

    do r = 1,size(this%r)
      if (this%r(r)%r%location_is_inside (x)) then
        fx = this%f(r)%f%eval (x)
        return
      end if
    end do

    call LS_fatal ('given position not in a given region')

  end function eval


end module piecewise_scalar_func_type
