#include "f90_assert.fpp"

module bfgs_min_class_test

  use kinds, only: r8
  use bfgs_min_class
  implicit none
  private

  public :: bfgs_min_class_test_suite

  type, extends(bfgs_min) :: bfgs_test_func
  contains
    procedure :: f => bfgs_test_f
  end type bfgs_test_func

contains

  subroutine bfgs_min_class_test_suite()

    integer :: ierr
    real(r8) :: x(2)
    type(bfgs_test_func) :: f_test

    x = [0.5_r8, -0.5_r8]
    call f_test%find_minimum(x, ierr)

    print *, x
    print *, ierr, f_test%numitr
    print *, norm2(x - [0.1_r8,-0.4_r8])

  end subroutine bfgs_min_class_test_suite

  real(r8) function bfgs_test_f(this, x)

    class(bfgs_test_func), intent(inout) :: this
    real(r8), intent(in) :: x(:)

    ASSERT(size(x) == 2)

    bfgs_test_f = 0.2_r8 - cos(2*(x(1)-0.1_r8)**2 + (x(2)+0.4_r8)**2)

  end function bfgs_test_f

end module bfgs_min_class_test
