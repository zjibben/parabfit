!!
!! BFGS_MIN_CLASS
!!
!! This module provides an interface to the BFGS minimization algorithm,
!! as described by Kelley, Iterative Methods for Optimization (1999).
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! August 2017
!!

module bfgs_min_class

  use kinds, only: r8
  implicit none
  private

  type, abstract, public :: bfgs_min
    real(r8) :: tol = 1e-6_r8
    integer :: maxitr = 20
    integer :: line_search_max = 20
    integer :: numitr = 0
  contains
    procedure, non_overridable :: find_minimum
    procedure(func), deferred :: f
    procedure :: gradf ! in general this should be overridable. here estimate it with differencing.
    procedure, private :: line_search
  end type bfgs_min

  abstract interface
    function func (this, x) result(fx)
      import r8, bfgs_min
      class(bfgs_min), intent(in) :: this
      real(r8), intent(in) :: x(:)
      real(r8) :: fx
    end function func
  end interface

contains

  subroutine find_minimum(this, x, status)

    use array_utils, only: outer_product

    class(bfgs_min), intent(inout) :: this
    real(r8), intent(inout) :: x(:)
    integer, intent(out) :: status

    integer :: i
    real(r8) :: l, f, gradf(size(x)), d(size(x)), s(size(x)), y(size(x)), xold(size(x)), &
        hess_inv(size(x),size(x)), identity(size(x),size(x)), &
        gradfnew(size(x)), &
        xhist(size(x), this%maxitr+1) ! DEBUGGING

    identity = 0
    do i = 1,size(x)
      identity(i,i) = 1
    end do
    status = 0
    f = this%f(x)
    ! WARN: initial guess for a good difference step size here is assumed
    gradf = this%gradf(x, 1e-7_r8)
    hess_inv = identity

    xhist(:,1) = x
    ! print '(a,2es13.3)', 'x: ', x
    ! print '(a,2es13.3)', 'f: ', f
    ! print '(a,3es13.3)', 'g: ', gradf, norm2(gradf)

    if (norm2(gradf) < this%tol) return

    do i = 1,this%maxitr
      ! calculate the search direction
      d = - matmul(hess_inv, gradf)

      ! get new position
      xold = x
      call this%line_search(l, s, f, x, d, gradf, status)
      xhist(:,i+1) = x
      if (status==1) exit

      ! update the hessian inverse and gradient
      gradfnew = this%gradf(x, 1e-7_r8) ! 1e-7_r8
      y = gradfnew - gradf
      gradf = gradfnew
      !gradf = gradf + y
      if (dot_product(y,s) > 0) then
        hess_inv = matmul(hess_inv, identity - outer_product(y,s) / dot_product(y,s))
        hess_inv = matmul(identity - outer_product(s,y) / dot_product(y,s), hess_inv)
        hess_inv = hess_inv + outer_product(s,s) / dot_product(y,s)
      end if

      ! print *, i
      ! print '(a,2es13.3)', 'd: ', d
      ! print '(a,3es13.3)', 's: ', s, norm2(s)
      ! print '(a,2es13.3)', 'x: ', x
      ! print '(a,2es13.3)', 'f: ', f
      ! print '(a,3es13.3)', 'g: ', gradf, norm2(gradf)
      ! print '(a,2es13.3)', 'ys: ', dot_product(y,s)
      ! print *

      if (norm2(gradf) < this%tol .or. norm2(x-xold) < 1e-10_r8) exit
    end do

    this%numitr = i
    if (this%numitr > this%maxitr) then
      !print *, "too many iterations in bfgs", i
      status = 1
    end if

    ! do i = 1,this%numitr+1
    !   print *, '[',xhist(1,i), ', ', xhist(2,i),'],'
    ! end do

    ! do i = 2,this%numitr+1
    !   print *, norm2(xhist(:,i) - xhist(:,i-1))
    ! end do

  end subroutine find_minimum

  subroutine line_search(this, l, s, f, x, d, gradf, status)

    class(bfgs_min), intent(in) :: this
    real(r8), intent(out) :: l, s(:)
    real(r8), intent(inout) :: f, x(:)
    real(r8), intent(in) :: d(:), gradf(:)
    integer, intent(out) :: status

    real(r8), parameter :: blow = 0.1_r8, bhigh = 0.5_r8
    integer :: i
    real(r8) :: goalval, fnext, xnext(size(x)), c(4), fprev, lprev, b(2), Ainv(2,2)

    status = 0
    goalval = 1e-4_r8 * dot_product(gradf,d)
    ! print *, 'line-search g:    ', goalval

    ! initial guess for step length
    l = min(1.0_r8, 100/(1+norm2(gradf)))
    s = l*d
    xnext = x + s
    fnext = this%f(xnext)
    if (fnext <= f + l*goalval) then
      !print '(a,es13.3)', 'df: ', fnext - f
      x = xnext
      f = fnext
      return
    end if
    ! print *, 'line-search l:    ', l

    ! quadratic guess for step length
    c(1) = f ! xi(0)
    c(2) = dot_product(gradf,d) ! xi'(0)
    c(3) = (fnext - c(1) - c(2)) / l**2

    lprev = l
    fprev = fnext
    l = min(max(- c(2) / (2*c(3)), l*blow), l*bhigh)
    s = l*d
    xnext = x + s
    fnext = this%f(xnext)
    if (fnext <= f + l*goalval) then
      !print '(a,es13.3)', 'df: ', fnext - f
      x = xnext
      f = fnext
      return
    end if
    ! print *, 'line-search l:    ', l
    ! print *, 'line-search g:    ', goalval

    ! cubic guesses for step length
    ! print '(a,5es13.3)', 'line-search d:    ', d
    ! print '(a,5es13.3)', 'line-search g:    ', gradf
    do i = 1,this%line_search_max
      b = [fnext - c(1) - c(2)*l, fprev - c(1) - c(2)*lprev]
      Ainv(:,1) = [lprev**3, -lprev**2]
      Ainv(:,2) = [-l**3, l**2]
      Ainv = Ainv / (lprev**3*l**2 - l**3*lprev**2)
      c(3:4) = matmul(Ainv,b)

      ! local minimum of cubic
      lprev = l
      fprev = fnext
      l = min(max((-c(3) + sqrt(c(3)**2 - 3*c(4)*c(2))) / (3*c(4)), l*blow), l*bhigh)

      s = l*d
      xnext = x + s
      fnext = this%f(xnext)

      ! print '(a,i6)',      'line-search i:    ', i
      ! print '(a,5es13.3)', 'line-search c:    ', c
      ! print '(a,5es13.3)', 'line-search f:    ', fnext, f, fnext - f
      ! print '(a,5es13.3)', 'line-search x:    ', norm2(s)
      ! print '(a,6es13.3)', 'line-search diff: ', fnext - f, l*goalval, l, &
      !     (-c(3) + sqrt(c(3)**2 - 3*c(4)*c(2))) / (3*c(4)) / lprev, &
      !     l/lprev
      if (fnext <= f + l*goalval) exit
    end do

    if (i > this%line_search_max) then
      !print *, "too many backtracks in line search"
      status = 1
    end if
    ! print *
    ! print '(a,es13.3)', 'df: ', fnext - f
    x = xnext
    f = fnext

  end subroutine line_search

  ! estimate the gradient of f using standard differencing
  function gradf(this, x, dx)

    class(bfgs_min), intent(in) :: this
    real(r8), intent(in) :: x(:), dx
    real(r8) :: gradf(size(x))

    integer :: d
    real(r8) :: xt(size(x))

    xt = x
    do d = 1,size(x)
      ! positive side
      xt(d) = xt(d) + dx
      gradf(d) = this%f(xt)

      ! negative side
      xt(d) = xt(d) - 2*dx
      gradf(d) = gradf(d) - this%f(xt)

      gradf(d) = gradf(d) / (2*dx)

      ! restore xt
      xt(d) = xt(d) + dx
    end do

  end function gradf

end module bfgs_min_class
