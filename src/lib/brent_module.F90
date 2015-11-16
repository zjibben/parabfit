!!
!! BRENT_MODULE
!!
!! This module provides an interface to a Brent's algorithm.
!! Note this could probably be replaced with some canned package,
!! like the GNU Scientific Library or something else.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! Last revised 4 Nov 2012.
!!

module brent_module
  use kinds,  only: r8
  use consts, only: alittle,cutvof
  implicit none
  private

  public :: brent
  
  integer,  parameter :: iter_max = 100
  real(r8), parameter :: cgold = 1.0_r8 - 2.0_r8/(1.0_r8 + sqrt(5.0_r8))

  type, abstract, public :: brent_func
  contains
    procedure(eval), deferred :: eval
  end type brent_func

  abstract interface
    real(r8) function eval (this, x)
      import r8, brent_func
      class(brent_func), intent(in) :: this
      real(r8),    intent(in) :: x
    end function eval
  end interface

contains

  ! perform Brent's method (see Numerical Recipes)
  function brent (ax,bx,cx,f) result(x)
    real(r8),          intent(in) :: ax,bx,cx
    class(brent_func), intent(in) :: f

    real(r8) :: a,b
    real(r8) :: fu,fv,fw,fx, u,v,w,x,xm, p,q,r,tol1,tol2, etemp,d,e
    integer  :: iter

    a = min(ax,cx)
    b = max(ax,cx)
    x = bx; w = bx; v = bx;
    fx = f%eval (x); fw = fx; fv = fx
    e = 0.0_r8; d = 0.0_r8

    do iter = 1,iter_max
      xm = 0.5_r8 * (a+b)
      tol1 = cutvof*abs(x) + alittle
      tol2 = 2.0_r8 * tol1
      if (abs(x-xm) <= tol2-0.5_r8*(b-a)) return

      if (abs(e) > tol1) then
        r = (x-w)*(fx-fv)
        q = (x-v)*(fx-fw)
        p = (x-v)*q - (x-w)*r
        q = 2.0_r8*(q-r)
        if (q > 0.0_r8) p = -p
        q = abs(q)
        etemp = e
        e = d
        if (abs(p) >= abs(0.5_r8*q*etemp) .or. p <= q*(a-x) .or. p >= q*(b-x)) then
          e = merge(a-x,b-x, x >= xm)
          d = cgold * e
        else
          d = p/q
          u = x+d
          if (u-a < tol2 .or. b-u < tol2) d = sign(tol1,xm-x)
        end if
      else
        e = merge(a-x,b-x, x >= xm)
        d = cgold * e
      end if

      u = merge(x+d,x+sign(tol1,d), abs(d) >= tol1)
      fu = f%eval (u)

      if (fu <= fx) then
        if (u >= x) then
          a = x
        else
          b = x
        end if
        call shft3(v,w,x,u)
        call shft3(fv,fw,fx,fu)
      else
        if (u < x) then
          a = u
        else
          b = u
        end if
        if (fu <= fw .or. w==x) then
          v = w
          w = u
          fv = fw
          fw = fu
        else if (fu <= fv .or. v==x .or. v==w) then
          v = u
          fv = fu
        end if
      end if
    end do

    ! too many iterations if we get here
    write(*,*) 'too many brent iterations!'

  end function brent

  pure subroutine shft3 (a,b,c,d)
    real(r8), intent(out)   :: a
    real(r8), intent(inout) :: b,c
    real(r8), intent(in)    :: d

    a = b
    b = c
    c = d
  end subroutine shft3

end module brent_module
