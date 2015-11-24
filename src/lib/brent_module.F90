!!
!! BRENT_MODULE
!!
!! This module provides an interface to a Brent's algorithm.
!! Note this could probably be replaced with some canned package,
!! like the GNU Scientific Library, netlib, or something else.
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
  
  real(r8), parameter :: cgold = 1.0_r8 - 2.0_r8/(1.0_r8 + sqrt(5.0_r8))

  type, abstract, public :: brent_func
  contains
    procedure(eval), deferred :: eval
  end type brent_func

  abstract interface
    real(r8) function eval (this, x)
      import r8, brent_func
      class(brent_func), intent(in) :: this
      real(r8),          intent(in) :: x
    end function eval
  end interface

contains

  ! perform Brent's method
  function brent (f,x_min,x_mid,x_max,tol,iter_max) result(x)
    class(brent_func), intent(in) :: f
    real(r8),          intent(in) :: x_min,x_mid,x_max,tol
    integer,           intent(in) :: iter_max

    real(r8) :: a,b, fu,fv,fw,fx, u,v,w,x,xm, p,q,r,tol1,tol2, d,e
    integer  :: iter

    a = x_min; b = x_max
    x = x_mid; w = x_mid; v = x_mid
    fx = f%eval (x); fw = fx; fv = fx
    e = 0.0_r8; d = 0.0_r8

    do iter = 1,iter_max
      xm = 0.5_r8 * (a+b)
      tol1 = tol*abs(x) + alittle
      tol2 = 2.0_r8 * tol1
      if (abs(x-xm) <= tol2-0.5_r8*(b-a)) return

      ! interpolate the polynomial
      if (abs(e) > tol1) then
        r = (x-w)*(fx-fv)
        q = (x-v)*(fx-fw)
        p = sign(1.0_r8, r-q)*((x-v)*q - (x-w)*r)
        q = 2.0_r8*abs(q-r)
      end if
      
      if (abs(e) > tol1 .and. abs(p) < 0.5_r8*q*abs(e) .and. a-x < p/q .and. p/q < b-x) then
        e = d
        d = p/q
        u = x+d
        if (u-a < tol2 .or. b-u < tol2) d = sign(tol1,xm-x)
      else
        e = merge(a, b, x >= xm) - x ! distance from x to farthest bracket
        d = cgold * e
      end if

      ! evaluate at the next point
      u = x + merge(d, sign(tol1,d), abs(d) >= tol1)
      fu = f%eval (u)

      if (fu <= fx) then
        ! u is a new minima location
        ! x is now a left or right bracket
        if (u >= x) then
          a = x
        else
          b = x
        end if
        ! cycle the points
        call shft3 (v,w,x,u)
        call shft3 (fv,fw,fx,fu)
      else
        ! u is now a left or right bracket
        if (u < x) then
          a = u
        else
          b = u
        end if
        if (fu <= fw .or. w==x) then
          call shft2 (v,w,u)
          call shft2 (fv,fw,fu)
        else if (fu <= fv .or. v==x .or. v==w) then
          v = u
          fv = fu
        end if
      end if
    end do

    ! too many iterations if we get here
    write(*,*) 'too many brent iterations! error: ',fx

  end function brent
  
  pure subroutine shft3 (a,b,c,d)
    real(r8), intent(out)   :: a
    real(r8), intent(inout) :: b,c
    real(r8), intent(in)    :: d

    a = b
    b = c
    c = d
  end subroutine shft3

  pure subroutine shft2 (a,b,c)
    real(r8), intent(out)   :: a
    real(r8), intent(inout) :: b
    real(r8), intent(in)    :: c

    a = b
    b = c
  end subroutine shft2

end module brent_module
