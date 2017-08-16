!!
!! BRENT_ROOT_CLASS
!!
!! This module provides an interface to a Brent's root-finding algorithm.
!! Note this could be replaced with some canned package,
!! like the GNU Scientific Library, netlib, or something else.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! July 2016
!!

#include "f90_assert.fpp"

module brent_root_class

  use kinds, only: r8
  use consts, only: alittle
  implicit none
  private

  type, abstract, public :: brent_root
    real(r8) :: eps
    integer  :: maxitr, numitr = 0
  contains
    procedure, non_overridable :: find_root
    procedure(func), deferred :: f
  end type brent_root

  abstract interface
    function func (this, x) result (fx)
      import brent_root, r8
      class(brent_root), intent(in) :: this
      real(r8), intent(in) :: x
      real(r8) :: fx
    end function func
  end interface

contains

  subroutine find_root (this, xmin, xmax, root, stat)

    class(brent_root), intent(inout) :: this
    real(r8), intent(in) :: xmin, xmax
    real(r8), intent(out) :: root
    integer, intent(out) :: stat

    real(r8) :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
    integer :: i

    a = xmin
    b = xmax
    c = xmax
    fa = this%f(a)
    fb = this%f(b)
    fc = fb

    stat = 0

    root=huge(1.0_r8)
    do i = 1,this%maxitr
      if ((fb > 0 .and. fc > 0) .or. (fb < 0 .and. fc < 0)) then
        c = a
        fc = fa
        e = b-a
        d = e
      end if

      if (abs(fc) < abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
      end if

      tol1 = 2*alittle*abs(b) + this%eps / 2
      xm = (c-b) / 2

      if (abs(xm) <= tol1 .or. fb==0) then
        root=b
        this%numitr = i
        return
      end if

      if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
        s = fb/fa
        if (a==c) then
          p = 2 * xm * s
          q = 1.0 - s
        else
          q = fa / fc
          r = fb / fc
          p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1))
          q = (q-1)*(r-1)*(s-1)
        end if
        if (p > 0) q = -q
        p = abs(p)
        if (2*p < min(3*xm*q-abs(tol1*q), abs(e*q))) then
          e=d
          d=p/q
        else
          d=xm
          e=d
        end if
      else
        d=xm
        e=d
      end if
      a=b
      fa=fb
      if (abs(d) > tol1) then
        b = b + d
      else
        b = b + sign(tol1,xm)
      end if
      fb = this%f(b)
    end do

    stat = 1

  end subroutine find_root

end module brent_root_class
