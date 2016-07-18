!!
!! BRENT_MODULE
!!
!! This module provides an interface to Brent's minimization algorithm.
!! Note this could be replaced with some canned package,
!! like the GNU Scientific Library, netlib, or something else.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! March 2016
!!

module brent_min_class

  use kinds,  only: r8
  use consts, only: alittle
  implicit none
  private

  type, abstract, public :: brent_min
    real(r8) :: eps
    integer  :: maxitr, numitr=0
  contains
    procedure, non_overridable :: find_minimum
    procedure(func), deferred :: f
  end type brent_min

  abstract interface
    function func (this, x) result(fx)
      import r8, brent_min
      class(brent_min), intent(in) :: this
      real(r8), intent(in) :: x
      real(r8) :: fx
    end function func
  end interface

contains

  ! perform Brent's method
  subroutine find_minimum (this,x_min,x_mid,x_max, x, stat)

    use logging_services
    
    class(brent_min), intent(in)  :: this
    real(r8),         intent(in)  :: x_min,x_mid,x_max
    real(r8),         intent(out) :: x
    integer,          intent(out) :: stat

    real(r8), parameter :: cgold = 1.0_r8 - 2.0_r8/(1.0_r8 + sqrt(5.0_r8))
    real(r8) :: a,b, fu,fv,fw,fx, u,v,w,xm, p,q,r,tol1,tol2, d,e
    integer  :: iter

    a = x_min; b = x_max
    x = x_mid; w = x_mid; v = x_mid
    fx = this%f (x); fw = fx; fv = fx
    e = 0.0_r8; d = 0.0_r8

    stat = 0

    do iter = 1,this%maxitr
      xm = 0.5_r8 * (a+b)
      tol1 = this%eps*abs(x) + alittle
      tol2 = 2.0_r8 * tol1
      if (abs(x-xm) <= tol2-0.5_r8*(b-a) .and. abs(fx) <= this%eps) return !.and. abs(fx) <= this%eps

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
      fu = this%f (u)

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

    ! ! too many iterations if we get here
    ! write(*,'(a,2es14.4)') 'error:  ',fx,x
    ! write(*,'(a,3es14.4)') 'bounds: ',x_min,x_mid,x_max
    ! call LS_fatal('too many brent iterations!')
    !write(*,*) 'WARNING: did not converge brent iterations'
    stat = 1

  end subroutine find_minimum
  
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

end module brent_min_class
