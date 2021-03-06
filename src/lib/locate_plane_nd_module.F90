!!
!! LOCATE_PLANE_ND_MODULE
!!
!! This module provides a plane reconstruction
!! subroutine for the nested dissection method.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! March 2016
!!

#include "f90_assert.fpp"

module locate_plane_nd_module

  use kinds, only: r8
  use consts, only: ndim
  use plane_type
  use polyhedron_type
  use brent_root_class
  use logging_services
  implicit none
  private

  public :: locate_plane_nd

  ! define type for error function
  type, extends(brent_root) :: vof_error_func
    real(r8)         :: tvol,norm(ndim),parvol
    type(polyhedron) :: poly
  contains
    procedure :: init => func_init
    procedure :: eval => func_eval
    procedure :: signed_eval => func_signed_eval
    procedure :: f => func_signed_eval
    !procedure :: f => func_eval
  end type vof_error_func

contains

  ! given a polyhedron, a normal, and a volume, calculate a plane constant
  ! and return a plane type
  ! Alternatively, should we return the two polyhedra? If they are already
  ! calculated, it would save resplitting the input polyhedron.
  type(plane) function locate_plane_nd (poly, norm, vol, cell_volume)

    use consts, only: cutvof,ndim
    use polyhedron_type
    use plane_type

    type(polyhedron), intent(inout) :: poly
    real(r8),         intent(in) :: norm(:), vol, cell_volume

    real(r8)             :: rho_min,rho_mid,rho_max
    integer              :: ierr
    type(vof_error_func) :: vof_error

    ASSERT(size(norm)==ndim)

    ! initialize error function
    call vof_error%init (norm, poly, vol, cell_volume)

    ! get bounds for Brent's method
    call rho_bracket (rho_min, rho_mid, rho_max, norm, poly, vof_error)

    ! start Brent's method (note ~30 iterations seem to be necessary to pass current unit tests)
    locate_plane_nd%normal = norm
    vof_error%eps = cutvof / 1e10_r8; vof_error%maxitr = 50
    call vof_error%find_root (rho_min, rho_max, locate_plane_nd%rho, ierr)

  end function locate_plane_nd

  ! loop through all the polyhedron vertices, checking if a plane
  ! intersecting each one is above or below the target volume,
  ! thereby bracketing the allowed range of plane constants
  subroutine rho_bracket (rho_min,rho_mid,rho_max, norm, poly, volume_error)
    use logging_services

    real(r8),                intent(out) :: rho_min, rho_mid, rho_max
    real(r8),                intent(in)  :: norm(:)
    type(polyhedron),        intent(inout)  :: poly
    type(vof_error_func), intent(inout)  :: volume_error

    real(r8)                             :: err_min,err_max,err, rho
    integer                              :: i, ierr
    integer, parameter                   :: iter_max = 10

    type(plane) :: P

    err_min = -huge(1.0_r8); rho_min = -huge(1.0_r8)
    err_max =  huge(1.0_r8); rho_max =  huge(1.0_r8)

    ! find the outer bounds
    do i = 1,poly%parent%nVerts
      rho = dot_product(poly%parent%x(:,i),norm)

      if (rho <= rho_min .or. rho >= rho_max) cycle

      err = volume_error%signed_eval (rho)

      if (0.0_r8 < err .and. err < err_max) then
        err_max = err
        rho_max = rho
      else if (err_min < err .and. err < 0.0_r8) then
        err_min = err
        rho_min = rho
      else if (err == 0) then
        err_min = err
        err_max = err

        rho_min = rho
        rho_max = rho
        rho_mid = rho
        return
      end if
    end do

    ! make sure the bounds were set
    if (rho_min==-huge(1.0_r8) .or. rho_max==huge(1.0_r8)) then
      P%normal = norm
      print '(a,3f10.3)', 'normal: ', norm
      call poly%print_data ()
      !write(*,*) 'volume ',poly%volume ()
      do i = 1,poly%parent%nVerts
        rho = dot_product(poly%parent%x(:,i),norm)
        err = volume_error%signed_eval(rho)
        P%rho = rho
        write(*,'(a,i3,3es14.4)') 'vertex, rho, error, vol: ', i, rho, err, &
            poly%volume_behind_plane (P,ierr)
      end do
      write(*,'(a,2es14.4)') 'rho min,max: ', rho_min, rho_max
      call LS_fatal ('rho bounds not set')
    end if

    ! Brent's method requires a guess minimum, find one
    ! first check a middle point found by weighting by the inverse of the error magnitudes
    rho_mid = (rho_min/abs(err_min)+rho_max/abs(err_max))/(1.0_r8/abs(err_min) + 1.0_r8/abs(err_max))
    err     = huge(1.0_r8)
    i = 1
    do while (i<=iter_max .and. &
         ((rho_mid >= rho_max .or. rho_mid <= rho_min) .and. &
         (abs(err) > abs(err_min) .or. abs(err) > abs(err_max))))
      err = volume_error%signed_eval (rho)
      rho_mid = merge(&
           (rho_mid/abs(err)+rho_min/abs(err_min))/(1.0_r8/abs(err) + 1.0_r8/abs(err_min)),&
           (rho_mid/abs(err)+rho_max/abs(err_max))/(1.0_r8/abs(err) + 1.0_r8/abs(err_max)),&
           err > 0.0_r8)
      i = i+1
    end do

    ! make sure the midpoint was set
    if (i>iter_max .and. (rho_mid >= rho_max .or. rho_mid <= rho_min)) then
      write(*,*) rho_min,rho_mid,rho_max
      call LS_fatal('too many bracketing iterations!')
    end if

  end subroutine rho_bracket

  ! error function procedures

  subroutine func_init (this, norm,poly,tvol,parvol)
    use polyhedron_type

    class(vof_error_func), intent(out) :: this
    real(r8),                 intent(in)  :: norm(:), tvol, parvol
    type(polyhedron),         intent(in)  :: poly

    this%norm   = norm
    this%poly   = poly
    this%tvol   = tvol
    this%parvol = parvol

  end subroutine func_init

  real(r8) function func_eval (this, x)
    class(vof_error_func), intent(inout) :: this
    real(r8), intent(in) :: x
    func_eval = abs(this%signed_eval (x))
  end function func_eval

  real(r8) function func_signed_eval (this, x)

    use plane_type

    class(vof_error_func), intent(inout) :: this
    real(r8), intent(in) :: x

    type(plane) :: P
    integer :: ierr

    P%rho = x; P%normal = this%norm
    func_signed_eval = (this%poly%volume_behind_plane (P,ierr) - this%tvol) / this%parvol

    if (ierr /= 0) call LS_fatal ("func_signed_eval failed")

  end function func_signed_eval

end module locate_plane_nd_module
