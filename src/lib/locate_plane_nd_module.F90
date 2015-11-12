!!
!! LOCATE_PLANE_ND_MODULE
!!
!! This module provides a plane reconstruction
!! subroutine for the nested dissection method.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! Last revised 4 Nov 2012.
!!

module locate_plane_nd_module
  use kinds, only: r8
  use plane_type
  use polyhedron_type
  use brent_module
  implicit none
  private

  public :: locate_plane_nd

  ! define type for error function
  type, extends(brent_func) :: volume_error_func
    real(r8)         :: tvol,norm(3)
    type(polyhedron) :: poly
  contains
    procedure :: init => func_init
    procedure :: eval => func_eval
    procedure :: signed_eval => func_signed_eval
  end type volume_error_func
  
contains

  ! given a polyhedron, a normal, and a volume, calculate a plane constant
  ! and return a plane type
  ! Alternatively, should we return the two polyhedra? If they are already
  ! calculated, it would save resplitting the input polyhedron.
  function locate_plane_nd (poly, norm, vol) result(P)
    use polyhedron_type
    use plane_type
    
    type(polyhedron), intent(in) :: poly
    real(r8),         intent(in) :: norm(:), vol
    type(plane)                  :: P
    
    real(r8)                 :: rho_min,rho_mid,rho_max
    type(volume_error_func)  :: volume_error
    
    ! initialize error function
    call volume_error%init (norm, poly, vol)

    ! get bounds for Brent's method
    call rho_bracket (rho_min, rho_mid, rho_max, norm, poly, volume_error)

    ! start Brent's method
    P%normal = norm
    P%rho = brent (rho_min, rho_mid, rho_max, volume_error)
    
  end function locate_plane_nd

  ! loop through all the polyhedron vertices, checking if a plane
  ! intersecting each one is above or below the target volume,
  ! thereby bracketing the allowed range of plane constants
  subroutine rho_bracket (rho_min,rho_mid,rho_max, norm, poly, volume_error)
    real(r8),                intent(out) :: rho_min, rho_mid, rho_max
    real(r8),                intent(in)  :: norm(:)
    type(polyhedron),        intent(in)  :: poly ! we could instead just pass in poly%x
    type(volume_error_func), intent(in)  :: volume_error

    real(r8) :: err_min,err_max,err, rho
    integer  :: i

    err_min = -huge(1.0_r8); rho_min = -huge(1.0_r8)
    err_max =  huge(1.0_r8); rho_max =  huge(1.0_r8)

    ! find the outer bounds
    do i = 1,poly%nVerts
      rho = sum(poly%x(:,i)*norm)

      if (rho < rho_min .or. rho > rho_max) cycle

      err = volume_error%signed_eval (rho)
      
      if (0.0_r8 < err .and. err < err_max) then
        err_max = err
        rho_max = rho
      else if (err_min < err .and. err < 0.0_r8) then
        err_min = err
        rho_min = rho
      end if
    end do

    ! Check if the bounds were not set? None of them are huge, right?

    ! Brent's method requires a guess minimum, find one
    rho_mid = huge(1.0_r8)
    rho = 0.5_r8 * (rho_min+rho_max) ! first check the middle
    do while (rho_mid >= rho_max .or. rho_mid <= rho_min)
      err = volume_error%eval (rho)
      if (err < abs(err_min) .and. err < abs(err_max)) then
        rho_mid = rho
      else
        if (err > 0.0_r8) then
          rho = 0.5_r8 * (rho + rho_min)
        else
          rho = 0.5_r8 * (rho + rho_max)
        end if
      end if
    end do

  end subroutine rho_bracket

  ! error function procedures

  subroutine func_init (this, norm_in,poly_in,tvol_in)
    use polyhedron_type

    class(volume_error_func), intent(out) :: this
    real(r8),                 intent(in)  :: norm_in(:), tvol_in
    type(polyhedron),         intent(in)  :: poly_in

    this%norm = norm_in
    this%poly = poly_in
    this%tvol = tvol_in
    
  end subroutine func_init

  real(r8) function func_eval (this, x)
    use plane_type

    class(volume_error_func), intent(in) :: this
    real(r8),                 intent(in) :: x

    type(plane) :: P
    
    P%rho = x; P%normal = this%norm
    func_eval = this%poly%volume_behind_plane (P) - this%tvol
  end function func_eval

  real(r8) function func_signed_eval (this, rho)
    class(volume_error_func), intent(in) :: this
    real(r8),                 intent(in) :: rho

    func_signed_eval = abs(this%eval (rho))
  end function func_signed_eval

  
  
end module locate_plane_nd_module
