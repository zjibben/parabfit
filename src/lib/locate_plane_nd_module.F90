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

  public :: locate_plane_nd, locate_plane_nd_unit_test_suite

  ! define type for error function
  type, extends(brent_func) :: volume_error_func
    real(r8)         :: tvol,norm(3)
    type(polyhedron) :: poly
  contains
    procedure        :: init => func_init
    procedure        :: eval => func_eval
    procedure        :: signed_eval => func_signed_eval
  end type volume_error_func
  
contains

  subroutine locate_plane_nd_unit_test_suite ()
    use array_utils, only: isZero
    use hex_types,   only: cube_v, hex_f, hex_e

    logical          :: success
    type(polyhedron) :: cube
    type(plane)      :: P
    real(r8)         :: vol

    ! just working with a cube for now
    ! need either more later, or tests on nested polyhedrons
    call cube%init (cube_v, hex_f, hex_e)

    write(*,*)
    write(*,*) 'LOCATE PLANE NESTED DISSECTION'
    write(*,*) '===================================================='

    ! cube half filled along x
    P%normal = [1.0_r8, 0.0_r8, 0.0_r8]
    P%rho    = 0.5_r8
    vol      = 0.5_r8
    success  = locate_plane_nd_unit_test (cube, P, vol)
    write(*,*) 'cube half filled along x?              ',success

    ! cube half filled along y
    P%normal = [0.0_r8, 1.0_r8, 0.0_r8]
    P%rho    = 0.5_r8
    vol      = 0.5_r8
    success  = locate_plane_nd_unit_test (cube, P, vol)
    write(*,*) 'cube half filled along y?              ',success

    ! cube half filled along z
    P%normal = [0.0_r8, 0.0_r8, 1.0_r8]
    P%rho    = 0.5_r8
    vol      = 0.5_r8
    success  = locate_plane_nd_unit_test (cube, P, vol)
    write(*,*) 'cube half filled along z?              ',success

    ! cube 8/10ths filled along x
    P%normal = [1.0_r8, 0.0_r8, 0.0_r8]
    P%rho    = 0.8_r8
    vol      = 0.8_r8
    success  = locate_plane_nd_unit_test (cube, P, vol)
    write(*,*) 'cube 8/10ths filled along x?           ',success
    
    ! cube 8/10ths filled along -x
    P%normal = -[1.0_r8, 0.0_r8, 0.0_r8]
    P%rho    = -0.2_r8
    vol      = 0.8_r8
    success  = locate_plane_nd_unit_test (cube, P, vol)
    write(*,*) 'cube 8/10ths filled along -x?          ',success

    ! cube half filled along xy diagonal
    P%normal = [1.0_r8, 1.0_r8, 0.0_r8] / sqrt(2.0_r8)
    P%rho    = 1.0_r8/sqrt(2.0_r8)
    vol      = 0.5_r8
    success  = locate_plane_nd_unit_test (cube, P, vol)
    write(*,*) 'cube half filled along xy diagonal?    ',success
    
    ! cube one eighth filled along -xy diagonal
    P%normal = -[1.0_r8, 1.0_r8, 0.0_r8] / sqrt(2.0_r8)
    P%rho    = -1.5_r8/sqrt(2.0_r8)
    vol      = 0.125_r8
    success  = locate_plane_nd_unit_test (cube, P, vol)
    write(*,*) 'cube 1/8th filled along -xy diagonal?  ',success

    ! cube half filled along xyz diagonal
    P%normal = [1.0_r8, 1.0_r8, 1.0_r8] / sqrt(3.0_r8)
    P%rho    = 1.5_r8/sqrt(3.0_r8)
    vol      = 0.5_r8
    success  = locate_plane_nd_unit_test (cube, P, vol)
    write(*,*) 'cube half filled along xyz diagonal?   ',success

    ! cube 1/8th filled along xyz diagonal
    P%normal = [1.0_r8, 1.0_r8, 1.0_r8] / sqrt(3.0_r8)
    ! P%rho    = 0.5_r8/sqrt(3.0_r8)
    ! vol      = 1.0_r8/48.0_r8
    P%rho    = (6.0_r8*(0.125_r8))**(1.0_r8/3.0_r8) / sqrt(3.0_r8)
    vol      = 0.125_r8
    success  = locate_plane_nd_unit_test (cube, P, vol)
    write(*,*) 'cube 1/8th filled along xyz diagonal?  ',success
    
    write(*,*) '===================================================='
    write(*,*)

  end subroutine locate_plane_nd_unit_test_suite

  logical function locate_plane_nd_unit_test (poly, P, vol)
    use consts,      only: cutvof

    type(polyhedron), intent(in) :: poly
    type(plane),      intent(in) :: P
    real(r8),         intent(in) :: vol

    type(plane) :: P_out
    
    P_out = locate_plane_nd (poly, P%normal, vol)
    locate_plane_nd_unit_test = abs(P%rho - P_out%rho) < cutvof
    
  end function locate_plane_nd_unit_test

  ! given a polyhedron, a normal, and a volume, calculate a plane constant
  ! and return a plane type
  ! Alternatively, should we return the two polyhedra? If they are already
  ! calculated, it would save resplitting the input polyhedron.
  type(plane) function locate_plane_nd (poly, norm, vol)
    use polyhedron_type
    use plane_type
    
    type(polyhedron), intent(in) :: poly
    real(r8),         intent(in) :: norm(:), vol
    
    real(r8)                 :: rho_min,rho_mid,rho_max
    type(volume_error_func)  :: volume_error
    
    ! initialize error function
    call volume_error%init (norm, poly, vol)
    
    ! get bounds for Brent's method
    call rho_bracket (rho_min, rho_mid, rho_max, norm, poly, volume_error)
    
    ! start Brent's method
    locate_plane_nd%normal = norm
    locate_plane_nd%rho = brent (rho_min, rho_mid, rho_max, volume_error)

  end function locate_plane_nd

  ! loop through all the polyhedron vertices, checking if a plane
  ! intersecting each one is above or below the target volume,
  ! thereby bracketing the allowed range of plane constants
  subroutine rho_bracket (rho_min,rho_mid,rho_max, norm, poly, volume_error)
    use logging_services
    
    real(r8),                intent(out) :: rho_min, rho_mid, rho_max
    real(r8),                intent(in)  :: norm(:)
    type(polyhedron),        intent(in)  :: poly ! we could instead just pass in poly%x
    type(volume_error_func), intent(in)  :: volume_error

    real(r8)                             :: err_min,err_max,err, rho
    integer                              :: i
    integer, parameter                   :: iter_max = 10

    err_min = -huge(1.0_r8); rho_min = -huge(1.0_r8)
    err_max =  huge(1.0_r8); rho_max =  huge(1.0_r8)

    ! find the outer bounds
    do i = 1,poly%nVerts
      rho = sum(poly%x(:,i)*norm)

      if (rho <= rho_min .or. rho >= rho_max) cycle

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

    if (i>iter_max .and. (rho_mid >= rho_max .or. rho_mid <= rho_min)) then
      write(*,*) rho_min,rho_mid,rho_max
      call LS_fatal('too many bracketing iterations!')
    end if

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
    class(volume_error_func), intent(in) :: this
    real(r8),                 intent(in) :: x
    
    func_eval = abs(this%signed_eval (x))

  end function func_eval

  real(r8) function func_signed_eval (this, x)
    use plane_type

    class(volume_error_func), intent(in) :: this
    real(r8),                 intent(in) :: x

    type(plane) :: P

    P%rho = x; P%normal = this%norm
    func_signed_eval = this%poly%volume_behind_plane (P) - this%tvol

  end function func_signed_eval
  
end module locate_plane_nd_module
