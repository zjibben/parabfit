#include "f90_assert.fpp"

module locate_plane_test

  use kinds,     only: r8
  use consts,    only: ndim,nfc,nvc,alittle
  use locate_plane_module
  use logging_services
  implicit none
  private
  
  public :: locate_plane_unit_test_suite, locate_plane_unit_test

contains

  subroutine locate_plane_unit_test_suite ()
    type(locate_plane_hex) :: locate_plane
    real(r8)               :: posxyzn(3), posxyn(3), posxn(3)
    logical                :: test_result

    write(*,*)
    write(*,*) 'LOCATE PLANE'
    write(*,*) '===================================================='

    ! define normals
    posxyzn = [1.0_r8, 1.0_r8, 1.0_r8] / sqrt(3.0_r8) ! positive x-y-z
    posxyn  = [1.0_r8, 1.0_r8, 0.0_r8] / sqrt(2.0_r8) ! positive x-y direction
    posxn   = [1.0_r8, 0.0_r8, 0.0_r8]                ! positive x-direction

    ! cube half filled along x
    test_result = locate_plane_unit_test (posxn, 0.5_r8, 0.5_r8)
    write(*,*) 'passed cube half filled along x?     ', test_result

    ! cube half filled along xy diagonal
    test_result = locate_plane_unit_test (posxyn, 0.5_r8, 1.0_r8/sqrt(2.0_r8))
    write(*,*) 'passed cube half filled along xy?    ', test_result

    ! cube half filled along xyz diagonal
    test_result = locate_plane_unit_test (posxyzn, 0.5_r8, sqrt(3.0_r8)/2.0_r8)
    write(*,*) 'passed cube half filled along xyz?   ', test_result

    ! cube 8/10ths filled in x direction
    test_result = locate_plane_unit_test (posxn, 0.8_r8, 0.8_r8)
    write(*,*) 'passed cube 8/10ths filled along x?  ', test_result

    ! cube 8/10ths filled in -x direction
    test_result = locate_plane_unit_test (-posxn, 0.8_r8, -0.2_r8)
    write(*,*) 'passed cube 8/10ths filled along -x? ', test_result

    ! cube one eighth filled along -xy diagonal
    test_result = locate_plane_unit_test (-posxyn, 0.125_r8, -1.5_r8/sqrt(2.0_r8))
    write(*,*) 'passed cube 1/8th filled along -xy?  ', test_result

    write(*,*) '===================================================='
    write(*,*)

  end subroutine locate_plane_unit_test_suite

  logical function locate_plane_unit_test (normal, vof, rhoex)

    use consts,    only: cutvof
    use hex_types, only: cube_v

    real(r8), intent(in) :: normal(:), vof, rhoex

    integer :: iter
    type(locate_plane_hex) :: lph

    ASSERT(size(normal)==ndim)

    call lph%init (normal, vof, 1.0_r8, cube_v)
    
    call lph%locate_plane (iter)
    
    !write(*,'(2(a,f14.10),L)') 'plane constant = ',lph%P%rho,',   correct plane constant = ',rhoex
    locate_plane_unit_test = abs(lph%P%rho-rhoex)<cutvof

  end function locate_plane_unit_test

end module locate_plane_test
