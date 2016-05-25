module locate_plane_nd_test

  use kinds, only: r8
  use locate_plane_nd_module
  use plane_type
  use polyhedron_type
  use brent_module
  use logging_services
  implicit none
  private

  public :: locate_plane_nd_unit_test_suite

contains

  subroutine locate_plane_nd_unit_test_suite ()

    use array_utils, only: isZero
    use hex_types,   only: cube_v, hex_f, hex_e

    logical          :: success
    type(polyhedron) :: cube
    type(plane)      :: P
    real(r8)         :: vol
    integer          :: ierr

    ! just working with a cube for now
    ! need either more later, or tests on nested polyhedrons
    call cube%init (ierr, cube_v, hex_f, hex_e)

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
    use array_utils, only: isZero

    type(polyhedron), intent(inout) :: poly
    type(plane),      intent(in) :: P
    real(r8),         intent(in) :: vol

    type(plane) :: P_out

    P_out = locate_plane_nd (poly, P%normal, vol, poly%volume())
    locate_plane_nd_unit_test = isZero (P%rho - P_out%rho, cutvof)

  end function locate_plane_nd_unit_test

end module locate_plane_nd_test
