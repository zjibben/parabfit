module volume_track_test

  use kinds,  only: r8
  use volume_track_module
  implicit none
  private

  public :: volume_track_unit_test_suite

contains

  subroutine volume_track_unit_test_suite ()

    use locate_plane_module, only: locate_plane_hex
    use hex_types, only: cube_v

    integer                :: iter
    type(locate_plane_hex) :: locate_plane
    real(r8)               :: posXflow(6), posxyzn(3), posxyn(3), posxn(3)
    logical                :: test_result

    write(*,*)
    write(*,*) 'VOLUME TRACK'
    write(*,*) '===================================================='

    ! define face velocities [??, ??, -x, +x, ??, ??]
    posXflow = [ 0.0_r8, 0.0_r8, -1.0_r8, 1.0_r8, 0.0_r8, 0.0_r8 ]

    ! define normals
    posxyzn = [1.0_r8, 1.0_r8, 1.0_r8] / sqrt(3.0_r8) ! positive x-y-z
    posxyn  = [1.0_r8, 1.0_r8, 0.0_r8] / sqrt(2.0_r8) ! positive x-y plane
    posxn   = [1.0_r8, 0.0_r8, 0.0_r8]                ! positive x-direction

    ! cube half filled along xyz diagonal
    ! exact rho: sqrt(3.0_r8)/2.0_r8
    call locate_plane%init (posxyzn, 0.5_r8, 1.0_r8, cube_v)
    call locate_plane%locate_plane (iter)
    
    !test_result = material_flux_unit_test (locate_plane, 0.5_r8*posXflow, 7.0_r8/48.0_r8)

    ! cube half filled along xy diagonal
    ! test_result = locate_plane_unit_test (posxyn, 0.5_r8, 1.0_r8/sqrt(2.0_r8))

    !test_result = material_flux_unit_test (locate_plane, 0.5_r8*posXflow, 0.125_r8)

    ! cube one eighth filled along xy diagonal
    ! test_result = locate_plane_unit_test (-posxyn, 0.125_r8, 1.0_r8/(2.0_r8*sqrt(2.0_r8)))

    !locate_plane%normal = - locate_plane%normal
    ! call material_flux_unit_test (locate_plane,  0.5_r8*posXflow, 0.0_r8)
    ! call material_flux_unit_test (locate_plane, -0.5_r8*posXflow, 0.125_r8)

    ! cube 8/10ths filled in x direction (exact rho: 0.8_r8)
    call locate_plane%init (posxn, 0.8_r8, 1.0_r8, cube_v)
    call locate_plane%locate_plane (iter)

    test_result = material_flux_unit_test (locate_plane,  0.25_r8*posXflow, &
        [0.0_r8, 0.0_r8, 0.0_r8, 0.05_r8, 0.0_r8, 0.0_r8])
    write(*,*) 'passed cube 8/10ths filled in x, fluxed in +x?  ', test_result

    test_result = material_flux_unit_test (locate_plane, -0.25_r8*posXflow, &
        [0.0_r8, 0.0_r8, 0.25_r8, 0.0_r8, 0.0_r8, 0.0_r8])
    write(*,*) 'passed cube 8/10ths filled in x, fluxed in -x?  ', test_result


    ! cube 8/10ths filled in -x direction (exact rho: -0.2_r8)
    call locate_plane%init (-posxn, 0.8_r8, 1.0_r8, cube_v)
    call locate_plane%locate_plane (iter)

    test_result = material_flux_unit_test (locate_plane,  0.25_r8*posXflow, &
        [0.0_r8, 0.0_r8, 0.0_r8, 0.25_r8, 0.0_r8, 0.0_r8])
    write(*,*) 'passed cube 8/10ths filled in -x, fluxed in +x? ', test_result

    test_result = material_flux_unit_test (locate_plane, -0.25_r8*posXflow, &
        [0.0_r8, 0.0_r8, 0.05_r8, 0.0_r8, 0.0_r8, 0.0_r8])
    write(*,*) 'passed cube 8/10ths filled in -x, fluxed in -x? ', test_result

    write(*,*) '===================================================='
    write(*,*)

  end subroutine volume_track_unit_test_suite

  logical function material_flux_unit_test (plane_cell, face_fluxing_velocity, fluxex)
    use kinds,       only: r8
    use hex_types,   only: cell_data
    use array_utils, only: isZero
    use locate_plane_module

    type(locate_plane_hex), intent(in) :: plane_cell
    real(r8),               intent(in) :: fluxex(:), face_fluxing_velocity(:)

    type(cell_data) :: cell
    real(r8)        :: mat_vol_flux(6), flux_vol_sum(6)

    call cell%init (plane_cell%node)

    mat_vol_flux = 0.0_r8
    flux_vol_sum = 0.0_r8

    mat_vol_flux = material_volume_flux (flux_vol_sum, plane_cell, cell, .true., face_fluxing_velocity, 1.0_r8, 0.5_r8)

    !write(*,*) 'mat_vol_flux', mat_vol_flux,'exact: ',fluxex
    material_flux_unit_test = all(isZero (mat_vol_flux-fluxex))

  end function material_flux_unit_test

end module volume_track_test
