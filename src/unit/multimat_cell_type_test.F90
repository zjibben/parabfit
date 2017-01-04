module multimat_cell_type_test

  use kinds, only: r8
  use multimat_cell_type
  use logging_services
  implicit none
  private

  public :: multimat_cell_unit_test_suite

contains

  subroutine multimat_cell_unit_test_suite ()

    use hex_types, only: cube_v, hex_f, hex_e

    type(multimat_cell) :: cube
    logical             :: success
    integer             :: ierr
    real(r8)            :: posXflow(6), posZflow(6), posxyzn(3), posxyn(3), posxn(3), &
        posyn(3), tmp
    real(r8), allocatable :: vof(:), intnorm(:,:)


    write(*,*)
    write(*,*) 'MULTIMAT CELL TYPE'
    write(*,*) '===================================================='

    ! partitioning
    write(*,*) 'CELL PARTITIONING'
    ! cube with 3-matls, one very small
    call cube%init (ierr, reshape([&
        2.8125E-01_r8,    5.6250E-01_r8,    2.8125E-01_r8,&
        3.1250E-01_r8,    5.6250E-01_r8,    2.8125E-01_r8,&
        3.1250E-01_r8,    5.9375E-01_r8,    2.8125E-01_r8,&
        2.8125E-01_r8,    5.9375E-01_r8,    2.8125E-01_r8,&
        2.8125E-01_r8,    5.6250E-01_r8,    3.1250E-01_r8,&
        3.1250E-01_r8,    5.6250E-01_r8,    3.1250E-01_r8,&
        3.1250E-01_r8,    5.9375E-01_r8,    3.1250E-01_r8,&
        2.8125E-01_r8,    5.9375E-01_r8,    3.1250E-01_r8],&
        [3,8]), hex_f, hex_e)
    vof = [8.7485114822E-01_r8,    1.2514879878E-01_r8,    5.2994441145E-08_r8]
    intnorm = reshape([&
        9.8094153325E-01_r8,   -1.1866987389E-01_r8,   -1.5385437718E-01_r8,&
        -9.8094153821E-01_r8,    1.1866966066E-01_r8,    1.5385451003E-01_r8,&
        -4.6914971675E-01_r8,    7.9578595328E-01_r8,   -3.8291416771E-01_r8], [3,3])

    !call cube%partition (vof, intnorm)
    success = partition_unit_test (cube, vof, intnorm)
    ! success = fluxing_unit_test (cube, 2.692080062668622E-003*2.2845404866E-01_r8*posZflow, &
    !     reshape([&
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8, &
    !      0.0_r8, 0.0_r8, 0.0_r8], [3,6]) )

    write(*,*) '3-matl cube? ', success

    ! 3-matl cube2
    call cube%init (ierr, reshape([&
        2.96875E-01_r8,    3.75000E-01_r8,    5.31250E-01_r8,&
        3.12500E-01_r8,    3.75000E-01_r8,    5.31250E-01_r8,&
        3.12500E-01_r8,    3.90625E-01_r8,    5.31250E-01_r8,&
        2.96875E-01_r8,    3.90625E-01_r8,    5.31250E-01_r8,&
        2.96875E-01_r8,    3.75000E-01_r8,    5.46875E-01_r8,&
        3.12500E-01_r8,    3.75000E-01_r8,    5.46875E-01_r8,&
        3.12500E-01_r8,    3.90625E-01_r8,    5.46875E-01_r8,&
        2.96875E-01_r8,    3.90625E-01_r8,    5.46875E-01_r8], [3,8]), hex_f, hex_e)
    vof = [0.811404648085707_r8,  0.188595062059073_r8,   2.898552200505922E-007_r8]
    intnorm = reshape([&
        3.5592178718E-01_r8,   -9.3450685990E-01_r8,    4.0755626314E-03_r8,&
        -3.5592493730E-01_r8,    9.3450564851E-01_r8,   -4.0782237116E-03_r8,&
        5.7303886156E-01_r8,    5.7634688730E-01_r8,    5.8262400280E-01_r8], [3,3])
    success = partition_unit_test (cube, vof, intnorm)
    write(*,*) '3-matl cube2? ', success


    ! fluxing
    write(*,*) 'FLUXING'
    call cube%init (ierr, cube_v, hex_f, hex_e)

    ! define face velocities [+y, -y, -x, +x, -z, +z]
    posXflow = [ 0.0_r8, 0.0_r8, -1.0_r8, 1.0_r8,  0.0_r8, 0.0_r8 ]
    posZflow = [ 0.0_r8, 0.0_r8,  0.0_r8, 0.0_r8, -1.0_r8, 1.0_r8 ]

    ! define normals
    posxyzn = [1.0_r8, 1.0_r8, 1.0_r8] / sqrt(3.0_r8) ! positive x-y-z
    posxn   = [1.0_r8, 0.0_r8, 0.0_r8]                ! positive x-direction
    posyn   = [0.0_r8, 1.0_r8, 0.0_r8]                ! positive y-direction
    posxyn  = (posxn + posyn) / sqrt(2.0_r8) ! positive x-y plane

    ! cube 8/10ths filled in x direction
    deallocate(intnorm)
    allocate(intnorm(3,2))
    vof = [0.8_r8, 0.2_r8]
    intnorm(:,1) = posxn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)

    success = fluxing_unit_test (cube, 0.25_r8*posXflow, reshape([&
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.05_r8, 0.2_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 8/10ths filled in x, fluxed in +x?  ', success

    success = fluxing_unit_test (cube, -0.25_r8*posXflow, reshape([&
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.25_r8, 0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 8/10ths filled in x, fluxed in -x?  ', success

    ! cube 8/10ths filled in -x direction
    vof = [0.8_r8, 0.2_r8]
    intnorm(:,1) = -posxn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)

    success = fluxing_unit_test (cube, 0.25_r8*posXflow, reshape([&
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.25_r8, 0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 8/10ths filled in -x, fluxed in +x? ', success

    success = fluxing_unit_test (cube, -0.25_r8*posXflow, reshape([&
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.05_r8, 0.2_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 8/10ths filled in -x, fluxed in -x? ', success

    ! cube 1/8th filled in y direction
    vof = [0.125_r8, 0.875_r8]
    intnorm(:,1) = posyn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)

    success = fluxing_unit_test (cube, 0.5_r8*posXflow, reshape([&
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0625_r8,  0.4375_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 1/8th filled in y, fluxed in +x?    ', success

    ! cube 1/8th filled along xy diagonal
    vof = [0.125_r8, 0.875_r8]
    intnorm(:,1) = posxyn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)

    success = fluxing_unit_test (cube, 0.5_r8*posXflow, reshape([&
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.5_r8, &
        0.0_r8,  0.0_r8, &
        0.0_r8,  0.0_r8], [2,6]) )
    write(*,*) 'passed cube 1/8th filled in xy, fluxed in +x?   ', success

    success = fluxing_unit_test (cube, -0.5_r8*posXflow, reshape([&
        0.0_r8,   0.0_r8,   &
        0.0_r8,   0.0_r8,   &
        0.125_r8, 0.375_r8, &
        0.0_r8,   0.0_r8,   &
        0.0_r8,   0.0_r8,   &
        0.0_r8,   0.0_r8], [2,6]) )
    write(*,*) 'passed cube 1/8th filled in xy, fluxed in -x?   ', success

    ! cube 1/8th filled along xyz diagonal
    vof = [0.125_r8, 0.875_r8]
    intnorm(:,1) = posxyzn; intnorm(:,2) = -intnorm(:,1)
    call cube%partition (vof, intnorm)

    tmp = (0.25_r8 - (1.0_r8 - (6.0_r8*(0.125_r8))**(1.0_r8/3.0_r8)))**3.0_r8/6.0_r8
    success = fluxing_unit_test (cube, 0.25_r8*posXflow, reshape([&
        0.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, &
        tmp,    0.25_r8 - tmp, &
        0.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8], [2,6]) )
    write(*,*) 'passed cube 1/8th filled in xyz, fluxed in +x?  ', success

    write(*,*) '===================================================='
    write(*,*)

  end subroutine multimat_cell_unit_test_suite

  logical function fluxing_unit_test (cell, fluxing_velocity, volflux_ex)

    use consts,      only: cutvof,nfc
    use array_utils, only: isZero

    type(multimat_cell), intent(inout) :: cell
    real(r8),            intent(in)    :: fluxing_velocity(:),volflux_ex(:,:)

    real(r8) :: outward_volflux(size(volflux_ex,dim=1),nfc)
    integer  :: nmat, ierr

    nmat = size(volflux_ex,dim=1)

    write(*,*) 'WARNING: hardcoding face areas'
    outward_volflux = cell%outward_volflux (1.0_r8, fluxing_velocity, &
        [1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8],ierr)
    if (ierr /= 0) then
      call LS_fatal ("outward volflux failed")
    end if
    fluxing_unit_test = all(isZero (outward_volflux-volflux_ex,cutvof))

    ! if (.not.fluxing_unit_test) then
    !   write(*,*) 'failed: '
    !   write(*,'(6es20.10)') outward_volflux
    !   write(*,'(6es20.10)') volflux_ex
    ! end if

  end function fluxing_unit_test

  logical function partition_unit_test (cell, vof, intnorm)

    use consts, only: cutvof
    use array_utils, only: isZero

    type(multimat_cell), intent(inout) :: cell
    real(r8),            intent(in)    :: vof(:), intnorm(:,:)

    real(r8) :: vof_result
    integer :: m

    call cell%partition (vof, intnorm)

    partition_unit_test = .true.
    do m = 1,size(vof)
      vof_result = cell%mat_poly(m)%volume() / cell%volume()
      write(*,*) 'mv', m, vof_result, vof(m)
      if (.not.isZero(vof(m) - vof_result, cutvof)) then
        write(*,*) 'partition test failed on material: ',m
        write(*,'(a,4es20.10)') 'vof, vof_result, error: ',vof(m), vof_result, abs(vof(m)-vof_result)
        write(*,'(a,4es20.10)') 'remaining vof: ',1.0_r8 - sum(vof(:m))
        write(*,*) 'failure may simply be due to insufficient number of iterations'
        partition_unit_test = .false.
      end if
    end do

  end function partition_unit_test

end module multimat_cell_type_test
