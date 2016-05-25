module polyhedron_type_test

  use kinds, only: r8
  use polyhedron_type
  implicit none
  private

  public :: polyhedron_unit_test

contains

  subroutine polyhedron_unit_test ()

    use plane_type
    use polygon_type
    use hex_types,   only: hex_f, hex_e, cube_v
    use array_utils, only: isZero

    type(polyhedron) :: cube,pyramid,pyramid3,cutcube,tmp(2)
    type(polygon)    :: intpoly
    type(plane)      :: P
    real(r8)         :: volume,tmpr1,tmpr2
    logical          :: success
    integer          :: ierr
    real(r8)         :: pyr3_v(3,4) = reshape([ &
        0.0_r8, 0.0_r8, 0.0_r8, & ! vertex positions
        0.0_r8, 1.0_r8, 0.0_r8, &
        1.0_r8, 0.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 1.0_r8],&
        shape(pyr3_v))
    integer          :: pyr3_f(3,4) = reshape([ &
        1,2,3, & ! face vertices
        1,4,2, &
        1,3,4, &
        4,3,2],&
        shape(pyr3_f))
    integer          :: pyr3_e(2,6) = reshape([ &
        1,2, & ! edge vertices
        2,3, &
        3,4, &
        1,3, &
        4,1, &
        4,2],&
        shape(pyr3_e))
    real(r8)         :: pyr_v(3,5) = reshape([ & ! pyramid with square bottom
        0.0_r8, 0.0_r8, 0.0_r8, & ! vertex positions
        1.0_r8, 0.0_r8, 0.0_r8, &
        1.0_r8, 1.0_r8, 0.0_r8, &
        0.0_r8, 1.0_r8, 0.0_r8, &
        0.5_r8, 0.5_r8, 0.5_r8],&
        shape(pyr_v))
    integer          :: pyr_f(4,5) = reshape([ &
        1,2,3,4, & ! face vertices
        1,2,5,0, &
        2,3,5,0, &
        3,4,5,0, &
        4,1,5,0],&
        shape(pyr_f))
    integer          :: pyr_e(2,8) = reshape([ &
        1,2, & ! edge vertices
        2,3, &
        3,4, &
        4,1, &
        1,5, &
        2,5, &
        3,5, &
        4,5],&
        shape(pyr_e))
    real(r8)         :: cutcube_v(3,10) = reshape([ &
        0.9_r8, 0.0_r8, 0.0_r8, & ! cube cut with triangle face
        0.0_r8, 0.9_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 0.9_r8, &
        1.0_r8, 0.0_r8, 0.0_r8, &
        1.0_r8, 1.0_r8, 0.0_r8, &
        0.0_r8, 1.0_r8, 0.0_r8, &
        0.0_r8, 0.0_r8, 1.0_r8, &
        1.0_r8, 0.0_r8, 1.0_r8, &
        1.0_r8, 1.0_r8, 1.0_r8, &
        0.0_r8, 1.0_r8, 1.0_r8],&
        shape(cutcube_v))
    integer :: cutcube_f(5,7) = reshape([ &
        3,2,1,0,0, &
        7,3,1,4,8, &
        4,5,9,8,0, &
        10,9,5,6,0, &
        10,6,2,3,7, &
        1,2,6,5,4, &
        10,7,8,9,0],&
        shape(cutcube_f))
    integer :: cutcube_e(2,15) = reshape([ &
        1,2, &
        2,3, &
        3,1, &
        1,4, &
        4,5, &
        5,6, &
        6,2, &
        3,7, &
        7,8, &
        8,9, &
        10,9, &
        10,7, &
        10,6, &
        9,5, &
        8,4],  &
        shape(cutcube_e))

    write(*,*)
    write(*,*) 'POLYHEDRON'
    write(*,*) '===================================================='

    ! calculate the volume of a unit cube (1)
    write(*,*) 'SHAPE VOLUMES'
    call cube%init (ierr, cube_v, hex_f, hex_e)
    volume = cube%volume ()
    success = isZero (volume-1.0_r8)
    write(*,*) 'cube volume?                     ', success
    if (.not.success) write(*,*) 'volume: ',volume

    ! calculate the volume of a pyramid (1/6)
    call pyramid%init (ierr, pyr_v, pyr_f, pyr_e)
    volume = pyramid%volume ()
    success = isZero (volume-1.0_r8/6.0_r8)
    write(*,*) 'pyramid volume?                  ', success
    if (.not.success) write(*,*) 'volume: ',volume

    ! calculate the volume of a pyramid
    call pyramid%init (ierr, reshape([&
        3.7500000000000000000000000e-01_r8, 3.5937500000000000000000000e-01_r8, 2.0312500000000000000000000e-01_r8,&
        3.7500000000000000000000000e-01_r8, 3.7187838412018936473657504e-01_r8, 2.0312500000000000000000000e-01_r8,&
        3.7500000000000000000000000e-01_r8, 3.5937500000000000000000000e-01_r8, 2.0205613718049295068901472e-01_r8,&
        3.8156580412702972848748573e-01_r8, 3.5937500000000000000000000e-01_r8, 2.0312500000000000000000000e-01_r8], [3,4]), reshape([&
        4, 1, 3,&
        3, 1, 2,&
        2, 1, 4,&
        4, 3, 2], [3,4]), reshape([&
        1, 3,&
        1, 4,&
        1, 2,&
        2, 3,&
        3, 4,&
        4, 2], [2,6]))
    volume = pyramid%volume ()
    success = volume > 0.0_r8
    write(*,*) 'pyramid2 volume?                 ', success
    if (.not.success) write(*,*) 'volume: ',volume

    ! calculate the volume of a pyramid
    call pyramid%init (ierr, reshape([&
        2.5000000000000000000000000e-01_r8, 2.5000000000000000000000000e-01_r8, 0.0000000000000000000000000e+00_r8,&
        2.5000000000000000000000000e-01_r8, 2.6287194342498776400418592e-01_r8, 0.0000000000000000000000000e+00_r8,&
        2.5000000000000000000000000e-01_r8, 2.5000000000000000000000000e-01_r8, -1.1618077816026584070385752e-04_r8,&
        2.6287373171168715302314922e-01_r8, 2.5000000000000000000000000e-01_r8, 0.0000000000000000000000000e+00_r8], [3,4]), reshape([&
        4, 1, 3,&
        3, 1, 2,&
        2, 1, 4,&
        4, 3, 2], [3,4]), reshape([&
        1, 3,&
        1, 4,&
        1, 2,&
        2, 3,&
        3, 4,&
        4, 2], [2,6]))
    volume = pyramid%volume ()
    success = volume > 0.0_r8
    write(*,*) 'pyramid3 volume?                 ', success
    if (.not.success) write(*,*) 'volume: ',volume

    ! calculate the volume of a pyramid (1/6)
    call pyramid3%init (ierr, pyr3_v, pyr3_f, pyr3_e)
    volume = pyramid3%volume ()
    write(*,*) 'pyramid3 volume?                 ', isZero (volume-1.0_r8/6.0_r8)

    ! calculate the volume of a "cutcube" (1-0.9**3/6)
    call cutcube%init (ierr, cutcube_v, cutcube_f, cutcube_e)
    volume = cutcube%volume ()
    write(*,*) 'cutcube volume?                  ', isZero (volume-(1.0_r8-0.9_r8**3/6.0_r8))

    ! calculate the volume of a "cutcube"
    call cutcube%init (ierr, reshape([&
        7.03125E-01_r8,    2.50000E-01_r8,    6.3092397329198735000943543E-01_r8,&
        6.87500E-01_r8,    2.65625E-01_r8,    6.3835968598238201909822465E-01_r8,&
        7.03125E-01_r8,    2.65625E-01_r8,    6.3231363915688509891310787E-01_r8,&
        7.03125E-01_r8,    2.65625E-01_r8,    6.3231362814843028452571616E-01_r8,&
        6.87500E-01_r8,    2.65625E-01_r8,    6.3835968087108196922230263E-01_r8,&
        6.87500E-01_r8,    2.50000E-01_r8,    6.3697002011748427019455221E-01_r8,&
        7.03125E-01_r8,    2.50000E-01_r8,    6.3092396677678574956615876E-01_r8], [3,7]), reshape([&
        5,   2,   3,   4,&
        7,   1,   6,   0,&
        6,   2,   5,   0,&
        4,   3,   1,   7,&
        3,   2,   6,   1,&
        7,   6,   5,   4], [4,6]), reshape([&
        1,   7,&
        3,   4,&
        2,   5,&
        1,   6,&
        2,   6,&
        2,   3,&
        3,   1,&
        4,   5,&
        5,   6,&
        6,   7,&
        7,   4], [2,11]))
    volume = cutcube%volume ()
    write(*,*) 'cutcube2 volume?                 ', volume > 0.0_r8

    ! calculate the volume of a "cutcube"
    call cutcube%init (ierr, reshape([&
        6.2500000000000000000000000E-01_r8,    4.9205247210652008904574473E-01_r8,    2.6562500000000000000000000E-01_r8,&
        6.2500000000000000000000000E-01_r8,    4.8437500000000000000000000E-01_r8,    2.5576988559870150741204498E-01_r8,&
        6.2500000000000000000000000E-01_r8,    4.9205228021985458752851628E-01_r8,    2.6562500000000000000000000E-01_r8,&
        6.2500000000000000000000000E-01_r8,    4.8437500000000000000000000E-01_r8,    2.5577007468008733370723462E-01_r8,&
        6.0962392793390562939492838E-01_r8,    4.8437500000000000000000000E-01_r8,    2.6523569565858362562238426E-01_r8,&
        6.0937500000000000000000000E-01_r8,    4.8455889840500104837062167E-01_r8,    2.6562500000000000000000000E-01_r8], [3,6]), reshape([&
        5,   2,   4,   0,&
        4,   2,   1,   3,&
        3,   1,   6,   0,&
        6,   1,   2,   5,&
        6,   5,   4,   3], [4,5]), reshape([&
        2,   4,&
        1,   3,&
        1,   6,&
        2,   5,&
        2,   1,&
        3,   4,&
        4,   5,&
        5,   6,&
        6,   3], [2,9]))
    volume = cutcube%volume ()
    write(*,*) 'cutcube3 volume?                 ', volume > 0.0_r8

    ! create a plane, and return coordinates it intersects with polyhedron edges
    write(*,*) 'SHAPE SPLITTING'
    P%normal = [ 1.0_r8, 0.0_r8, 0.0_r8 ]
    P%rho    = 0.5_r8

    ! intpoly = cube%intersection_verts (P)

    ! write(*,*) 'intersection points'
    ! do i = 1,intpoly%nVerts
    !   write(*,*) intpoly%x(:,i)
    ! end do

    ! split the cube vertically down the center
    !write(*,*) 'cube split volumes'
    call cube%split (P,tmp,ierr)
    success = isZero (tmp(1)%volume ()-0.5_r8) .and. isZero (tmp(2)%volume ()-0.5_r8)
    write(*,*) 'vertical cut?                    ',success

    ! split the cube at an angle
    P%normal = [ 1.0_r8, 1.0_r8, 0.0_r8 ] / sqrt(2.0_r8)
    P%rho    = 1.5_r8 / sqrt(2.0_r8)
    call cube%split (P,tmp,ierr)
    success = isZero (tmp(1)%volume ()-0.125_r8) .and. isZero (tmp(2)%volume ()-0.875_r8)
    write(*,*) 'xy-angle off-center cut?         ',success

    ! split the cube at an angle through the center
    P%normal = [ 1.0_r8, 1.0_r8, 0.0_r8 ] / sqrt(2.0_r8)
    P%rho    = 1.0_r8 / sqrt(2.0_r8)
    call cube%split (P,tmp,ierr)
    success = isZero (tmp(1)%volume ()-0.5_r8) .and. isZero (tmp(2)%volume ()-0.5_r8)
    write(*,*) 'center xy-angle cut?             ',success

    ! split the cube at an angle through the center
    P%normal = [ 1.0_r8, 1.0_r8, 1.0_r8 ] / sqrt(3.0_r8)
    P%rho    = 1.5_r8 / sqrt(3.0_r8)
    call cube%split (P,tmp,ierr)
    success = isZero (tmp(1)%volume ()-0.5_r8) .and. isZero (tmp(2)%volume ()-0.5_r8)
    write(*,*) 'center xyz-angle cut?            ',success

    ! split the cube at an angle through an offset
    P%normal = [ 1.0_r8, 1.0_r8, 1.0_r8 ] / sqrt(3.0_r8)
    P%rho    = 0.5_r8/sqrt(3.0_r8)
    call cube%split (P,tmp,ierr)
    success = isZero (tmp(1)%volume ()-47.0_r8/48.0_r8) .and. isZero (tmp(2)%volume ()-1.0_r8/48.0_r8)
    write(*,*) 'off-center xyz-angle cut?        ',success

    ! split the pyramid in the x direction
    P%normal = -[ 1.0_r8, 0.0_r8, 0.0_r8 ]
    P%rho    = -0.8_r8
    call pyramid3%split (P,tmp,ierr)
    success = isZero (tmp(1)%volume ()-0.992_r8/6.0_r8) .and. isZero (tmp(2)%volume ()-4e-3_r8/3.0_r8)
    write(*,*) 'pyramid3 cut?                    ',success

    ! split the cutcube
    P%normal = -[ 1.0_r8, 0.0_r8, 0.0_r8 ]
    P%rho    = -0.8_r8
    call cutcube%split (P,tmp,ierr)
    tmpr1 = 1.0_r8-0.9_r8**3/6.0_r8
    tmpr2 = 0.2_r8 - 0.1_r8**3/6.0_r8
    success = isZero (tmp(1)%volume ()-(tmpr1-tmpr2)) .and. isZero (tmp(2)%volume ()-tmpr2)
    write(*,*) 'cutcube cut?                     ',success

    ! split the cutcube2
    P%normal = [0.0_r8, -1.0_r8, 0.0_r8]
    P%rho    = -3.43734026719638874336E-01_r8
    call cutcube%init (ierr, reshape([&
        5.7537989065754502338023713E-01_r8,    3.28125E-01_r8,    5.00000E-01_r8,&
        5.7517531782297459663766404E-01_r8,    3.28125E-01_r8,    4.84375E-01_r8,&
        5.7213201658870560528669102E-01_r8,    3.43750E-01_r8,    4.84375E-01_r8,&
        5.7233658942327592100696165E-01_r8,    3.43750E-01_r8,    5.00000E-01_r8,&
        5.7213201762257304139325242E-01_r8,    3.43750E-01_r8,    4.84375E-01_r8,&
        5.7517531975197921934039869E-01_r8,    3.28125E-01_r8,    4.84375E-01_r8,&
        5.7537989252387589100834475E-01_r8,    3.28125E-01_r8,    5.00000E-01_r8], [3,7]), reshape([&
        3,   4,   5,   0,&
        1,   2,   6,   7,&
        2,   3,   5,   6,&
        1,   7,   4,   0,&
        3,   2,   1,   4,&
        7,   6,   5,   4], [4,6]), reshape([&
        2,   6,&
        3,   5,&
        1,   7,&
        1,   2,&
        2,   3,&
        3,   4,&
        1,   4,&
        4,   5,&
        5,   6,&
        6,   7,&
        7,   4], [2,11]))
    call cutcube%split (P,tmp,ierr)
    tmpr1 = 0.0_r8
    tmpr2 = 0.0_r8
    success = ierr==0
    write(*,*) 'cutcube2 cut?                     ',success,tmp(1)%volume(),tmp(2)%volume()

    ! split the cutcube3
    P%normal = [-1.1282635566E-01_r8, -5.757662896250E-01_r8, 8.097921913669E-01_r8]
    P%rho    = 9.69687230833222828241E-03_r8
    call cutcube%init (ierr, reshape([&
        3.1250000000000E-01_r8,    4.53125000000000E-01_r8,    3.75000000000000E-01_r8,&
        3.2812500000000E-01_r8,    4.53125000000000E-01_r8,    3.75000000000000E-01_r8,&
        3.2812500000000E-01_r8,    4.68750000000000E-01_r8,    3.75000000000000E-01_r8,&
        3.1250000000000E-01_r8,    4.68750000000000E-01_r8,    3.75000000000000E-01_r8,&
        3.2812500000000E-01_r8,    4.68750000000000E-01_r8,    3.90625000000000E-01_r8,&
        3.2812500000000E-01_r8,    4.68257797705384E-01_r8,    3.90625000000000E-01_r8,&
        3.2812500000000E-01_r8,    4.53125000000000E-01_r8,    3.79865503229404E-01_r8,&
        3.1250000000000E-01_r8,    4.53125000000000E-01_r8,    3.77688516622590E-01_r8,&
        3.1250000000000E-01_r8,    4.68750000000000E-01_r8,    3.88797971748706E-01_r8,&
        3.2561322556469E-01_r8,    4.68750000000000E-01_r8,    3.90625000000000E-01_r8], [3,10]), &
        reshape([&
        3,    4,   9,  10,   5,&
        1,    2,   7,   8,   0,&
        1,    8,   9,   4,   0,&
        2,    3,   5,   6,   7,&
        1,    4,   3,   2,   0,&
        5,   10,   6,   0,   0,&
        10,   9,   8,   7,   6], [5,7]), reshape([&
        1,    2,&
        2,    3,&
        3,    4,&
        4,    1,&
        1,    8,&
        2,    7,&
        3,    5,&
        4,    9,&
        5,    6,&
        5,   10,&
        6,    7,&
        7,    8,&
        8,    9,&
        9,   10,&
        10,   6], [2,15]))
    call cutcube%split (P,tmp,ierr)
    tmpr1 = 0.0_r8
    tmpr2 = 0.0_r8
    success = ierr==0
    write(*,*) 'cutcube3 cut?                     ',success,tmp(1)%volume(),tmp(2)%volume()

    ! split the cutcube4
    P%normal = [-1.93230506727E-01_r8, 4.44887220431205E-01_r8, -8.74492614243836E-01_r8]
    P%rho    = -2.71937202667160926595E-01_r8
    call cutcube%init (ierr, reshape([&
        3.437500000000000000E-01_r8, 2.96875000000000E-01_r8, 3.906250000000000000E-01_r8,&
        3.593750000000000000E-01_r8, 2.96875000000000E-01_r8, 3.906250000000000000E-01_r8,&
        3.593750000000000000E-01_r8, 3.12500000000000E-01_r8, 3.906250000000000000E-01_r8,&
        3.593750000000000000E-01_r8, 3.12500000000000E-01_r8, 3.905375988910991802E-01_r8,&
        3.589794501821763073E-01_r8, 3.12500000000000E-01_r8, 3.906250000000000000E-01_r8,&
        3.437500000000000000E-01_r8, 3.05885355602952E-01_r8, 3.906250000000000000E-01_r8,&
        3.437500000000000000E-01_r8, 2.96875000000000E-01_r8, 3.860410971129954460E-01_r8,&
        3.593750000000000000E-01_r8, 2.96875000000000E-01_r8, 3.825885804765565278E-01_r8], [3,8]), &
        reshape([&
        5,   3,   4,   0,   0,&
        8,   2,   1,   7,   0,&
        7,   1,   6,   0,   0,&
        4,   3,   2,   8,   0,&
        6,   1,   2,   3,   5,&
        8,   7,   6,   5,   4], [5,6]), reshape([&
        1,   7,&
        2,   8,&
        3,   4,&
        1,   2,&
        2,   3,&
        3,   5,&
        1,   6,&
        4,   5,&
        5,   6,&
        6,   7,&
        7,   8,&
        8,   4], [2,12]))
    call cutcube%split (P,tmp,ierr)
    tmpr1 = tmp(1)%volume()
    tmpr2 = tmp(2)%volume()
    success = ierr==0
    write(*,*) 'cutcube4 cut?                     ',success,tmpr1,tmpr2

    ! ! split the cutcube5
    ! P%normal = [3.55921787180000004369E-01_r8, -9.34506859900000041996E-01_r8, 4.07556263139999958023E-03_r8]
    ! P%rho    = -2.42551881899006921417E-01_r8
    ! call cutcube%init (ierr, reshape([&
    !     2.96875E-01_r8,    3.75000E-01_r8,    5.31250E-01_r8,&
    !     3.12500E-01_r8,    3.75000E-01_r8,    5.31250E-01_r8,&
    !     3.12500E-01_r8,    3.90625E-01_r8,    5.31250E-01_r8,&
    !     2.96875E-01_r8,    3.90625E-01_r8,    5.31250E-01_r8,&
    !     2.96875E-01_r8,    3.75000E-01_r8,    5.46875E-01_r8,&
    !     3.12500E-01_r8,    3.75000E-01_r8,    5.46875E-01_r8,&
    !     3.12500E-01_r8,    3.90625E-01_r8,    5.46875E-01_r8,&
    !     2.96875E-01_r8,    3.90625E-01_r8,    5.46875E-01_r8], [3,8]), hex_f, hex_e)
    ! call cutcube%split(P,tmp,ierr)
    ! success = ierr==0
    ! write(*,*) 'cutcube5 cut?                     ',success,tmp(1)%volume(),tmp(2)%volume()
    ! call tmp(1)%print_data()

    ! split the cutcube5
    P%normal = [-3.559249373021E-01_r8, 9.345056485104E-01_r8, -4.07822371163224E-03_r8]
    P%rho    = 2.42549038109635933802E-01_r8
    call cutcube%init (ierr, reshape([&
        3.1250000000000000000000000E-01_r8,    3.7500000000000000000000000E-01_r8,    5.3125000000000000000000000E-01_r8,&
        2.9687500000000000000000000E-01_r8,    3.7500000000000000000000000E-01_r8,    5.4687500000000000000000000E-01_r8,&
        3.1250000000000000000000000E-01_r8,    3.7500000000000000000000000E-01_r8,    5.4687500000000000000000000E-01_r8,&
        2.9687500000000000000000000E-01_r8,    3.7500525766029885188501680E-01_r8,    5.4687500000000000000000000E-01_r8,&
        2.9687500000000000000000000E-01_r8,    3.7500000000000000000000000E-01_r8,    5.4566944384477711338377048E-01_r8,&
        2.9704011309117928085754556E-01_r8,    3.7500000000000000000000000E-01_r8,    5.3125000000000000000000000E-01_r8,&
        3.1250000000000000000000000E-01_r8,    3.8088814359134498532810653E-01_r8,    5.3125000000000000000000000E-01_r8,&
        3.1250000000000000000000000E-01_r8,    3.8095628719611474011230712E-01_r8,    5.4687500000000000000000000E-01_r8], [3,8]), &
        reshape([&
        6,   1,   3,   2,   5,&
        5,   2,   4,   0,   0,&
        8,   3,   1,   7,   0,&
        7,   1,   6,   0,   0,&
        4,   2,   3,   8,   0,&
        8,   7,   6,   5,   4], [5,6]), reshape([&
        1,   6,&
        1,   7,&
        2,   5,&
        1,   3,&
        2,   3,&
        3,   8,&
        2,   4,&
        4,   5,&
        5,   6,&
        6,   7,&
        7,   8,&
        8,   4], [2,12]))
    tmpr1 = cutcube%volume()
    call cutcube%split (P,tmp,ierr)
    tmpr1 = 0.0_r8
    tmpr2 = 0.0_r8
    success = ierr==0 .and. tmp(1)%volume() > 0.0_r8
    write(*,*) 'cutcube5 cut?                     ',success !,tmp(1)%volume(),tmp(2)%volume()

    write(*,*) '===================================================='
    write(*,*)

  end subroutine polyhedron_unit_test

end module polyhedron_type_test
