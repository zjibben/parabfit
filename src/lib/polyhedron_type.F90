!!
!! polyhedron_type
!!
!! This module defines an arbitrary polyhedron type, along with routines for
!! calculating volume, splitting polyhedra, locating intersections, etc.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2015
!!
!! References:
!!     1. Hopcroft and Kahn. A Paradigm for Robust Geometric Algorithms. Algorithmica, 1992
!! 

#include "f90_assert.fpp"

module polyhedron_type

  use kinds,  only: r8
  use consts, only: ndim
  use logging_services
  use polygon_type
  implicit none
  private

  public :: polyhedron_unit_test

  type, public :: polyhedron
    !private
    integer               :: nVerts, nEdges, nFaces       ! number of vertices, edges, and faces
    real(r8), allocatable :: x(:,:),face_normal(:,:)      ! vertex positions and face outward normals
    integer,  allocatable :: face_vid(:,:), edge_vid(:,:) ! face and edge IDs
    real(r8)              :: vol,x0(ndim),xl(ndim)        ! should be private, but need inheritance
    ! ideally vol should be private, but Fortran then also hides it from child types
  contains
    procedure, private :: init_polyhedron
    procedure, private :: init_polyhedron_copy
    procedure, private :: init_polyhedron_null
    generic            :: init => init_polyhedron, init_polyhedron_copy, init_polyhedron_null
    procedure          :: volume
    procedure          :: intersection_verts
    procedure          :: split
    procedure          :: volume_behind_plane
    procedure          :: print_data
    !procedure, private :: edge_containing_vertices
    procedure, private :: polyhedron_on_side_of_plane
    procedure, private :: update_face_normals
    procedure, private :: is_valid
    procedure, private :: remove_dangling_vertices
    procedure, private :: scale
    procedure, private :: descale
    !final :: polyhedron_delete
  end type polyhedron

contains
  
  ! subroutine polyhedron_delete (this)
  !   type(polyhedron) :: this
  !   if (allocated(this%x)) deallocate(this%x)
  !   if (allocated(this%face_normal)) deallocate(this%face_normal)
  !   if (allocated(this%face_vid)) deallocate(this%face_vid)
  !   if (allocated(this%edge_vid)) deallocate(this%edge_vid)
  ! end subroutine polyhedron_delete

  subroutine polyhedron_unit_test ()
    use plane_type
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

  subroutine init_polyhedron (this, ierr, x, face_v, edge_v, vol, face_normal)

    class(polyhedron),  intent(out) :: this
    integer,            intent(out) :: ierr
    real(r8),           intent(in)  :: x(:,:)
    integer,            intent(in)  :: face_v(:,:), edge_v(:,:)
    real(r8), optional, intent(in)  :: vol, face_normal(:,:)
    
    integer :: f,nV
    
    this%nVerts = size(x,     dim=2)
    this%nEdges = size(edge_v,dim=2)
    this%nFaces = size(face_v,dim=2)
    ierr = 0

    if (allocated(this%x))           deallocate(this%x)
    if (allocated(this%face_vid))    deallocate(this%face_vid)
    if (allocated(this%edge_vid))    deallocate(this%edge_vid)
    if (allocated(this%face_normal)) deallocate(this%face_normal)
    
    allocate(this%x(ndim,this%nVerts),& !this%face(this%nFaces),&
         this%face_vid(size(face_v,dim=1),this%nFaces),&
         this%edge_vid(size(edge_v,dim=1),this%nEdges),&
         this%face_normal(ndim,this%nFaces))
    
    this%x = x
    this%edge_vid = edge_v
    this%face_vid = face_v
    
    this%vol = merge(vol, 0.0_r8, present(vol))
    
    if (present(face_normal)) then
      this%face_normal = face_normal
    else
      this%face_normal = 0.0_r8
      call this%update_face_normals ()
    end if

    !call this%remove_dangling_vertices ()
    ierr = this%is_valid()
    if (ierr/=0) then
      write(*,*) 'tried to create invalid polyhedron!'
      call this%print_data()
      write(*,*)
    end if
    
    ! if the faces are of type polygon
    ! do f = 1,this%nFaces
    !   nV = count(face_v(:,f) /= 0) ! number of vertices on this face
    !   call this%face(f)%init (x(:,face_v(1:nV,f)))
    ! end do
    
    ! note that by taking the cross product between edges described in a
    ! counter-clockwise manner, we guarantee the normal to be outward facing
    
  end subroutine init_polyhedron

  function is_valid (this) result(ierr)

    class(polyhedron), intent(in) :: this
    integer :: ierr

    integer :: v

    ierr = 0

    if (this%nVerts < 4) ierr = 1 ! must have >= 4 vertices
    if (this%nFaces < 3) ierr = 1 ! must have >= 3 faces
    if (any(count(this%face_vid /= 0, dim=1) < 3)) ierr = 1 ! each face must have >= 3 vertices

    ! ! no dangling vertices
    ! ! every vertex must be connected to at least three other vertices
    ! do v = 1,this%nVerts
    !   if (count(this%edge_vid==v)<3) then
    !     write(*,*) 'dangling vertex ',v
    !     ierr = 1
    !   end if
    ! end do

  end function is_valid

  subroutine update_face_normals (this, force)

    use array_utils,   only: normalize, isZero
    use cell_geometry, only: cross_product

    class(polyhedron), intent(inout) :: this
    logical, intent(in), optional :: force

    integer :: f,v,nV
    logical :: forceh

    forceh = merge(force, .false., present(force))
    if (forceh) this%face_normal = 0.0_r8 ! force the normals to be recalculated

    do f = 1,this%nFaces
      if (all(isZero(this%face_normal(:,f)))) then
        nV = count(this%face_vid(:,f)/=0)
        v = 3
        do while (all(isZero (this%face_normal(:,f))) .and. v<=nV)
          this%face_normal(:,f) = normalize(cross_product (&
              this%x(:,this%face_vid(2,f)) - this%x(:,this%face_vid(1,f)), &
              this%x(:,this%face_vid(v,f)) - this%x(:,this%face_vid(1,f))))
          v = v + 1
        end do
        if (v>nV .and. all(isZero(this%face_normal(:,f)))) then
          call this%print_data ()
          call LS_fatal ("could not calculate polyhedron face normal")
        end if
      end if
    end do

  end subroutine update_face_normals

  subroutine remove_dangling_vertices (this)

    class(polyhedron), intent(inout) :: this

    integer :: v,vv

    ! no dangling vertices
    ! every vertex must be connected to at least three other vertices
    ! vertices connected to only two other vertices may be removed, and
    ! its connected vertices joined
    ! WARNING: this assumes the dangling vertex is still attached to two other vertices,
    !          and doesn't consider vertices with only one edge
    v = 1
    do while (v <= this%nVerts)
      if (count(this%edge_vid==v)<3) then
        ! delete this vertex
        do vv = v+1,this%nVerts
          this%x(:,vv-1) = this%x(:,vv)
        end do
        this%nVerts = this%nVerts - 1

        ! delete the two edges it is connected to, and add a new one between those two vertices
        call delete_edges_containing_vertex (this,v)

        ! remove the vertex from the face it is connected to
        call delete_faces_containing_vertex (this,v)
      else
        v = v+1
      end if
    end do

  contains

    subroutine delete_edges_containing_vertex (this,v)

      class(polyhedron), intent(inout) :: this
      integer, intent(in) :: v

      integer :: e,ee,vn(2),j

      e = 1; j=1
      do while (e <= this%nEdges)
        if (any(this%edge_vid(:,e)==v)) then
          ! store the other vertex connected to this edge
          vn(j) = this%edge_vid(1,e)
          if (vn(j)==v) vn(j) = this%edge_vid(2,e)
          j = j+1

          ! delete edges containing this vertex
          do ee = e+1,this%nEdges
            this%edge_vid(:,ee-1) = this%edge_vid(:,ee)
          end do
          this%nEdges = this%nEdges - 1
        else
          e = e+1
        end if
      end do

      ! add the vertex connecting [v1,v2]
      this%nEdges = this%nEdges + 1
      this%edge_vid(:,this%nEdges) = vn

      ! update the edge structure
      where (this%edge_vid>v) this%edge_vid = this%edge_vid - 1

    end subroutine delete_edges_containing_vertex

    subroutine delete_faces_containing_vertex (this,v)

      use array_utils, only: first_true_loc

      class(polyhedron), intent(inout) :: this
      integer, intent(in) :: v

      integer :: f,ff,fv,nV

      f = 1
      do while (f <= this%nFaces)
        if (any(this%face_vid(:,f)==v)) then
          if (count(this%face_vid(:,f)/=0)==3) then
            ! if the face only has three vertices, delete the entire face
            do ff = f+1,this%nFaces
              this%face_vid(:,ff-1) = this%face_vid(:,ff)
              this%face_normal(:,ff-1) = this%face_normal(:,ff)
            end do
            this%nFaces = this%nFaces - 1
          else
            nV = count(this%face_vid(:,f)/=0)
            fv = first_true_loc (this%face_vid(:,f)==v)
            do ff = fv+1,nV
              this%face_vid(ff-1,f) = this%face_vid(ff,f)
            end do
            this%face_vid(nV,f) = 0
            f = f+1
          end if
        else
          f = f+1
        end if
      end do

      ! update the face structure
      where (this%face_vid>v) this%face_vid = this%face_vid - 1

    end subroutine delete_faces_containing_vertex

  end subroutine remove_dangling_vertices

  subroutine init_polyhedron_null (this)
    class(polyhedron), intent(out) :: this
    this%nVerts = 0
    this%nEdges = 0
    this%nFaces = 0
    this%vol = 0.0_r8
  end subroutine init_polyhedron_null

  subroutine init_polyhedron_copy (this,poly)
    class(polyhedron), intent(out) :: this
    class(polyhedron), intent(in)  :: poly

    this%nVerts = poly%nVerts
    this%nEdges = poly%nEdges
    this%nFaces = poly%nFaces
    this%vol    = poly%vol

    if (allocated(this%x))           deallocate(this%x)
    if (allocated(this%edge_vid))    deallocate(this%edge_vid)
    if (allocated(this%face_vid))    deallocate(this%face_vid)
    if (allocated(this%face_normal)) deallocate(this%face_normal)

    this%x = poly%x
    this%edge_vid = poly%edge_vid
    this%face_vid = poly%face_vid
    this%face_normal = poly%face_normal

  end subroutine init_polyhedron_copy

  subroutine scale (this)

    class(polyhedron), intent(inout) :: this

    integer :: v

    this%x0 = minval(this%x,dim=2)
    this%xl = maxval(this%x,dim=2) - this%x0
    do v = 1,this%nVerts
      this%x(:,v) = (this%x(:,v)-this%x0)/this%xl
    end do
    call this%update_face_normals (force=.true.)

  end subroutine scale

  subroutine descale (this)

    class(polyhedron), intent(inout) :: this

    integer :: v

    do v = 1,this%nVerts
      this%x(:,v) = this%x(:,v)*this%xl + this%x0
    end do
    call this%update_face_normals (force=.true.)
    this%vol = this%vol * product(this%xl)

  end subroutine descale

  ! This function calculates the volume of a polehedron following the
  ! algorithm presented by [1]. It calculates the sum of the surface
  ! integral of x over all faces, which is equal to the volume by the
  ! divergence theorem.
  !
  ! note 1: this is done so we can find the volume of very tiny ones without hitting precision limits
  !         in some cases, very tiny polyhedra (without this trick) would produce negative volumes
  !         because their volume was on the order of floating point errors. Scaling the polyhedron
  !         up to some normalized size before calculating the volume, then scaling back, seems to
  !         counter this rather well, but seems a bit tricky. It might instead be worth just setting
  !         the volume of a polyhedron to zero if it is calculated as below zero and very close to
  !         zero. Note this occurs even though the polyhedra splitting routine is "robust" and will
  !         not allow vertices to be within some distance alpha (=1e-9). These tiny distances can
  !         still produce volumes on the order of e-27, far below the double precision limit of e-16.
  real(r8) function volume (this)

    use consts,          only: alittle
    use array_utils,     only: isZero, normalize
    use ieee_arithmetic, only: ieee_is_nan
    use cell_geometry,   only: cross_product

    class(polyhedron), intent(inout) :: this
    !integer, intent(out) :: ierr
    
    integer :: f,nV,v
    real(r8) :: tmp(ndim), n(ndim), t(ndim)
    type(polygon) :: face

    !ierr = 0
    volume = merge(this%vol, 0.0_r8, .not.ieee_is_nan(this%vol))
    if (.not.allocated(this%x) .or. volume > 0.0_r8) return

    ! scale the polyhedron (see note 1)
    call this%scale ()
    
    
    ! sum up the integral of n_x*x over all faces (could just as easily be any other direction)
    volume = 0.0_r8

    ! do f = 1,this%nFaces
    !   ! generate a polygon from this face
    !   nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
    !   call face%init (this%x(:,this%face_vid(1:nV,f)), this%face_normal(:,f))

    !   ! calculate this face's contribution
    !   if (.not.isZero (face%norm(1))) volume = volume + face%norm(1) * face%intXdA (1)
    ! end do
    
    !print *, 'calculating volume...'
    do f = 1,this%nFaces
      nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
      tmp = 0.0_r8
      do v = 1,nV
        tmp = tmp + &
            cross_product (this%x(:,this%face_vid(v,f)), this%x(:,this%face_vid(modulo(v,nV)+1,f)))
      end do
      volume = volume + &
          dot_product(this%face_normal(:,f),this%x(:,this%face_vid(1,f))) * &
          dot_product(this%face_normal(:,f),tmp)
      ! print '(a,4es14.4)', 'volume: ',volume
      ! print *, this%face_vid(:,f)
      ! print *, tmp
      ! print *, this%face_normal(:,f)
      ! print *, this%x(:,this%face_vid(1,f))
      ! do v = 1,nV
      !   print *, 'x: ',this%x(:,this%face_vid(1,v))
      ! end do
    end do
    volume = volume/6.0_r8
    this%vol = volume

    ! do f = 1,this%nFaces
    !   nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
    !   do v = 1,nV
    !     t = normalize(this%x(:,this%face_vid(modulo(v,nV)+1,f))-this%x(:,this%face_vid(v,f)))
    !     n = cross_product(this%face_normal(:,f),t)
    !     volume = volume + &
    !         dot_product(this%x(:,this%face_vid(v,f)),this%face_normal(:,f)) * &
    !         dot_product(this%x(:,this%face_vid(v,f)),t) * &
    !         dot_product(this%x(:,this%face_vid(v,f)),n)

    !     t = normalize(this%x(:,this%face_vid(modulo(v-2,nV)+1,f))-this%x(:,this%face_vid(v,f)))
    !     n = -cross_product(this%face_normal(:,f),t)
    !     volume = volume + &
    !         dot_product(this%x(:,this%face_vid(v,f)),this%face_normal(:,f)) * &
    !         dot_product(this%x(:,this%face_vid(v,f)),t) * &
    !         dot_product(this%x(:,this%face_vid(v,f)),n)
    !   end do
    ! end do
    ! volume = volume/6.0_r8

    ! rescale the polyhedron
    call this%descale ()
    volume = this%vol

    if (this%vol < 0.0_r8) then
      if (isZero(this%vol, 1e5_r8*alittle)) then
        ! if this polyhedron has a volume of almost zero, make it zero
        ! this seems to be necessary for very tiny polyhedrons,
        ! where floating point errors may make the volume calculation drop below zero
        volume = 0.0_r8
        this%vol = 0.0_r8
        deallocate(this%x)
        this%nVerts = 0
      else
        write(*,*) "calculated negative polyhedron volume!"
        call this%print_data ()
        write(*,*)
        !ierr = 1
        !call LS_fatal ("calculated negative polyhedron volume!")
      end if
    end if

  end function volume

  !
  ! Given an equation of a plane and a polyhedron, return a polygon from the
  ! points where the plane intersects the polyhedron edges
  !
  type(polygon) function intersection_verts (this,P,v_assoc_pe)

    use plane_type
    use array_utils, only: reverse, containsPoint, pointIndex

    class(polyhedron), intent(in)  :: this
    class(plane),      intent(in)  :: P
    integer, optional, intent(out) :: v_assoc_pe(:)

    integer  :: e,Nintersections
    real(r8) :: x(ndim,this%nEdges),intx(ndim)

    if (present(v_assoc_pe)) v_assoc_pe = -1
    
    ! loop through all edges
    Nintersections = 0
    do e = 1,this%nEdges
      ! check if the P intersects this edge
      if (P%intersects(this%x(:,this%edge_vid(:,e)))) then
        ! if it does, find the point where they intersect
        intx = P%intersection_point (this%x(:,this%edge_vid(:,e)))

        if ((all(intx==this%x(:,this%edge_vid(1,e))) .or. all(intx==this%x(:,this%edge_vid(2,e)))) &
            .and. Nintersections>0 .and. containsPoint(intx, x(:,1:Nintersections))) then
          ! if the point was already found, note that this edge intersects with that found point
          ! this particularly comes into effect when we intersect a vertex
          if (present(v_assoc_pe)) v_assoc_pe(e) = pointIndex(intx, x(:,1:Nintersections))
        else
          Nintersections = Nintersections + 1
          x(:,Nintersections) = intx
          if (present(v_assoc_pe)) v_assoc_pe(e) = Nintersections
        end if
      end if
    end do
    
    ! pass the intersection points to the polygon constructor
    if (Nintersections>2) then
      call intersection_verts%init (x(:,1:Nintersections))
      
      ! this probably doesn't need to be called every time this function is used
      if (present(v_assoc_pe)) then
        call intersection_verts%order (v_assoc_pe)
      else
        call intersection_verts%order ()
      end if
      
      ! make sure the vertices are ordered counter-clockwise
      if (sum(intersection_verts%norm*P%normal)<0.0_r8) then
        ! reverse both x and v_assoc_pe
        intersection_verts%x = reverse (intersection_verts%x)
        
        if (present(v_assoc_pe)) then
          do e = 1,size(v_assoc_pe)
            if (v_assoc_pe(e)>0) v_assoc_pe(e) = Nintersections - v_assoc_pe(e)+1
          end do
        end if
        
        call intersection_verts%update_plane_normal ()
      end if
    ! else
    !   print *, Nintersections, " intersection points"
    !   call this%print_data()
    !   call P%print_data()
    !   call LS_fatal ("not enough intersection points")
    end if
    
  end function intersection_verts

  ! return a list of the edge ids for edges intersected by the plane
  function intersected_edges (this,P) result(inte)
    use plane_type

    class(polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P

    integer :: e
    logical :: inte(this%nEdges)
    
    do e = 1,this%nEdges
      inte(e) = P%intersects(this%x(:,this%edge_vid(:,e)))
    end do

  end function intersected_edges

  ! return two polyhedrons produced by dividing a given polyhedron with a plane
  ! the first element returned is in front of the plane
  ! the second element returned is behind the plane
  !
  ! note 1: this part is not part of the Hopcroft and Kahn algorithm
  !         When splitting a very thin polyhedron, the intersection
  !         vertices may end up within alpha of each other. In this
  !         case, say both polyhedra are null, since we are within
  !         the cutvof anyways.
  subroutine split (this,P,split_poly,ierr)

    use consts, only: alpha
    use plane_type
    
    class(polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P
    type(polyhedron), intent(out) :: split_poly(:)
    integer,          intent(out) :: ierr

    integer       :: v, v_assoc_pe(this%nEdges),side(this%nVerts)
    type(polygon) :: intpoly
    real(r8)      :: dist, tmp1, tmp2

    ASSERT(size(split_poly)==2)
    ierr = 0
    
    ! check which side of the plane each vertex lies
    ! vertices within distance alpha of the plane are considered to lie on it
    !  dist  >  alpha => side =  1
    !  dist  < -alpha => side = -1
    ! |dist| <  alpha => side =  0
    do v = 1,this%nVerts
      dist = P%signed_distance (this%x(:,v)) ! calculate the signed distance from the plane
      side(v) = merge(int(sign(1.0_r8, dist)), 0, abs(dist) > alpha) ! decide where it lies
    end do
    
    if (.not.any(side<0)) then
      split_poly(1) = this
      call split_poly(2)%init ()
    else if (.not.any(side>0)) then
      split_poly(2) = this
      call split_poly(1)%init ()
    else
      intpoly = this%intersection_verts (P,v_assoc_pe)

      if (intpoly%nVerts > 2) then
        split_poly(1) = this%polyhedron_on_side_of_plane (P,  1, side, intpoly, v_assoc_pe, ierr)
        if (ierr/=0) call LS_fatal ("polyhedron split failed: invalid child")
        split_poly(2) = this%polyhedron_on_side_of_plane (P, -1, side, intpoly, v_assoc_pe, ierr)
        if (ierr/=0) call LS_fatal ("polyhedron split failed: invalid child")
      else ! see note 1
        call split_poly(1)%init ()
        call split_poly(2)%init ()
      end if
    end if

    ! if any of the polyhedrons have a face with less than 3 vertices, throw an error
    do v = 1,2
      if (allocated(split_poly(v)%face_vid)) then
        if (any(count(split_poly(v)%face_vid /= 0, dim=1) < 3)) then
          call split_poly(v)%print_data (normalized=.true.)
          write(*,*)
          ierr = 1
        end if
      end if
    end do
    if (ierr/=0) call LS_fatal ("polyhedron split failed--one of the children has an invalid face")

    ! if either of the volumes are less than zero, throw an error
    ! TODO: this is a more expensive check. might be worth skipping when sufficiently confident.
    tmp1 = split_poly(1)%volume()
    tmp2 = split_poly(2)%volume()
    if (tmp1 < 0.0_r8 .or. tmp2 < 0.0_r8) then
      write(*,*)
      write(*,*) 'parent:'
      call this%print_data ()
      write(*,*)

      write(*,*) 'child1:'
      call split_poly(1)%print_data ()
      write(*,*)

      write(*,*) 'child2:'
      call split_poly(2)%print_data ()
      write(*,*)

      write(*,*) 'other data:'
      write(*,'(15i3)') side
      call P%print_data()

      write(*,*) 
      call intpoly%print_data ()
      write(*,*)
      
      write(*,*) 'problematic vols: ',tmp1,tmp2
      ierr = 1
      !call LS_fatal ('polyhedron split failed: invalid volume')
    end if

  end subroutine split

  ! return the volume behind (opposite normal) a plane and inside the polyhedron
  real(r8) function volume_behind_plane (this,P,ierr)

    use consts,          only: alpha
    use ieee_arithmetic, only: ieee_is_nan
    use plane_type

    class(polyhedron), intent(in) :: this
    class(plane),      intent(in) :: P
    integer,           intent(out) :: ierr

    real(r8)         :: dist
    integer          :: v, side(this%nVerts), v_assoc_pe(this%nEdges)
    type(polyhedron) :: behind, split(2)
    type(polygon)    :: intpoly

    ierr = 0

    ! check which side of the plane each vertex lies
    ! vertices within distance alpha of the plane are considered to lie on it
    !  dist  >  alpha -> side =  1
    !  dist  < -alpha -> side = -1
    ! |dist| <  alpha -> side =  0
    do v = 1,this%nVerts
      dist = P%signed_distance (this%x(:,v)) ! calculate the signed distance from the plane
      side(v) = merge(int(sign(1.0_r8, dist)), 0, abs(dist) > alpha) ! decide where it lies
    end do
    
    if (.not.any(side>0)) then
      behind = this
    else if (any(side>0) .and. any(side<0)) then
      intpoly = this%intersection_verts (P,v_assoc_pe)
      if (intpoly%nVerts > 2) then
        behind = this%polyhedron_on_side_of_plane (P, -1, side, intpoly, v_assoc_pe,ierr)
        if (ierr /= 0) then
          call dumpData ()
          write(*,*) 'volume_behind_plane: polyhedron split failed'
          volume_behind_plane = 0.0_r8
          return
          call LS_fatal ("polyhedron split failed")
        end if
      else
        volume_behind_plane = 0.0_r8
        return
      end if
    else
      volume_behind_plane = 0.0_r8
      return
    end if
    
    ! if any of the polyhedrons have a face with less than 3 vertices, throw an error
    if (allocated(behind%face_vid)) then
      if (any(count(behind%face_vid /= 0, dim=1) < 3)) then
        call dumpData ()
        call LS_fatal ("polyhedron split failed: invalid face")
      end if
    end if

    ! calculate the volume of the polyhedron behind the plane
    volume_behind_plane = behind%volume ()

    if (ieee_is_nan(volume_behind_plane) .or. volume_behind_plane < 0.0_r8) then
      call dumpData ()
      write(*,*) 'problematic vol: ',volume_behind_plane
      call LS_fatal ('polyhedron split failed: invalid volume')
    end if
    
  contains

    subroutine dumpData ()
      write(*,*)
      write(*,*) 'parent:'
      call this%print_data ()
      write(*,*)
      
      write(*,*) 'child:'
      call behind%print_data ()
      write(*,*)
      
      write(*,*) 'other data:'
      write(*,'(15i3)') side
      call P%print_data()
      
      write(*,*) 
      call intpoly%print_data ()
      write(*,*)
    end subroutine dumpData

  end function volume_behind_plane

  ! WARNING: need to update this routine to pass face normals down to the child polyhedron
  ! Reference [1]
  ! note 1: In this case, the polyhedron may not have been clearly all within the half-space.
  !         Some vertices on the face may have landed on each side of the plane, but if
  !         more than three landed on the plane itself, we consider this equivalent to
  !         the entire face landing on the plane.
  type(polyhedron) function polyhedron_on_side_of_plane (this,P,valid_side,side,intpoly,v_assoc_pe,ierr)

    use plane_type

    class(polyhedron), intent(in) :: this
    type(plane),       intent(in) :: P
    integer,           intent(in) :: side(:)       ! gives which side of the plane vertices lie on
    integer,           intent(in) :: valid_side    ! +/- 1, indicating the side of the plane we want
    type(polygon),     intent(in) :: intpoly       ! polygon of intersection coordinates
    integer,           intent(in) :: v_assoc_pe(:) ! intersection polygon vertex id for a given parent edge id
    integer,           intent(out) :: ierr

    integer :: pcf, nVerts, nParVerts, nEdges, nFaces, tmp, invalid_side, &
         p2c_vid(this%nVerts), & ! parent to child vertex id conversion table for cases they correspond
         edge_vid(2,this%nEdges+intpoly%nVerts), & ! can have intpoly%nVerts more edges than parent
         face_vid(size(this%face_vid,dim=1)+2,this%nFaces+1) ! could have 1 more face and faces could be 2 longer than parent
    
    ! note: an updated planar face can only include 1 more node than the original,
    !       but I'm not sure if there is a limit to how many nodes the entirely new face can have.
    !       For cubes the number is 2.
    real(r8)      :: x(3,this%nVerts+intpoly%nVerts)

    call polyhedron_on_side_of_plane%init ()

    invalid_side = -valid_side

    if (.not.any(side==-valid_side)) then ! the entire polyhedron is within the half-space
      polyhedron_on_side_of_plane = this
    else if (any(side==-valid_side) .and. any(side==valid_side)) then ! the polyhedron is split
      pcf = plane_coinciding_face(this,side,valid_side)
      if (pcf > 0) then
        ! the polyhedron really landed entirely inside or entirely outside the half-space (see note 1)
        if (dot_product(this%face_normal(:,pcf), P%normal) > 0.0_r8) &
            polyhedron_on_side_of_plane = this
      else
        ! make a list of vertices for the new polyhedron
        call generate_new_verts (x,p2c_vid,nVerts,nParVerts, this,side,v_assoc_pe,valid_side)

        ! find new edges, knowing the original structure and the interface polygon
        call find_edges (edge_vid,nEdges, this,side,valid_side,p2c_vid,nParVerts,nVerts,v_assoc_pe)

        ! construct a set of faces from the edge information
        ! note these will not be in a particular order, like pececillo would expect for hexes
        call find_faces (face_vid,nFaces,ierr, &
            this,side,valid_side,p2c_vid,nParVerts,nVerts,intpoly%nVerts,v_assoc_pe)
        if (ierr /= 0) call fatal()
        
        ! initialize final polyhedron
        tmp = maxval(count(face_vid(:,:) /= 0,dim=1)) ! the maximum number of vertices on a face
        if (nVerts < 3) then
          call this%print_data()
          write(*,*)
          call P%print_data()
          call LS_fatal ("not enough vertices for a polygon!")
        end if
        call polyhedron_on_side_of_plane%init (ierr, x(:,1:nVerts), face_vid(1:tmp,1:nFaces), &
            edge_vid(:,1:nEdges))
        
        !write(*,*) 'poly', nVerts, tmp, nFaces, nEdges
        if (ierr /= 0) call fatal()
      end if
    end if

  contains

    subroutine fatal ()
      print *, 'parent: '
      call this%print_data()
      write(*,*)

      write(*,*) 'other data:'
      write(*,'(15i3)') side
      call P%print_data()

      write(*,*)
      call intpoly%print_data ()
      write(*,*)

      call LS_fatal ("polyhedron_on_side_of_plane failed: invalid child polyhedron")
    end subroutine fatal

    ! finds the polyhedron face which coincides with the face, if one exists
    integer function plane_coinciding_face (this,side,valid_side)

      class(polyhedron), intent(in) :: this
      integer, intent(in) :: side(:), valid_side

      integer :: nV

      do plane_coinciding_face = 1,this%nFaces
        nV = count(this%face_vid(:,plane_coinciding_face) > 0)
        ! the face coincides with the plane if at least 3 of its vertices lie on the plane
        if (count(side(this%face_vid(1:nV,plane_coinciding_face))==0)>=3) return
      end do
      plane_coinciding_face = -1

    end function plane_coinciding_face

    ! note 1: these are added entirely at the end, rather than including the parent vertices
    !         that lie on the plane in the above loop because having all these vertices together
    !         makes it easy for constructing the new face later.
    subroutine generate_new_verts (x,p2c_vid,nVerts,nParVerts, this,side,v_assoc_pe,valid_side)

      use array_utils, only: magnitude,isZero

      real(r8),         intent(out) :: x(:,:)
      integer,          intent(out) :: p2c_vid(:), nVerts, nParVerts
      type(polyhedron), intent(in)  :: this
      integer,          intent(in)  :: side(:), v_assoc_pe(:), valid_side

      integer :: v,iv

      nVerts = 0; p2c_vid = 0

      ! first get the vertices fully on the valid side of the intersection plane
      do v = 1,this%nVerts
        if (side(v)==valid_side) then
          nVerts = nVerts+1
          x(:,nVerts) = this%x(:,v)
          p2c_vid(v) = nVerts ! parent to child vertex id
        end if
      end do
      nParVerts = nVerts ! number of vertices acquired from parent

      ! then get the vertices from the plane-polyhedron intersection (see note 1)
      x(:,nParVerts+1:nParVerts+intpoly%nVerts) = intpoly%x
      nVerts = nVerts + intpoly%nVerts

      if (nVerts < 3) call LS_fatal("not enough vertices to make a polyhedron")

      ! update the parent to child vertex id table with vertices that lie on the plane
      do v = 1,this%nVerts
        if (side(v)==0) then
          ! search for the interface polygon vertex which coincides with this vertex
          do iv = 1,intpoly%nVerts
            !if (isZero(magnitude(this%x(:,v) - intpoly%x(:,iv)))) then
            if (all(isZero(this%x(:,v) - intpoly%x(:,iv)))) then
              p2c_vid(v) = nParVerts + iv
              exit
            end if
          end do
          if (p2c_vid(v)==0) call LS_fatal ("parent point not found in intersecting polygon")
        end if
      end do

    end subroutine generate_new_verts

    subroutine find_edges (edge_vid,nEdges, this,side,valid_side,p2c_vid,nParVerts,nVerts,v_assoc_pe)

      use array_utils, only: containsPair

      integer,          intent(out) :: edge_vid(:,:), nEdges
      type(polyhedron), intent(in)  :: this
      integer,          intent(in)  :: side(:), valid_side, p2c_vid(:), nParVerts, nVerts, &
          v_assoc_pe(:)

      integer :: e, v, n_on_valid_side, new_edge(2)

      nEdges = 0; edge_vid = 0

      ! add edges that were a part of the parent
      do e = 1,this%nEdges
        ! how many of the edge vertices are on the valid side of the plane?
        n_on_valid_side = count(side(this%edge_vid(:,e))==valid_side)

        if (n_on_valid_side>0) then
          if (n_on_valid_side==2) then ! this entire edge is on the valid side of the plane
            new_edge = p2c_vid(this%edge_vid(:,e))
          else if (n_on_valid_side==1) then ! this edge was divided (or one of the vertices lies on the plane)
            ! the edge consists of the vertex on this side of the plane and the intersection point
            if (side(this%edge_vid(1,e))==valid_side) then
              new_edge = [p2c_vid(this%edge_vid(1,e)),nParVerts+v_assoc_pe(e)]
            else
              new_edge = [p2c_vid(this%edge_vid(2,e)),nParVerts+v_assoc_pe(e)]
            end if
          end if

          ! if this edge isn't already listed, add it
          ! the edge might already be listed in cases where we have very close vertices (O(alpha)),
          ! which are then combined in the new polyhedron.
          if (.not.containsPair(new_edge, edge_vid(:,1:nEdges))) then
            nEdges = nEdges + 1
            edge_vid(:,nEdges) = new_edge
          end if
        end if
      end do

      ! points on the plane make up edges with each other
      do v = nParVerts+1,nVerts-1
        nEdges = nEdges + 1
        edge_vid(:,nEdges) = [v,v+1]
      end do
      nEdges = nEdges + 1
      edge_vid(:,nEdges) = [nVerts,nParVerts+1] ! complete the loop

    end subroutine find_edges

    ! note 1: There may not be a valid edge between the previously found vertex and this one.
    !         This particularly happens when there are multiple vertices on this face which
    !         also lie on the plane. In that case, we loop through those points, adding them
    !         until we find an edge between one and the next vertex.
    subroutine find_faces (face_vid,nFaces,ierr, this,side,valid_side,p2c_vid,nParVerts,nVerts,&
        nPolyVerts,v_assoc_pe)

      use array_utils, only: containsPair,xrange

      integer,          intent(out) :: face_vid(:,:),nFaces,ierr
      type(polyhedron), intent(in)  :: this
      integer,          intent(in)  :: side(:), valid_side, p2c_vid(:), nParVerts, nVerts, &
          nPolyVerts, v_assoc_pe(:)

      integer :: v,nV,f, e, edge_cont_verts(this%nVerts,this%nVerts)

      ! first make a lookup table for finding the edge associated with a pair of vertices
      edge_cont_verts = 0
      do e = 1,this%nEdges
        edge_cont_verts(this%edge_vid(1,e),this%edge_vid(2,e)) = e
        edge_cont_verts(this%edge_vid(2,e),this%edge_vid(1,e)) = e
      end do

      nFaces = 0; face_vid = 0; ierr = 0

      ! loop over all of the parent's faces, adding, modifying, or ignoring it's faces as needed
      do f = 1,this%nFaces
        ! if any vertices for this original face are on the valid side of the plane,
        ! then the face structure can be preserved or modified.
        ! Otherwise, it is thrown out since the entire face is behind the plane.
        nV = count(this%face_vid(:,f) /= 0) ! number of vertices on this face
        if (any(side(this%face_vid(1:nV,f))==valid_side)) then
          nFaces = nFaces + 1
          if (all(side(this%face_vid(1:nV,f))/=-valid_side)) then
            ! if all vertices are on the valid side of the face, this face is preserved exactly
            face_vid(1:nV,nFaces) = p2c_vid(this%face_vid(1:nV,f))
          else
            ! if some vertices are on the valid side of the face, this face is modified
            ! write(*,*) 
            ! write(*,*) 'f: ',nFaces
            ! write(*,*) this%face_vid(1:nV,f)
            ! write(*,*) side(this%face_vid(1:nV,f))
            ! write(*,'(20i3)') v_assoc_pe
            call modified_parent_face (face_vid(:,nFaces),ierr, this%face_vid(1:nV,f), p2c_vid, &
                edge_vid(:,1:nEdges), side, valid_side, nParVerts, nPolyVerts, v_assoc_pe, &
                edge_cont_verts)
            
            ! check that the face is valid by making sure each pair of points corresponds to an edge
            nV = count(face_vid(:,nFaces)/=0)
            do v = 1,nV
              if (.not.containsPair(face_vid([v,mod(v,nV)+1],nFaces),edge_vid(:,1:nEdges))) then
                write(*,*) "find_faces: face contains edge that doesn't exist"
                write(*,*) 'vertices: ',face_vid([v,mod(v,nV)+1],nFaces)
                ierr = 1
                exit
              end if
            end do
            
            if (ierr /= 0) then
              nV = count(this%face_vid(:,f) /= 0)
              write(*,*) 'invalid face'
              write(*,'(a,12i3)') 'face: ',face_vid(:,nFaces)
              write(*,'(a,  i3)') 'valid_side: ',valid_side
              write(*,'(a,15i3)') 'parent vertex sides ',side
              write(*,'(a,15i3)') 'parent face vertices ',this%face_vid(1:nV,f)
              write(*,'(a,15i3)') 'parent face vertex sides ',side(this%face_vid(1:nV,f))
              write(*,'(a,15i3)') 'p2c_vid ',p2c_vid(this%face_vid(1:nV,f))
              write(*,'(a,15i3)') 'p2c_vid ',p2c_vid
              do e = 1,nEdges
                write(*,'(a,i4,a,2i4)') 'edge ',e,': ',edge_vid(:,e)
              end do
              do e = 1,nVerts
                write(*,'(a,i3,a,3es35.25)') 'x ',e,':  ',x(:,e)
              end do
              write(*,*)
              return
            end if

            ! ! if this face is invalid, delete it.
            ! ! this line is here to counter a common, but easily fixable problem:
            ! ! invalid (and unnecessary) faces generated when two edges are very
            ! ! close together (<2alpha).
            ! ! WARNING: this is also a prime target for subtle errors, and should
            ! !          really be figured out and handled in a more elegant manner
            ! if (count(face_vid(:,nFaces)>0) < 3) nFaces = nFaces - 1
          end if
        end if
      end do

      ! add the new face from points lying on the plane
      nFaces = nFaces + 1
      if (valid_side>0) then ! needs to be opposite the input order
        face_vid(1:nVerts-nParVerts,nFaces) = xrange (nVerts,nParVerts+1)
      else
        face_vid(1:nVerts-nParVerts,nFaces) = xrange (nParVerts+1,nVerts)
      end if

    end subroutine find_faces

    ! calculate the face by modifying the parent face structure with the new edge from the plane
    !
    ! note 1: in cases where several vertices on the parent face lie on the plane,
    !         this vertex may not be connected to the previous one found, but instead
    !         one of the other previous ones. If this edge does not exist and the
    !         last point was also on the plane, delete the last found point
    subroutine modified_parent_face (face_vid, ierr, par_face_vid, p2c_vid, edge_vid, side, &
        valid_side, nParVerts, nPolyVerts, v_assoc_pe, edge_cont_verts)

      use array_utils, only: first_true_loc, last_true_loc, containsPair

      integer, intent(out) :: face_vid(:), ierr
      integer, intent(in) :: par_face_vid(:), p2c_vid(:), edge_vid(:,:), side(:), valid_side, &
          nParVerts, nPolyVerts, v_assoc_pe(:), edge_cont_verts(:,:)

      integer :: v, cvid, new_v, ecv, pv, nV

      face_vid = 0; ierr = 0
      nV = size(par_face_vid)
      cvid = 1

      ! start with the first valid vertex after the section intersected by the plane
      v = first_true_loc (side(par_face_vid)==valid_side)
      if (v==1) v = modulo(last_true_loc (side(par_face_vid)/=valid_side),nV)+1
      !write(*,*) 'v', v, nV, modulo(v-2,nV)+1, modulo(v-1,nV)+1
      
      ! the edge [v-1,v] *must* intersect with the plane
      ecv = edge_cont_verts(par_face_vid(modulo(v-2,nV)+1),par_face_vid(modulo(v-1,nV)+1))
      if (ecv < 1) then
        write(*,*) 'find faces failed: edge missing'
        write(*,*) 'vertices: ',par_face_vid(modulo(v-2,nV)+1),par_face_vid(modulo(v-1,nV)+1)
        ierr = 1
        return
      end if
      face_vid(cvid) = nParVerts + v_assoc_pe(ecv)
      cvid = cvid+1
      !write(*,*) '1',face_vid(cvid-1), ecv, nParVerts, v_assoc_pe(ecv)

      ! add all the parent's valid vertices in sequence, stop when we hit one that is invalid
      ! this stopping point indicates the next edge [v-1,v] which *must* be intersected by the plane
      do while (side(par_face_vid(modulo(v-1,nV)+1))==valid_side)
        face_vid(cvid) = p2c_vid(par_face_vid(modulo(v-1,nV)+1))
        v = v+1; cvid = cvid+1
        !write(*,*) '2',face_vid(cvid-1)
      end do

      ! the second new vertex is on this edge [v-1,v]
      ecv = edge_cont_verts(par_face_vid(modulo(v-2,nV)+1),par_face_vid(modulo(v-1,nV)+1))
      if (ecv < 1) then
        write(*,*) 'find faces failed: edge missing'
        write(*,*) 'vertices: ',par_face_vid(modulo(v-2,nV)+1),par_face_vid(modulo(v-1,nV)+1)
        ierr = 1
        return
      end if
      face_vid(cvid) = nParVerts + v_assoc_pe(ecv)
      cvid = cvid+1
      !write(*,*) '3',face_vid(cvid-1), ecv, nParVerts, v_assoc_pe(ecv)
      
      ! add additional vertices on this face which lie on the plane
      ! stop when a complete face is formed
      !write(*,*) 'edge: ',[face_vid(1),face_vid(cvid-1)]
      v = v+1; pv = 1
      do while (.not.containsPair([face_vid(1),face_vid(cvid-1)],edge_vid) .and. pv<=nV)
        if (side(par_face_vid(modulo(v-1,nV)+1))==0) then
          face_vid(cvid) = p2c_vid(par_face_vid(modulo(v-1,nV)+1))
          cvid = cvid+1
          !write(*,*) '4',face_vid(cvid-1)
        end if
        v = v+1; pv = pv+1
      end do
      if (pv==nV .and. .not.containsPair([face_vid(1),face_vid(cvid-1)],edge_vid)) then
        write(*,*) 'final pair not found: ',[face_vid(1),face_vid(cvid-1)]
        ierr = 1
      end if

    end subroutine modified_parent_face

  end function polyhedron_on_side_of_plane

  subroutine print_data (this,normalized)

    class(polyhedron), intent(in) :: this
    logical, optional, intent(in) :: normalized

    integer :: v,e,f
    real(r8) :: x0(3),xl(3)
    logical :: normalizedh

    normalizedh = merge(normalized, .false., present(normalized))
    
    write(*,*) 'POLYHEDRON DATA:'
    
    if (allocated(this%x)) then
      if (normalizedh) then
        x0 = minval(this%x,dim=2)
        xl = maxval(this%x,dim=2) - x0
        do v = 1,this%nVerts
          write(*,'(a,i3,a,3es25.15)') 'x ',v,':  ',(this%x(:,v)-x0)/xl
        end do
      else
        do v = 1,this%nVerts
          write(*,'(a,i3,a,3es35.25)') 'x ',v,':  ',this%x(:,v)
        end do
      end if
      write(*,*)
    end if
    
    if (allocated(this%edge_vid)) then
      do e = 1,this%nEdges
        write(*,'(a,i3,a,2i4)') 'edge ',e,':  ',this%edge_vid(:,e)
      end do
      write(*,*)
    end if
    
    if (allocated(this%face_vid)) then
      do f = 1,this%nFaces
        write(*,'(a,i3,a,10i4)') 'face ',f,':  ',this%face_vid(:,f)
      end do
      write(*,*)
    end if

    if (allocated(this%face_normal)) then
      do f = 1,this%nFaces
        write(*,'(a,i3,a,3es35.25)') 'norm ',f,':  ',this%face_normal(:,f)
      end do
      write(*,*)
    end if

    write(*,'(a,es20.10)') 'volume ',this%vol

  end subroutine print_data

end module polyhedron_type
