#include "f90_assert.fpp"

module analytic_surface_type_test

  use kinds, only: r8
  use unstr_mesh_type
  use mesh_geom_type
  use analytic_surface_type
  use paraboloid_type
  use logging_services
  implicit none
  private

  public :: analytic_surface_test_suite, curvature_test_suite

contains

  subroutine analytic_surface_test_suite ()

    integer :: i, ncell
    real(r8), allocatable :: lnorm(:)

    print '(a)'
    print '(a)', 'ANALYTIC_SURFACE'
    print '(a)', '===================================================='

    ! call plane_test ()
    ! call parabola_test ()
    ! call messy_test ()
    ! call messy_test2 ()
    ! call messy_test3 ()
    ! call messy_test4 ()
    call messy_test5()
    ! print *

    print '(a)', '===================================================='
    print '(a)'

  end subroutine analytic_surface_test_suite

  subroutine curvature_test_suite ()

    use, intrinsic :: iso_fortran_env, only: output_unit
    use timer_tree_type

    integer :: i, ncell, fh1, fh2
    real(r8) :: lnormFT(3), lnormHF(3)

    !call set_timer_type (realtime_timing)

    print '(a)'
    print '(a)', 'CURVATURE'
    print '(a)', '===================================================='

    open(newunit=fh1, file="ft_conv_dump.txt")
    open(newunit=fh2, file="hf_conv_dump.txt")

    write (fh1, '(a)') '# dx l1 l2 l3'
    write (fh2, '(a)') '# dx l1 l2 l3'

    do i = 1,4
      ncell = 10 * 2**i
      !call mesh_2d_test (ncell, 'cylinder.json', lnormFT, lnormHF)
      call mesh_3d_test (ncell, 'ellipsoid.json', lnormFT, lnormHF)
      !call mesh_unstr_test (ncell, 'ellipsoid.json', lnormFT, lnormHF)
      write (fh1, '(4es15.5)') 1.0_r8 / ncell, lnormFT
      write (fh2, '(4es15.5)') 1.0_r8 / ncell, lnormHF
    end do

    close(fh1); close(fh2)

    print '(a)', '===================================================='
    print '(a)'

    call write_timer_tree (output_unit, indent=3)

  end subroutine curvature_test_suite

  subroutine plane_test ()

    real(r8), allocatable :: x(:,:)
    real(r8) :: dx
    integer :: N,ind,i,j
    type(analytic_surface) :: surf

    dx = 0.1_r8
    N = 9

    allocate(x(3,N*N))
    do i = 1,N
      do j = 1,N
        ind = i + (j-1)*N
        x(1,ind) = (i-N/2+1)*dx
        x(2,ind) = (j-N/2+1)*dx
        x(3,ind) = x(1,ind) + x(2,ind)
      end do
    end do

    call surf%bestFit (x)

#ifndef NAGFOR
    ! NAG doesn't support derived type i/o as of v6.1
    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature([0.0_r8,0.0_r8,0.0_r8])
#endif

  end subroutine plane_test

  subroutine parabola_test ()

    real(r8), allocatable :: x(:,:)
    real(r8) :: dx
    integer :: N,ind,i,j
    type(analytic_surface) :: surf

    dx = 0.1_r8
    N = 9

    allocate(x(3,N*N))
    do i = 1,N
      do j = 1,N
        ind = i + (j-1)*N
        x(1,ind) = (i-N/2+1)*dx
        x(2,ind) = (j-N/2+1)*dx
        x(3,ind) = x(1,ind)**2 + x(2,ind)**2 ! parabola
      end do
    end do

    call surf%bestFit (x)

#ifndef NAGFOR
    ! NAG doesn't support derived type i/o as of v6.1
    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature([0.0_r8,0.0_r8,0.0_r8])
#endif

  end subroutine parabola_test

  subroutine messy_test ()

    real(r8), allocatable :: x(:,:)
    real(r8) :: dx
    integer :: N,ind,i,j
    type(analytic_surface) :: surf

    dx = 0.1_r8
    N = 9

    x = reshape([ 3.3765E-03,    2.5013E-01,    0.0000E+00,&
                  1.6882E-03,    2.5031E-01,    1.3532E-02,&
                  1.6882E-03,    2.5031E-01,   -1.3532E-02,&
                 -4.9481E-02,    2.4442E-01,   -6.2500E-02,&
                 -2.9952E-02,    2.4721E-01,   -4.8968E-02,&
                 -2.9952E-02,    2.4721E-01,   -7.6032E-02,&
                  2.3358E-02,    2.4814E-01,   -6.2500E-02,&
                  4.2929E-02,    2.4534E-01,   -4.8968E-02,&
                  4.2929E-02,    2.4534E-01,   -7.6032E-02,&
                  9.3750E-02,    2.3076E-01,   -7.8125E-02,&
                  1.0728E-01,    2.2496E-01,   -5.4688E-02,&
                  8.0218E-02,    2.3656E-01,   -5.4688E-02,&
                 -6.2148E-03,    2.5023E-01,   -6.2500E-02,&
                 -3.1074E-03,    2.5058E-01,   -7.6032E-02,&
                 -3.1074E-03,    2.5058E-01,   -4.8968E-02,&
                  3.3765E-03,    2.5013E-01,   -6.2500E-02,&
                  1.6882E-03,    2.5031E-01,   -4.8968E-02,&
                  1.6882E-03,    2.5031E-01,   -7.6032E-02,&
                 -4.9481E-02,    2.4442E-01,    0.0000E+00,&
                 -2.9952E-02,    2.4721E-01,    1.3532E-02,&
                 -2.9952E-02,    2.4721E-01,   -1.3532E-02,&
                  2.3358E-02,    2.4814E-01,    0.0000E+00,&
                  4.2929E-02,    2.4534E-01,    1.3532E-02,&
                  4.2929E-02,    2.4534E-01,   -1.3532E-02,&
                  9.3750E-02,    2.3076E-01,   -1.5625E-02,&
                  1.0728E-01,    2.2496E-01,    7.8125E-03,&
                  8.0218E-02,    2.3656E-01,    7.8125E-03,&
                 -6.2148E-03,    2.5023E-01,    0.0000E+00,&
                 -3.1074E-03,    2.5058E-01,   -1.3532E-02,&
                 -3.1074E-03,    2.5058E-01,    1.3532E-02,&
                 -4.9481E-02,    2.4442E-01,    6.2500E-02,&
                 -2.9952E-02,    2.4721E-01,    7.6032E-02,&
                 -2.9952E-02,    2.4721E-01,    4.8968E-02,&
                  2.3358E-02,    2.4814E-01,    6.2500E-02,&
                  4.2929E-02,    2.4534E-01,    7.6032E-02,&
                  4.2929E-02,    2.4534E-01,    4.8968E-02,&
                  9.3750E-02,    2.3076E-01,    4.6875E-02,&
                  1.0728E-01,    2.2496E-01,    7.0312E-02,&
                  8.0218E-02,    2.3656E-01,    7.0312E-02,&
                 -6.2148E-03,    2.5023E-01,    6.2500E-02,&
                 -3.1074E-03,    2.5058E-01,    4.8968E-02,&
                 -3.1074E-03,    2.5058E-01,    7.6032E-02,&
                  3.3765E-03,    2.5013E-01,    6.2500E-02,&
                  1.6882E-03,    2.5031E-01,    7.6032E-02,&
                  1.6882E-03,    2.5031E-01,    4.8968E-02], [3,45])

    call surf%bestFit (x)

#ifndef NAGFOR
    ! NAG doesn't support derived type i/o as of v6.1
    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature(sum(x(:,1:3), dim=2) / 3.0_r8)
#endif

  end subroutine messy_test

  subroutine messy_test2 ()

    real(r8), allocatable :: x(:,:)
    real(r8) :: dx
    integer :: N,ind,i,j
    type(analytic_surface) :: surf

    dx = 0.1_r8
    N = 9

    x = reshape([&
        -3.0328E-01,   -1.8111E-01,    0.0000E+00,&
        -3.0789E-01,   -1.7153E-01,   -1.3532E-02,&
        -3.0789E-01,   -1.7153E-01,    1.3532E-02,&
        -2.7379E-01,   -2.1738E-01,   -4.6875E-02,&
        -2.6349E-01,   -2.3032E-01,   -7.0312E-02,&
        -2.8409E-01,   -2.0444E-01,   -7.0312E-02,&
        -2.4913E-01,   -2.4742E-01,   -6.2500E-02,&
        -2.4783E-01,   -2.4871E-01,   -7.6032E-02,&
        -2.4783E-01,   -2.4871E-01,   -4.8968E-02,&
        -3.1624E-01,   -1.5081E-01,   -6.2500E-02,&
        -3.2185E-01,   -1.3791E-01,   -7.6032E-02,&
        -3.2185E-01,   -1.3791E-01,   -4.8968E-02,&
        -3.0328E-01,   -1.8111E-01,   -6.2500E-02,&
        -3.0789E-01,   -1.7153E-01,   -7.6032E-02,&
        -3.0789E-01,   -1.7153E-01,   -4.8968E-02,&
        -3.3724E-01,   -9.3750E-02,   -7.8125E-02,&
        -3.4086E-01,   -8.0218E-02,   -5.4688E-02,&
        -3.3362E-01,   -1.0728E-01,   -5.4688E-02,&
        -2.7379E-01,   -2.1738E-01,    1.5625E-02,&
        -2.6349E-01,   -2.3032E-01,   -7.8125E-03,&
        -2.8409E-01,   -2.0444E-01,   -7.8125E-03,&
        -2.4913E-01,   -2.4742E-01,    0.0000E+00,&
        -2.4783E-01,   -2.4871E-01,   -1.3532E-02,&
        -2.4783E-01,   -2.4871E-01,    1.3532E-02,&
        -3.2372E-01,   -1.3360E-01,    0.0000E+00,&
        -3.1811E-01,   -1.4651E-01,   -1.3532E-02,&
        -3.1811E-01,   -1.4651E-01,    1.3532E-02,&
        -3.3724E-01,   -9.3750E-02,   -1.5625E-02,&
        -3.4086E-01,   -8.0218E-02,    7.8125E-03,&
        -3.3362E-01,   -1.0728E-01,    7.8125E-03,&
        -2.7379E-01,   -2.1738E-01,    7.8125E-02,&
        -2.6349E-01,   -2.3032E-01,    5.4688E-02,&
        -2.8409E-01,   -2.0444E-01,    5.4687E-02,&
        -2.4913E-01,   -2.4742E-01,    6.2500E-02,&
        -2.4783E-01,   -2.4871E-01,    4.8968E-02,&
        -2.4783E-01,   -2.4871E-01,    7.6032E-02,&
        -3.1624E-01,   -1.5081E-01,    6.2500E-02,&
        -3.2185E-01,   -1.3791E-01,    4.8968E-02,&
        -3.2185E-01,   -1.3791E-01,    7.6032E-02,&
        -3.0328E-01,   -1.8111E-01,    6.2500E-02,&
        -3.0789E-01,   -1.7153E-01,    4.8968E-02,&
        -3.0789E-01,   -1.7153E-01,    7.6032E-02,&
        -3.3724E-01,   -9.3750E-02,    4.6875E-02,&
        -3.4086E-01,   -8.0218E-02,    7.0312E-02,&
        -3.3362E-01,   -1.0728E-01,    7.0312E-02], [3,45])

    !call surf%bestFit (x)
    !call surf%bestParaboloidFit (x)
    call surf%bestOneSheetFit (x)

#ifndef NAGFOR
    ! NAG doesn't support derived type i/o as of v6.1
    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature(sum(x(:,1:3), dim=2) / 3.0_r8)
#endif

  end subroutine messy_test2

  subroutine messy_test3 ()

    real(r8), allocatable :: x(:,:), weight(:)
    real(r8) :: dx
    integer :: N,ind,i,j
    !type(analytic_surface) :: surf
    type(paraboloid) :: surf

    dx = 0.1_r8
    N = 9

    x = reshape([&
        -1.5187991043E-01_r8,   -2.0048700559E-01_r8,    0.0000000000E+00_r8,&
        -1.5093995521E-01_r8,   -2.0121751398E-01_r8,   -1.0825317547E-02_r8,&
        -1.5093995521E-01_r8,   -2.0121751398E-01_r8,    1.0825317547E-02_r8,&
        -1.5187991043E-01_r8,   -2.0048700559E-01_r8,   -5.0000000000E-02_r8,&
        -1.5093995521E-01_r8,   -2.0121751398E-01_r8,   -6.0825317547E-02_r8,&
        -1.5093995521E-01_r8,   -2.0121751398E-01_r8,   -3.9174682453E-02_r8,&
        -1.2500000000E-01_r8,   -2.1599817276E-01_r8,   -6.2500000000E-02_r8,&
        -1.3582531755E-01_r8,   -2.1000567662E-01_r8,   -4.3750000000E-02_r8,&
        -1.1417468245E-01_r8,   -2.2199066890E-01_r8,   -4.3750000000E-02_r8,&
        -2.0048795062E-01_r8,   -1.5187626973E-01_r8,   -5.0000000000E-02_r8,&
        -2.0121987654E-01_r8,   -1.5093813487E-01_r8,   -6.0825317547E-02_r8,&
        -2.0121987654E-01_r8,   -1.5093813487E-01_r8,   -3.9174682453E-02_r8,&
        -1.7636193379E-01_r8,   -1.7641207895E-01_r8,   -6.2500000000E-02_r8,&
        -1.8659751671E-01_r8,   -1.6619820953E-01_r8,   -4.3750000000E-02_r8,&
        -1.6612635087E-01_r8,   -1.8662594838E-01_r8,   -4.3750000000E-02_r8,&
        -1.2500000000E-01_r8,   -2.1599817276E-01_r8,   -1.2500000000E-02_r8,&
        -1.3582531755E-01_r8,   -2.1000567662E-01_r8,    6.2500000000E-03_r8,&
        -1.1417468245E-01_r8,   -2.2199066890E-01_r8,    6.2500000000E-03_r8,&
        -2.0146385185E-01_r8,   -1.5062542324E-01_r8,    0.0000000000E+00_r8,&
        -2.0073192593E-01_r8,   -1.5156355811E-01_r8,   -1.0825317547E-02_r8,&
        -2.0073192593E-01_r8,   -1.5156355811E-01_r8,    1.0825317547E-02_r8,&
        -1.7636193379E-01_r8,   -1.7641207895E-01_r8,    1.2500000000E-02_r8,&
        -1.6612635087E-01_r8,   -1.8662594838E-01_r8,   -6.2500000000E-03_r8,&
        -1.8659751671E-01_r8,   -1.6619820953E-01_r8,   -6.2500000000E-03_r8,&
        -1.5062663681E-01_r8,   -2.0146101677E-01_r8,    5.0000000000E-02_r8,&
        -1.5156659202E-01_r8,   -2.0073050839E-01_r8,    6.0825317547E-02_r8,&
        -1.5156659202E-01_r8,   -2.0073050839E-01_r8,    3.9174682453E-02_r8,&
        -1.2500000000E-01_r8,   -2.1599817276E-01_r8,    3.7500000000E-02_r8,&
        -1.3582531755E-01_r8,   -2.1000567662E-01_r8,    5.6250000000E-02_r8,&
        -1.1417468245E-01_r8,   -2.2199066890E-01_r8,    5.6250000000E-02_r8,&
        -2.0048795062E-01_r8,   -1.5187626973E-01_r8,    5.0000000000E-02_r8,&
        -2.0121987654E-01_r8,   -1.5093813487E-01_r8,    6.0825317547E-02_r8,&
        -2.0121987654E-01_r8,   -1.5093813487E-01_r8,    3.9174682453E-02_r8,&
        -1.7636193379E-01_r8,   -1.7641207895E-01_r8,    3.7500000000E-02_r8,&
        -1.8659751671E-01_r8,   -1.6619820953E-01_r8,    5.6250000000E-02_r8,&
        -1.6612635087E-01_r8,   -1.8662594838E-01_r8,    5.6250000000E-02_r8], [3,36])

    allocate(weight(size(x,dim=2)))
    weight = 1

    call surf%bestFit (x, weight)
    !call surf%bestParaboloidFit (x)
    !call surf%bestOneSheetFit (x)

#ifndef NAGFOR
    ! NAG doesn't support derived type i/o as of v6.1
    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature(sum(x(:,1:3), dim=2) / 3)
#endif

  end subroutine messy_test3

  subroutine messy_test4 ()

    real(r8), allocatable :: x(:,:), weight(:)
    real(r8) :: dx
    integer :: N,ind,i,j
    !type(analytic_surface) :: surf
    type(paraboloid) :: surf

    dx = 0.1_r8
    N = 9

    x = reshape([&
        -7.1511535192E-03_r8,   -2.5012802513E-01_r8,    0.0000000000E+00_r8,&
        -3.5755767596E-03_r8,   -2.5032006282E-01_r8,   -1.0825317547E-02_r8,&
        -3.5755767596E-03_r8,   -2.5032006282E-01_r8,    1.0825317547E-02_r8,&
        -2.3837178397E-03_r8,   -2.5038407538E-01_r8,   -5.0000000000E-02_r8,&
        -5.9592945993E-03_r8,   -2.5019203769E-01_r8,   -6.0825317547E-02_r8,&
        -5.9592945993E-03_r8,   -2.5019203769E-01_r8,   -3.9174682453E-02_r8,&
        6.8663166500E-03_r8,   -2.5013333602E-01_r8,   -5.0000000000E-02_r8,&
        3.4331583250E-03_r8,   -2.5033334005E-01_r8,   -6.0825317547E-02_r8,&
        3.4331583250E-03_r8,   -2.5033334005E-01_r8,   -3.9174682453E-02_r8,&
        -7.5000000000E-02_r8,   -2.3808064461E-01_r8,   -6.2500000000E-02_r8,&
        -8.5825317547E-02_r8,   -2.3460907085E-01_r8,   -4.3750000000E-02_r8,&
        -6.4174682453E-02_r8,   -2.4155221837E-01_r8,   -4.3750000000E-02_r8,&
        -3.8749783245E-02_r8,   -2.4672498711E-01_r8,   -5.0000000000E-02_r8,&
        -2.1874458113E-02_r8,   -2.4836249355E-01_r8,   -6.0825317547E-02_r8,&
        -2.1874458113E-02_r8,   -2.4836249355E-01_r8,   -3.9174682453E-02_r8,&
        1.5279179784E-02_r8,   -2.4878639736E-01_r8,   -5.0000000000E-02_r8,&
        3.2639589892E-02_r8,   -2.4696599340E-01_r8,   -6.0825317547E-02_r8,&
        3.2639589892E-02_r8,   -2.4696599340E-01_r8,   -3.9174682453E-02_r8,&
        6.8663166500E-03_r8,   -2.5013333602E-01_r8,    0.0000000000E+00_r8,&
        3.4331583250E-03_r8,   -2.5033334005E-01_r8,   -1.0825317547E-02_r8,&
        3.4331583250E-03_r8,   -2.5033334005E-01_r8,    1.0825317547E-02_r8,&
        -7.5000000000E-02_r8,   -2.3808064461E-01_r8,   -1.2500000000E-02_r8,&
        -8.5825317547E-02_r8,   -2.3460907085E-01_r8,    6.2500000000E-03_r8,&
        -6.4174682453E-02_r8,   -2.4155221837E-01_r8,    6.2500000000E-03_r8,&
        -3.8749783245E-02_r8,   -2.4672498711E-01_r8,    0.0000000000E+00_r8,&
        -2.1874458113E-02_r8,   -2.4836249355E-01_r8,   -1.0825317547E-02_r8,&
        -2.1874458113E-02_r8,   -2.4836249355E-01_r8,    1.0825317547E-02_r8,&
        3.8426393261E-02_r8,   -2.4635919209E-01_r8,    0.0000000000E+00_r8,&
        2.1065983154E-02_r8,   -2.4817959604E-01_r8,   -1.0825317547E-02_r8,&
        2.1065983154E-02_r8,   -2.4817959604E-01_r8,    1.0825317547E-02_r8,&
        -2.3837178397E-03_r8,   -2.5038407538E-01_r8,    5.0000000000E-02_r8,&
        -5.9592945993E-03_r8,   -2.5019203769E-01_r8,    6.0825317547E-02_r8,&
        -5.9592945993E-03_r8,   -2.5019203769E-01_r8,    3.9174682453E-02_r8,&
        6.8663166500E-03_r8,   -2.5013333602E-01_r8,    5.0000000000E-02_r8,&
        3.4331583250E-03_r8,   -2.5033334005E-01_r8,    6.0825317547E-02_r8,&
        3.4331583250E-03_r8,   -2.5033334005E-01_r8,    3.9174682453E-02_r8,&
        -7.5000000000E-02_r8,   -2.3808064461E-01_r8,    3.7500000000E-02_r8,&
        -8.5825317547E-02_r8,   -2.3460907085E-01_r8,    5.6250000000E-02_r8,&
        -6.4174682453E-02_r8,   -2.4155221837E-01_r8,    5.6250000000E-02_r8,&
        -3.8749783245E-02_r8,   -2.4672498711E-01_r8,    5.0000000000E-02_r8,&
        -2.1874458113E-02_r8,   -2.4836249355E-01_r8,    6.0825317547E-02_r8,&
        -2.1874458113E-02_r8,   -2.4836249355E-01_r8,    3.9174682453E-02_r8,&
        1.5279179784E-02_r8,   -2.4878639736E-01_r8,    5.0000000000E-02_r8,&
        3.2639589892E-02_r8,   -2.4696599340E-01_r8,    6.0825317547E-02_r8,&
        3.2639589892E-02_r8,   -2.4696599340E-01_r8,    3.9174682453E-02_r8], [3,45])

    allocate(weight(size(x,dim=2)))
    weight = 1

    call surf%bestFit (x, weight)
    !call surf%bestParaboloidFit (x)
    !call surf%bestOneSheetFit (x)

#ifndef NAGFOR
    ! NAG doesn't support derived type i/o as of v6.1
    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature(sum(x(:,1:3), dim=2) / 3)
#endif

  end subroutine messy_test4

  subroutine messy_test5 ()

    integer :: i
    real(r8), allocatable :: x(:,:), weight(:)
    real(r8) :: curvature, R
    type(paraboloid) :: surf

    R = 0.35_r8

    x = reshape([&
        6.4153003891E-02_r8,    3.4406901123E-01_r8,    0.0000000000E+00_r8, &
        6.7276826530E-02_r8,    3.4347202777E-01_r8,   -6.2500000000E-03_r8, &
        7.1875000000E-02_r8,    3.4253552490E-01_r8,   -6.2500000000E-03_r8, &
        5.9375000000E-02_r8,    3.4492210380E-01_r8,   -6.2500000000E-03_r8, &
        6.4153003885E-02_r8,    3.4406901123E-01_r8,   -6.2500000000E-03_r8, &
        6.7276826531E-02_r8,    3.4347202777E-01_r8,    0.0000000000E+00_r8, &
        7.1875000000E-02_r8,    3.4253552490E-01_r8,    0.0000000000E+00_r8, &
        5.9375000000E-02_r8,    3.4492210379E-01_r8,    0.0000000000E+00_r8, &
        6.7276826531E-02_r8,    3.4347202777E-01_r8,    6.2500000000E-03_r8, &
        7.1875000000E-02_r8,    3.4253552490E-01_r8,    6.2500000000E-03_r8, &
        5.9375000000E-02_r8,    3.4492210380E-01_r8,    6.2500000000E-03_r8, &
        6.4153003884E-02_r8,    3.4406901123E-01_r8,    6.2500000000E-03_r8], [3,12])

    ! x = reshape([&
    !     6.4153003891E-02_r8,    3.4406901123E-01_r8,    0.0000000000E+00_r8, &
    !     6.7276826531E-02_r8,    3.4347202777E-01_r8,    0.0000000000E+00_r8, &
    !     !7.1875000000E-02_r8,    3.4253552490E-01_r8,    0.0000000000E+00_r8, &
    !     5.9375000000E-02_r8,    3.4492210379E-01_r8,    0.0000000000E+00_r8, &
    !     6.4153003885E-02_r8,    3.4406901123E-01_r8,   -6.2500000000E-03_r8, &
    !     6.7276826530E-02_r8,    3.4347202777E-01_r8,   -6.2500000000E-03_r8, &
    !     !7.1875000000E-02_r8,    3.4253552490E-01_r8,   -6.2500000000E-03_r8, &
    !     5.9375000000E-02_r8,    3.4492210380E-01_r8,   -6.2500000000E-03_r8, &
    !     6.4153003884E-02_r8,    3.4406901123E-01_r8,    6.2500000000E-03_r8, &
    !     6.7276826531E-02_r8,    3.4347202777E-01_r8,    6.2500000000E-03_r8, &
    !     !7.1875000000E-02_r8,    3.4253552490E-01_r8,    6.2500000000E-03_r8, &
    !     5.9375000000E-02_r8,    3.4492210380E-01_r8,    6.2500000000E-03_r8 & !, &
    !     ], [3,9])

    do i = 1,size(x,dim=2)
      print *, 'd: ', abs(norm2(x(:2,i)) - R)

      if (abs(norm2(x(:2,i)) - R) > 4e-6) &
          x(2,i) = sqrt(R**2 - x(1,i)**2)
    end do

    allocate(weight(size(x,dim=2)))
    weight = 1

    call surf%bestFit (x, weight, [-1.89492268E-01_r8, -9.81882213E-01_r8, 9.46296280E-08_r8])
    curvature = surf%curvature(x(:,1))

#ifndef NAGFOR
    ! NAG doesn't support derived type i/o as of v6.1
    print '(dt,2(a,es15.5))', surf, ',     curvature: ', curvature, &
        ',   err: ', abs(curvature + 1/R) * R
#endif

  end subroutine messy_test5

  subroutine mesh_2d_test (mesh_size, shape_filename, lnormFT, lnormHF)

    use,intrinsic :: iso_c_binding, only: C_NEW_LINE
    use consts, only: cutvof
    use unstr_mesh_factory
    use surface_type
    use parameter_list_type
    use parameter_list_json
    use int_norm_module
    use interface_patch_type
    use vof_init
    use multimat_cell_type
    use array_utils, only: normalize, isZero
    use curvature_hf
    use vof_io
    use mixed_cell_subset_constructor
    use mesh_subset_type
    use vof_init_ex_circle

    integer, intent(in) :: mesh_size
    character(*), intent(in) :: shape_filename
    real(r8), intent(out) :: lnormFT(:), lnormHF(:)

    integer, parameter :: thck = 5
    character(:), allocatable :: errmsg, filename
    character(30) :: tmp
    type(unstr_mesh) :: mesh
    type(mesh_geom) :: gmesh
    type(mesh_subset) :: mixed_cells
    type(parameter_list), pointer :: plist
    real(r8), allocatable :: curvature_ex(:)
    real(r8) :: vof(2,mesh_size*mesh_size*thck), curvature(mesh_size*mesh_size*thck), &
        vofex(2,mesh_size*mesh_size*thck)
    real(r8), allocatable :: int_norm(:,:,:)
    integer :: infile

    ! create a regular 2D mesh
    mesh = new_unstr_mesh ([-0.5_r8, -0.5_r8, -thck*0.5_r8/mesh_size], &
        [0.5_r8, 0.5_r8, thck*0.5_r8 / mesh_size], [mesh_size,mesh_size,thck])
    call gmesh%init (mesh)

    ! fill the mesh with volume fractions for a circle
    ! right now this relies on an input file to describe the cylinder. In the future I'd like to
    ! initialize material_geometry_types without a parameter_list_type, or initialize a
    ! parameter_list_type without a JSON file
    open(newunit=infile,file=shape_filename,action='read',access='stream')
    call parameter_list_from_json_stream (infile, plist, errmsg)
    if (.not.associated(plist)) call LS_fatal ("error reading input file:" // C_NEW_LINE // errmsg)
    close(infile)

    plist => plist%sublist('initial-vof')

    ! write (tmp, '(a,i0,a)') "vof_2d_", mesh_size, ".dat"
    ! filename = trim(adjustl(tmp))
    ! if (file_exists(filename)) then
    !   call read_vof_field(filename, vof)
    ! else
    !   call vof_initialize (mesh, plist, vof, [1,2], 2)
    !   call store_vof_field(filename, vof)
    ! end if
    deallocate(plist)

    call vof_init_circle(mesh, 0.35_r8, vofex)
    !call compare_vof_fields(vof, vofex)
    vof = vofex

    ! get mixed cells subset
    mixed_cells = mixed_cell_subset(vof, mesh)

    ! get the interface reconstructions
    int_norm = interface_normal (vof, mesh, gmesh, .false.)

    ! calculate errors for FT and HF curvature methods
    lnormFT = ft_mesh_test ('reg', vof, int_norm, mesh, gmesh, mixed_cells, curvature_ex)
    lnormHF = hf_mesh_test (vof, int_norm, mesh, gmesh, curvature_ex)

    print '(i5, 2(a,3es10.2))', mesh_size, '  FT L1,L2,Linf = ',lnormFT, &
        ',  HF L1,L2,Linf = ',lnormHF

  end subroutine mesh_2d_test

  subroutine mesh_3d_test (mesh_size, shape_filename, lnormFT, lnormHF)

    use,intrinsic :: iso_c_binding, only: C_NEW_LINE
    use consts, only: cutvof
    use unstr_mesh_factory
    use surface_type
    use parameter_list_type
    use parameter_list_json
    use int_norm_module
    use interface_patch_type
    use vof_init
    use multimat_cell_type
    use array_utils, only: normalize, isZero
    use curvature_hf
    use vof_io
    use mixed_cell_subset_constructor
    use mesh_subset_type

    integer, intent(in) :: mesh_size
    character(*), intent(in) :: shape_filename
    real(r8), intent(out) :: lnormFT(:), lnormHF(:)

    character(:), allocatable :: errmsg, filename
    character(30) :: tmp
    type(unstr_mesh) :: mesh
    type(mesh_geom) :: gmesh
    type(mesh_subset) :: mixed_cells
    type(parameter_list), pointer :: plist
    real(r8), target :: vof(2,mesh_size**3)
    real(r8), allocatable :: int_norm(:,:,:), curvature_ex(:)
    integer :: infile

    ! create a regular 2D mesh
    mesh = new_unstr_mesh ([-0.5_r8, -0.5_r8, -0.5_r8], &
        [0.5_r8, 0.5_r8, 0.5_r8], [mesh_size,mesh_size,mesh_size])
    call gmesh%init (mesh)
    print '(a)', 'mesh initialized'

    ! fill the mesh with volume fractions for a circle
    ! right now this relies on an input file to describe the cylinder. In the future I'd like to
    ! initialize material_geometry_types without a parameter_list_type, or initialize a
    ! parameter_list_type without a JSON file
    open(newunit=infile,file=shape_filename,action='read',access='stream')
    call parameter_list_from_json_stream (infile, plist, errmsg)
    if (.not.associated(plist)) call LS_fatal ("error reading input file:" // C_NEW_LINE // errmsg)
    close(infile)

    plist => plist%sublist('initial-vof')

    write (tmp, '(a,i0,a)') "vof_3d_", mesh_size, ".dat"
    filename = trim(adjustl(tmp))
    if (file_exists(filename)) then
      call read_vof_field(filename, vof)
    else
      call vof_initialize (mesh, gmesh, plist, vof, [1,2], 2)
      call store_vof_field(filename, vof)
    end if
    deallocate(plist)
    print '(a)', 'vof initialized'

    ! get mixed cells subset
    mixed_cells = mixed_cell_subset(vof, mesh)

    ! get the interface reconstructions
    int_norm = interface_normal (vof, mesh, gmesh, .false.)

    ! calculate errors for FT and HF curvature methods
    lnormFT = ft_mesh_test ('reg', vof, int_norm, mesh, gmesh, mixed_cells, curvature_ex)
    lnormHF = hf_mesh_test (vof, int_norm, mesh, gmesh, curvature_ex)

    print '(i5, 2(a,3es10.2))', mesh_size, '  FT L1,L2,Linf = ',lnormFT, &
        ',  HF L1,L2,Linf = ',lnormHF

  end subroutine mesh_3d_test

  subroutine mesh_unstr_test (mesh_size, shape_filename, lnormFT, lnormHF)

    use,intrinsic :: iso_c_binding, only: C_NEW_LINE
    use consts, only: cutvof
    use unstr_mesh_factory
    use surface_type
    use parameter_list_type
    use parameter_list_json
    use int_norm_module
    use interface_patch_type
    use vof_init
    use multimat_cell_type
    use array_utils, only: normalize, isZero
    use curvature_hf
    use vof_io
    use mixed_cell_subset_constructor
    use mesh_subset_type

    integer, intent(in) :: mesh_size
    character(*), intent(in) :: shape_filename
    real(r8), intent(out) :: lnormFT(:), lnormHF(:)

    character(:), allocatable :: errmsg, filename, mesh_filename, cell_type
    character(30) :: tmp
    type(unstr_mesh) :: mesh
    type(mesh_geom) :: gmesh
    type(mesh_subset) :: mixed_cells
    type(parameter_list), pointer :: plist
    real(r8), allocatable, target :: vof(:,:)
    real(r8), allocatable :: curvature(:), curvature_ex(:)
    integer :: infile

    ! create a regular 2D mesh
    cell_type = 'tet' ! rnd, tet
    write (tmp, '(a,i0,3a)') "cube_", mesh_size, "_", cell_type, ".exo"
    mesh_filename = trim(adjustl(tmp))
    mesh = new_unstr_mesh (mesh_filename)
    call gmesh%init (mesh)
    call recalculate_mesh_volumes(mesh,gmesh)
    allocate(vof(2,mesh%ncell), curvature(mesh%ncell))
    print '(a,i0)', "mesh initialized ", mesh%ncell

    ! fill the mesh with volume fractions for a circle
    ! right now this relies on an input file to describe the cylinder. In the future I'd like to
    ! initialize material_geometry_types without a parameter_list_type, or initialize a
    ! parameter_list_type without a JSON file
    open(newunit=infile,file=shape_filename,action='read',access='stream')
    call parameter_list_from_json_stream (infile, plist, errmsg)
    if (.not.associated(plist)) call LS_fatal ("error reading input file:" // C_NEW_LINE // errmsg)
    close(infile)

    plist => plist%sublist('initial-vof')
    print '(a)', "parameter list read"

    write (tmp, '(a,i0,3a)') "vof_unstr_", mesh_size, "_", cell_type, ".dat"
    filename = trim(adjustl(tmp))
    if (file_exists(filename)) then
      call read_vof_field(filename, vof)
    else
      call vof_initialize (mesh, gmesh, plist, vof, [1,2], 2)
      call store_vof_field(filename, vof)
    end if
    !deallocate(plist)

    ! get mixed cells subset
    mixed_cells = mixed_cell_subset(vof, mesh)

    ! calculate errors for FT and HF curvature methods
    lnormFT = ft_unstr_mesh_test(cell_type, vof, mesh, gmesh, mixed_cells, curvature_ex)
    lnormHF = 0

    print '(i8,a,3es10.2)', mesh%ncell, '  FT L1,L2,Linf = ',lnormFT

  end subroutine mesh_unstr_test

  function ft_mesh_test (cell_type, vof, int_norm, mesh, gmesh, mixed_cells, curvature_ex) &
      result(lnorm)

    use consts, only: cutvof
    use surface_type
    use int_norm_module
    use interface_patch_type
    use multimat_cell_type
    use hex_types, only: hex_f, hex_e
    use array_utils, only: normalize, isZero
    use curvature_hf, only: heightFunction
    use lvira_normals
    use fit_normals
    use vof_io
    use mesh_subset_type

    character(*), intent(in) :: cell_type
    real(r8), intent(in) :: vof(:,:), int_norm(:,:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    type(mesh_subset), intent(in) :: mixed_cells
    real(r8), allocatable, intent(inout) :: curvature_ex(:)
    real(r8) :: lnorm(3)

    character(:), allocatable :: filename
    character(30) :: tmp
    integer :: i, nvofcell, ierr, imax, w, fh, c
    type(surface) :: intrec
    type(multimat_cell) :: cell
    real(r8) :: int_norm_local(3,2), err, curvature(mesh%ncell), wgt_scale, totvolume, xc(3)
    real(r8), allocatable :: int_norm_lvira(:,:,:), interface_centroid(:,:) !,weight_scales(:)

    !call heightFunction (throwaway, int_norm_hf, vof(1,:), int_norm(:,1,:), mesh, gmesh)

    write (tmp, '(a,i0,3a)') "normals_3d_", mesh%ncell, "_", cell_type, ".dat"
    filename = trim(adjustl(tmp))
    if (file_exists(filename)) then
      allocate(int_norm_lvira(3,2,mesh%ncell))
      call read_normals(filename, int_norm_lvira)
    else
      print '(a)', 'calculating normals ...'
      call interface_normals_lvira(int_norm_lvira, vof, mesh, gmesh, mixed_cells)
      call store_normals(filename, int_norm_lvira)
      print '(a)', 'done'
    end if
    !call interface_normals_fit(int_norm_lvira, vof, mesh, gmesh)

    ! 2d override. important on boundary cells
    ! print *, 'WARNING: 2D override'
    ! int_norm_hf(3,:) = 0
    ! do i = 1,mesh%ncell
    !   int_norm_hf(:,i) = normalize(int_norm_hf(:,i))
    ! end do

    ! do i = 1,mesh%ncell
    !   if (any(gmesh%cneighbor(:,i)<1)) then
    !     int_norm_hf(3,i) = 0
    !     int_norm_hf(:,i) = normalize(int_norm_hf(:,i))
    !   end if
    ! end do

    allocate(interface_centroid(3,mixed_cells%ncell))

    ! get the interface reconstructions
    lnorm = 0; nvofcell = 0; totvolume = 0
    do c = 1,mixed_cells%ncell
      i = mixed_cells%cell_id(c)

      call cell%init (ierr, i, mesh, gmesh, tesselate=.false.)
      ! call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, gmesh%outnorm(:,:,i), &
      !     mesh%volume(i), tesselate=.false.)
      if (ierr /= 0) call LS_fatal ('cell_outward_volflux failed: could not initialize cell')

      call cell%partition (vof(:,i), int_norm_lvira(:,:,i))

      call intrec%append (cell%interface_polygons(1), i)

      interface_centroid(:,c) = cell%interface_polygons(1)%centroid()

      err = abs(vof(1,i) - cell%mat_poly(1)%volume() / mesh%volume(i))
      totvolume = totvolume + mesh%volume(i)
      lnorm(1) = lnorm(1) + err * mesh%volume(i)
      lnorm(2) = lnorm(2) + err**2 * mesh%volume(i)
      lnorm(3) = max(lnorm(3),err)
    end do
    lnorm(1) = lnorm(1) / totvolume
    lnorm(2) = sqrt(lnorm(2) / totvolume)
    !print '(a,3es10.2)', 'VOF L1,L2,Linf = ',lnorm

    call intrec%write_ply("ft_test.ply")

    call compute_exact_curvature_field(curvature_ex, interface_centroid, int_norm_lvira, &
        mesh, gmesh, mixed_cells)
    ! allocate(curvature_ex(mesh%ncell))
    ! curvature_ex = 0

    ! get the curvature
    !weight_scales = [0.0_r8, 1.0_r8, 2.0_r8, 3.0_r8, 1.0_r8 / 2, 1.0_r8 / 3, 1.0_r8 / 4, 1.0_r8 / 8]
    open(newunit=fh, file="curv_err.txt")
    do w = 1,1
      !wgt_scale = (2.0_r8 ** (w - 1) - 1) * 3.0_r8 / 2.0_r8 ** 19
      wgt_scale = 0

      curvature = 0; lnorm = 0; nvofcell = 0; totvolume = 0
      !i = 18178
      do i = 1,mesh%ncell
        if (any(gmesh%cneighbor(:,i)<1)) cycle ! WARN: skipping boundaries. BCs might be automatic?
        !if (.not.isZero(gmesh%xc(3,i))) cycle ! WARN: 2d

        ! TODO: this really should be in any cell neighboring a cell containing the interface
        !if (vof(1,i) > cutvof .and. vof(1,i) < 1-cutvof) then !.and. .not.isZero(curvature(i))) then
        err = 0
        if (vof(1,i) > 1e-2_r8 .and. vof(1,i) < 1-1e-2_r8) then
          curvature(i) = abs(curvature_from_patch (intrec%local_patch(i,gmesh, vof(1,:)), &
              wgt_scale, int_norm_lvira(:,1,i), vof(1,:), mesh, gmesh, i, centroid=xc))
          if (isZero(curvature(i))) cycle

          ! curvature_ex(i) = ellipsoid_curvature_x(maxloc(abs(int_norm_lvira(:,1,i)),1), &
          !     xc, [0.7_r8, 0.5_r8, 0.3_r8])

          ! append to norms
          err = abs((curvature(i) - curvature_ex(i)) / curvature_ex(i))
          if (err > lnorm(3)) imax = i
          totvolume = totvolume + mesh%volume(i)
          lnorm(1) = lnorm(1) + err * mesh%volume(i)
          lnorm(2) = lnorm(2) + err**2 * mesh%volume(i)
          lnorm(3) = max(lnorm(3),err)

          ! if (err > 6.5e-2_r8) then
          ! !if (i==47291 .or. i==47131) then
          ! !if (i==47131) then
          ! !if (i==72891) then
          ! !if (i==64185) then
          ! ! if (i==68512) then

          ! !   call print_details(i)

          ! !   print *
          ! !   curvature(i) = abs(curvature_from_patch (intrec%local_patch(i,gmesh, vof(1,:)), &
          ! !       wgt_scale, int_norm_hf2(:,1,i), vof(1,:), mesh, gmesh, i, .true.))
          ! !   print *
          !   print '(i6, 3es14.4)', i, curvature(i), curvature_ex, err !, c_new_line
          !   print '(2es15.5)', vof(1,i), 1 - vof(1,i)
          !   print '(2es15.5)', gmesh%xc(1:2,i)
          !   print '(2es15.5)', minval(mesh%x(1,mesh%cnode(:,i))), minval(mesh%x(2,mesh%cnode(:,i)))
          !   print *, 'nvofs: ', vof(1,gmesh%cneighbor(:,i))
          !   ! print '(a, 3es14.4)', 'yn: ', int_norm(:,1,i)
          !   ! print '(a, 3es14.4)', 'hn: ', int_norm_hf2(:,1,i)
          !   print *
          !   !call LS_fatal ("large curvature error")
          ! end if
        end if
        !write (fh, '(2(es15.5,a),es15.5)') gmesh%xc(1,i),',', gmesh%xc(2,i),',', err
        write (fh, '(2(es15.5,a),es15.5)') &
            minval(mesh%x(1,mesh%cnode(:,i))),',', minval(mesh%x(2,mesh%cnode(:,i))),',', err
      end do
      lnorm(1) = lnorm(1) / totvolume
      lnorm(2) = sqrt(lnorm(2) / totvolume)
      close(fh)

      !print '(es10.2, a,3es10.2)', wgt_scale, '  FT L1,L2,Linf = ',lnorm
    end do

    !print *, 'maxloc: ', imax
    ! print '(i5, a,4es15.4)', 0.5_r8 - abs(vof(1,imax) - 0.5_r8)

  contains

    subroutine print_details (i)

      use polygon_type
      use array_utils, only: normalize, crossProduct

      integer, intent(in) :: i

      integer :: nc, sint, j
      type(polygon), allocatable :: interface_reconstruction(:)
      real(r8), allocatable :: centroid(:,:), normal(:,:)
      real(r8) :: x(3), R(3,3)

      !interface_reconstruction = intrec%local_patch(i,gmesh, vof(1,:))
      sint = size(interface_reconstruction)

      allocate(centroid(3,sint), normal(3,sint))

      nc = 0
      do j = 1,sint
        x = interface_reconstruction(j)%centroid2()

        if (isZero(x(3), 1e-6_r8)) then
          nc = nc + 1
          centroid(:,nc) = x
          normal(:,nc) = interface_reconstruction(j)%norm

          call interface_reconstruction(j)%print_data()
          print *
        end if
      end do
      print *

      print '(a,2es20.10)', 'curvature: ', curvature(i), err
      print *
      do j = 1,nc
        print '(a,3es20.10)', 'x: ', centroid(:,j)
      end do
      print *
      do j = 1,nc
        print '(a,3es20.10)', 'n: ', normal(:,j)
      end do


      R(3,:) = normal(:,1)
      R(2,:) = normalize(crossProduct([0.0_r8,0.0_r8,1.0_r8], normal(:,1)))
      R(1,:) = crossProduct(R(2,:), normal(:,1))

      print *
      do j = 1,nc
        print '(3es20.10)', matmul(R, centroid(:,j) - centroid(:,1))
      end do

    end subroutine print_details

  end function ft_mesh_test

  function ft_unstr_mesh_test (cell_type, vof, mesh, gmesh, mixed_cells, curvature_ex) result(lnorm)

    use consts, only: cutvof
    use surface_type
    use int_norm_module
    use interface_patch_type
    use multimat_cell_type
    use hex_types, only: hex_f, hex_e
    use array_utils, only: normalize, isZero
    use curvature_hf, only: heightFunction
    use lvira_normals
    use fit_normals
    use vof_io
    use polygon_type
    use mesh_subset_type
    use int_norm_module ! DEBUGGING

    character(*), intent(in) :: cell_type
    real(r8), intent(in) :: vof(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    type(mesh_subset), intent(in) :: mixed_cells
    real(r8), allocatable, intent(inout) :: curvature_ex(:)
    real(r8) :: lnorm(3)

    character(:), allocatable :: filename
    character(30) :: tmp
    integer :: i, ierr, imax, c
    type(surface) :: intrec
    type(multimat_cell) :: cell
    type(polygon_box), allocatable :: intrec_patch(:)
    real(r8) :: err, curvature(mesh%ncell), totvolume, xc(3)
    real(r8), allocatable :: int_norm_lvira(:,:,:), interface_centroid(:,:)

    write (tmp, '(a,i0,3a)') "normals_unstr_", mesh%ncell, "_", cell_type, ".dat"
    filename = trim(adjustl(tmp))
    if (file_exists(filename)) then
      allocate(int_norm_lvira(3,2,mesh%ncell))
      call read_normals(filename, int_norm_lvira)
    else
      call interface_normals_lvira(int_norm_lvira, vof, mesh, gmesh, mixed_cells)
      call store_normals(filename, int_norm_lvira)
    end if
    print '(a)', 'normals calculated'
    !call interface_normals_fit(int_norm_lvira, vof, mesh, gmesh)

    ! i = 435791
    ! int_norm_lvira(:,:,i) = interface_normal_cell(vof, i, mesh, gmesh, .false.)
    ! call interface_normal_lvira(int_norm_lvira(:,1,i), i, vof(1,:), mesh, gmesh)
    ! print '(a,3es20.10)', 'normal: ',int_norm_lvira(:,1,i)
    ! stop

    ! get the interface reconstructions
    lnorm = 0; totvolume = 0

    ! ! DEBUGGING ########
    ! allocate(int_norm_lvira(3,2,mesh%ncell))
    ! int_norm_lvira = 0
    ! !do i = 1300,1329
    ! i = 1329

    ! ! call interface_normal_lvira(int_norm_lvira(:,1,i), i, vof(1,:), mesh, gmesh)
    ! int_norm_lvira(:,1,i) = [1.0_r8, 0.0_r8, 0.0_r8]
    ! int_norm_lvira(:,2,i) = -int_norm_lvira(:,1,i)
    ! ! ##################

    allocate(interface_centroid(3,mixed_cells%ncell))

    do c = 1,mixed_cells%ncell
      i = mixed_cells%cell_id(c)

      !print *, 'here0'
      call cell%init (ierr, i, mesh, gmesh)
      ! call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, gmesh%outnorm(:,:,i), &
      !     mesh%volume(i))
      if (ierr /= 0) call LS_fatal ('cell_outward_volflux failed: could not initialize cell')
      !print *, 'here1'

      call cell%partition (vof(:,i), int_norm_lvira(:,:,i))
      !print *, 'here2'

      if (cell%interface_polygons(1)%n_elements < 1) then
        print *, 'vof ',i,vof(1,i)
        call LS_fatal("invalid cell partition")
      end if

      call intrec%append (cell%interface_polygons(1), i)
      !print *, 'here3'

      interface_centroid(:,c) = cell%interface_polygons(1)%centroid()

      err = abs(vof(1,i) - cell%mat_poly(1)%volume() / mesh%volume(i))
      totvolume = totvolume + mesh%volume(i)
      lnorm(1) = lnorm(1) + err * mesh%volume(i)
      lnorm(2) = lnorm(2) + err**2 * mesh%volume(i)
      lnorm(3) = max(lnorm(3),err)
    end do
    lnorm(1) = lnorm(1) / totvolume
    lnorm(2) = sqrt(lnorm(2) / totvolume)
    !print '(a,3es10.2)', 'VOF L1,L2,Linf = ',lnorm
    ! print *, 'here4'
    ! stop
    print '(a)', 'interfaces located'

    call intrec%write_ply("ft_test.ply")

    call compute_exact_curvature_field(curvature_ex, interface_centroid, int_norm_lvira, &
        mesh, gmesh, mixed_cells)
    ! allocate(curvature_ex(mesh%ncell))
    ! curvature_ex = 0

    ! get the curvature
    curvature = 0; lnorm = 0; totvolume = 0
    do i = 1,mesh%ncell
      if (any(gmesh%cneighbor(:,i)<1)) cycle ! WARN: skipping boundaries. BCs might be automatic?

      err = 0
      if (vof(1,i) > 1e-2_r8 .and. vof(1,i) < 1-1e-2_r8) then
        intrec_patch = intrec%local_patch(i,gmesh, vof(1,:))
        curvature(i) = abs(curvature_from_patch (intrec_patch, &
            0.0_r8, int_norm_lvira(:,1,i), vof(1,:), mesh, gmesh, i, centroid=xc))
        if (isZero(curvature(i))) cycle

        ! curvature_ex(i) = ellipsoid_curvature_x(maxloc(abs(int_norm_lvira(:,1,i)),1), &
        !     xc, [0.7_r8, 0.5_r8, 0.3_r8])

        ! append to norms
        err = abs((curvature(i) - curvature_ex(i)) / curvature_ex(i))
        if (err > lnorm(3)) imax = i
        totvolume = totvolume + mesh%volume(i)
        lnorm(1) = lnorm(1) + err * mesh%volume(i)
        lnorm(2) = lnorm(2) + err**2 * mesh%volume(i)
        lnorm(3) = max(lnorm(3),err)

        ! if (err > 2e0_r8) then
        !   intrec_patch = intrec%local_patch(i,gmesh, vof(1,:), verbose=.true.)
        !   curvature(i) = abs(curvature_from_patch (intrec_patch, &
        !       0.0_r8, int_norm_lvira(:,1,i), vof(1,:), mesh, gmesh, i, verbose=.true.))

        !   print '(i6, 3es14.4)', i, curvature(i), curvature_ex, err !, c_new_line
        !   call LS_fatal ("large curvature error")
        ! end if
      end if
    end do
    lnorm(1) = lnorm(1) / totvolume
    lnorm(2) = sqrt(lnorm(2) / totvolume)

  end function ft_unstr_mesh_test

  function hf_mesh_test (vof, int_norm, mesh, gmesh, curvature_ex) result(lnorm)

    use consts, only: cutvof
    use array_utils, only: normalize, isZero
    use curvature_hf

    real(r8), intent(in) :: vof(:,:), int_norm(:,:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8), intent(in) :: curvature_ex(:)
    real(r8) :: lnorm(3)

    integer :: i, nvofcell, ierr, imax
    real(r8) :: err, curvature(mesh%ncell)

    ! get the curvature
    curvature = abs(curvatureHF(vof(1,:), int_norm(:,1,:), mesh, gmesh))

    ! calculate error
    lnorm = 0; nvofcell = 0
    do i = 1,mesh%ncell
      if (any(gmesh%cneighbor(:,i)<1)) cycle ! WARN: skipping boundaries. BCs might be automatic?

      ! TODO: this really should be in any cell neighboring a cell containing the interface
      if (vof(1,i) > 1e-2_r8 .and. vof(1,i) < 1-1e-2_r8 .and. .not.isZero(curvature(i))) then
        ! append to norms
        ! curvature_ex = ellipsoid_curvature_x(maxloc(abs(int_norm(:,1,i)),1), gmesh%xc(:,i), &
        !     [0.7_r8, 0.5_r8, 0.3_r8])
        err = abs((curvature(i) - curvature_ex(i)) / curvature_ex(i))
        if (err > lnorm(3)) imax = i
        nvofcell = nvofcell + 1
        lnorm(1) = lnorm(1) + err
        lnorm(2) = lnorm(2) + err**2
        lnorm(3) = max(lnorm(3),err)

        ! if (err > 4e-1_r8) then
        !   print '(i6, 3es14.4)', i, curvature(i), curvature_ex, err !, c_new_line
        !   call LS_fatal ("large curvature error")
        ! end if
      end if
    end do
    lnorm(1) = lnorm(1) / nvofcell
    lnorm(2) = sqrt(lnorm(2) / nvofcell)

  end function hf_mesh_test

  subroutine compare_vof_fields(vof, vofex)

    real(r8), intent(in) :: vof(:,:), vofex(:,:)

    integer :: i

    do i = 1,size(vof, dim=2)
      if (abs(vof(1,i) - vofex(1,i)) > 1e-2) then
        print *, i, vof(1,i), vofex(1,i)
        call LS_Fatal ("vof initialization error")
      end if
    end do

  end subroutine compare_vof_fields

  subroutine recalculate_mesh_volumes(mesh,gmesh)

    use polyhedron_type
    use hex_types, only: hex_f, hex_e

    type(unstr_mesh), intent(inout) :: mesh
    type(mesh_geom), intent(in) :: gmesh

    integer :: i, ierr
    type(polyhedron) :: cell

    print '(a)', 'recalculating mesh volumes ... '
    do i = 1,mesh%ncell
      call cell%init (ierr, i, mesh, gmesh)
      !call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e)
      mesh%volume(i) = cell%volume()
    end do
    print '(a)', 'done'

  end subroutine recalculate_mesh_volumes

  real(r8) function sinusoid_curvature(x, coeff)
    real(r8), intent(in) :: x(:), coeff(:)
    sinusoid_curvature = -(coeff(2)*coeff(3)**2*(coeff(4)**2*coeff(5)**2*sin(coeff(5&
        &)*x(2))**2 + 1)*cos(coeff(3)*x(1)) + coeff(4)*coeff(5)**2*(coeff(2)**2*coef&
        &f(3)**2*sin(coeff(3)*x(1))**2 + 1)*cos(coeff(5)*x(2)))*(coeff(2)**2*coeff(3&
        &)**2*sin(coeff(3)*x(1))**2 + coeff(4)**2*coeff(5)**2*sin(coeff(5)*x(2))**2 &
        &+ 1)**(-1.5_r8)
  end function sinusoid_curvature

  real(r8) function ellipsoid_curvature_x(dir, x, coeff)

    integer, intent(in) :: dir
    real(r8), intent(in) :: x(:), coeff(:)

    select case(dir)
    case(1)
      ellipsoid_curvature_x = -coeff(1)*x(3)*(-coeff(1)**2*x(3)/(coeff(3)**4*(1 - x(3)&
          &**2/coeff(3)**2 - x(2)**2/coeff(2)**2)) - coeff(1)**2*x(3)**3/(coeff(3)**6*&
          &(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)**2) - coeff(1)**2*x(2)**2*x&
          &(3)/(coeff(2)**4*coeff(3)**2*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2&
          &)**2))/(coeff(3)**2*sqrt(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)*(co&
          &eff(1)**2*x(3)**2/(coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)*&
          &*2)) + coeff(1)**2*x(2)**2/(coeff(2)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/&
          &coeff(2)**2)) + 1)**(1.5_r8)) - coeff(1)/(coeff(3)**2*sqrt(1 - x(3)**2/coeff(3&
          &)**2 - x(2)**2/coeff(2)**2)*sqrt(coeff(1)**2*x(3)**2/(coeff(3)**4*(1 - x(3)&
          &**2/coeff(3)**2 - x(2)**2/coeff(2)**2)) + coeff(1)**2*x(2)**2/(coeff(2)**4*&
          &(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)) + 1)) - coeff(1)*x(3)**2/(&
          &coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)**(1.5_r8)*sqrt(coe&
          &ff(1)**2*x(3)**2/(coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**&
          &2)) + coeff(1)**2*x(2)**2/(coeff(2)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/c&
          &oeff(2)**2)) + 1)) - coeff(1)*x(2)*(-coeff(1)**2*x(2)*x(3)**2/(coeff(2)**2*&
          &coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)**2) - coeff(1)*&
          &*2*x(2)/(coeff(2)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)) - coe&
          &ff(1)**2*x(2)**3/(coeff(2)**6*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**&
          &2)**2))/(coeff(2)**2*sqrt(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)*(c&
          &oeff(1)**2*x(3)**2/(coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)&
          &**2)) + coeff(1)**2*x(2)**2/(coeff(2)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2&
          &/coeff(2)**2)) + 1)**(1.5_r8)) - coeff(1)/(coeff(2)**2*sqrt(1 - x(3)**2/coeff(&
          &3)**2 - x(2)**2/coeff(2)**2)*sqrt(coeff(1)**2*x(3)**2/(coeff(3)**4*(1 - x(3&
          &)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)) + coeff(1)**2*x(2)**2/(coeff(2)**4&
          &*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)) + 1)) - coeff(1)*x(2)**2/&
          &(coeff(2)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)**2)**(1.5_r8)*sqrt(co&
          &eff(1)**2*x(3)**2/(coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/coeff(2)*&
          &*2)) + coeff(1)**2*x(2)**2/(coeff(2)**4*(1 - x(3)**2/coeff(3)**2 - x(2)**2/&
          &coeff(2)**2)) + 1))
    case(2)
      ellipsoid_curvature_x = -coeff(2)*x(3)*(-coeff(2)**2*x(3)/(coeff(3)**4*(1 - x(3)&
          &**2/coeff(3)**2 - x(1)**2/coeff(1)**2)) - coeff(2)**2*x(3)**3/(coeff(3)**6*&
          &(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)**2) - coeff(2)**2*x(1)**2*x&
          &(3)/(coeff(1)**4*coeff(3)**2*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2&
          &)**2))/(coeff(3)**2*sqrt(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)*(co&
          &eff(2)**2*x(3)**2/(coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)*&
          &*2)) + 1 + coeff(2)**2*x(1)**2/(coeff(1)**4*(1 - x(3)**2/coeff(3)**2 - x(1)&
          &**2/coeff(1)**2)))**(1.5_r8)) - coeff(2)/(coeff(3)**2*sqrt(1 - x(3)**2/coeff(3&
          &)**2 - x(1)**2/coeff(1)**2)*sqrt(coeff(2)**2*x(3)**2/(coeff(3)**4*(1 - x(3)&
          &**2/coeff(3)**2 - x(1)**2/coeff(1)**2)) + 1 + coeff(2)**2*x(1)**2/(coeff(1)&
          &**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)))) - coeff(2)*x(3)**2/(&
          &coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)**(1.5_r8)*sqrt(coe&
          &ff(2)**2*x(3)**2/(coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**&
          &2)) + 1 + coeff(2)**2*x(1)**2/(coeff(1)**4*(1 - x(3)**2/coeff(3)**2 - x(1)*&
          &*2/coeff(1)**2)))) - coeff(2)*x(1)*(-coeff(2)**2*x(1)*x(3)**2/(coeff(1)**2*&
          &coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)**2) - coeff(2)*&
          &*2*x(1)/(coeff(1)**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)) - coe&
          &ff(2)**2*x(1)**3/(coeff(1)**6*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**&
          &2)**2))/(coeff(1)**2*sqrt(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)*(c&
          &oeff(2)**2*x(3)**2/(coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)&
          &**2)) + 1 + coeff(2)**2*x(1)**2/(coeff(1)**4*(1 - x(3)**2/coeff(3)**2 - x(1&
          &)**2/coeff(1)**2)))**(1.5_r8)) - coeff(2)/(coeff(1)**2*sqrt(1 - x(3)**2/coeff(&
          &3)**2 - x(1)**2/coeff(1)**2)*sqrt(coeff(2)**2*x(3)**2/(coeff(3)**4*(1 - x(3&
          &)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)) + 1 + coeff(2)**2*x(1)**2/(coeff(1&
          &)**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)))) - coeff(2)*x(1)**2/&
          &(coeff(1)**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)**2)**(1.5_r8)*sqrt(co&
          &eff(2)**2*x(3)**2/(coeff(3)**4*(1 - x(3)**2/coeff(3)**2 - x(1)**2/coeff(1)*&
          &*2)) + 1 + coeff(2)**2*x(1)**2/(coeff(1)**4*(1 - x(3)**2/coeff(3)**2 - x(1)&
          &**2/coeff(1)**2))))
    case(3)
      ellipsoid_curvature_x = -coeff(3)*x(2)*(-coeff(3)**2*x(2)/(coeff(2)**4*(1 - x(2)&
          &**2/coeff(2)**2 - x(1)**2/coeff(1)**2)) - coeff(3)**2*x(2)**3/(coeff(2)**6*&
          &(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)**2) - coeff(3)**2*x(1)**2*x&
          &(2)/(coeff(1)**4*coeff(2)**2*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2&
          &)**2))/(coeff(2)**2*sqrt(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)*(1 &
          &+ coeff(3)**2*x(2)**2/(coeff(2)**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff&
          &(1)**2)) + coeff(3)**2*x(1)**2/(coeff(1)**4*(1 - x(2)**2/coeff(2)**2 - x(1)&
          &**2/coeff(1)**2)))**(1.5_r8)) - coeff(3)/(coeff(2)**2*sqrt(1 - x(2)**2/coeff(2&
          &)**2 - x(1)**2/coeff(1)**2)*sqrt(1 + coeff(3)**2*x(2)**2/(coeff(2)**4*(1 - &
          &x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)) + coeff(3)**2*x(1)**2/(coeff(1)&
          &**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)))) - coeff(3)*x(2)**2/(&
          &coeff(2)**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)**(1.5_r8)*sqrt(1 +&
          & coeff(3)**2*x(2)**2/(coeff(2)**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(&
          &1)**2)) + coeff(3)**2*x(1)**2/(coeff(1)**4*(1 - x(2)**2/coeff(2)**2 - x(1)*&
          &*2/coeff(1)**2)))) - coeff(3)*x(1)*(-coeff(3)**2*x(1)*x(2)**2/(coeff(1)**2*&
          &coeff(2)**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)**2) - coeff(3)*&
          &*2*x(1)/(coeff(1)**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)) - coe&
          &ff(3)**2*x(1)**3/(coeff(1)**6*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**&
          &2)**2))/(coeff(1)**2*sqrt(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)*(1&
          & + coeff(3)**2*x(2)**2/(coeff(2)**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coef&
          &f(1)**2)) + coeff(3)**2*x(1)**2/(coeff(1)**4*(1 - x(2)**2/coeff(2)**2 - x(1&
          &)**2/coeff(1)**2)))**(1.5_r8)) - coeff(3)/(coeff(1)**2*sqrt(1 - x(2)**2/coeff(&
          &2)**2 - x(1)**2/coeff(1)**2)*sqrt(1 + coeff(3)**2*x(2)**2/(coeff(2)**4*(1 -&
          & x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)) + coeff(3)**2*x(1)**2/(coeff(1&
          &)**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)))) - coeff(3)*x(1)**2/&
          &(coeff(1)**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff(1)**2)**(1.5_r8)*sqrt(1 &
          &+ coeff(3)**2*x(2)**2/(coeff(2)**4*(1 - x(2)**2/coeff(2)**2 - x(1)**2/coeff&
          &(1)**2)) + coeff(3)**2*x(1)**2/(coeff(1)**4*(1 - x(2)**2/coeff(2)**2 - x(1)&
          &**2/coeff(1)**2))))
    end select

    ! ellipsoid_curvature_x = -x(3)*sqrt(coeff(1)*coeff(2)*coeff(3))*(-coeff(1)*coeff(&
    !     &2)*x(3)/(coeff(3)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/c&
    !     &oeff(1)**2)) + coeff(1)*coeff(2)*x(3)**3/(coeff(3)**5*(x(3)**2/coeff(3)**2 &
    !     &+ x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)**2) + coeff(1)*x(2)**2*x(3)/(c&
    !     &oeff(2)**3*coeff(3)*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/co&
    !     &eff(1)**2)**2) + coeff(2)*x(1)**2*x(3)/(coeff(1)**3*coeff(3)*(x(3)**2/coeff&
    !     &(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)**2))/(coeff(3)**2*sqrt(&
    !     &x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)*(coeff(1)*&
    !     &coeff(2)*x(3)**2/(coeff(3)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + &
    !     &x(1)**2/coeff(1)**2)) + coeff(1)*coeff(3)*x(2)**2/(coeff(2)**3*(x(3)**2/coe&
    !     &ff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(2)*coeff(3)*&
    !     &x(1)**2/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/c&
    !     &oeff(1)**2)))**(3/2)) - sqrt(coeff(1)*coeff(2)*coeff(3))/(coeff(3)**2*sqrt(&
    !     &x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)*sqrt(coeff&
    !     &(1)*coeff(2)*x(3)**2/(coeff(3)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**&
    !     &2 + x(1)**2/coeff(1)**2)) + coeff(1)*coeff(3)*x(2)**2/(coeff(2)**3*(x(3)**2&
    !     &/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(2)*coeff&
    !     &(3)*x(1)**2/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)*&
    !     &*2/coeff(1)**2)))) + x(3)**2*sqrt(coeff(1)*coeff(2)*coeff(3))/(coeff(3)**4*&
    !     &(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)**(3/2)*sq&
    !     &rt(coeff(1)*coeff(2)*x(3)**2/(coeff(3)**3*(x(3)**2/coeff(3)**2 + x(2)**2/co&
    !     &eff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(1)*coeff(3)*x(2)**2/(coeff(2)**3*&
    !     &(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(&
    !     &2)*coeff(3)*x(1)**2/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2&
    !     & + x(1)**2/coeff(1)**2)))) - x(2)*sqrt(coeff(1)*coeff(2)*coeff(3))*(coeff(1&
    !     &)*x(2)*x(3)**2/(coeff(2)*coeff(3)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2&
    !     &)**2 + x(1)**2/coeff(1)**2)**2) - coeff(1)*coeff(3)*x(2)/(coeff(2)**3*(x(3)&
    !     &**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(1)*co&
    !     &eff(3)*x(2)**3/(coeff(2)**5*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(&
    !     &1)**2/coeff(1)**2)**2) + coeff(3)*x(1)**2*x(2)/(coeff(1)**3*coeff(2)*(x(3)*&
    !     &*2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)**2))/(coeff(2)*&
    !     &*2*sqrt(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)*(c&
    !     &oeff(1)*coeff(2)*x(3)**2/(coeff(3)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(&
    !     &2)**2 + x(1)**2/coeff(1)**2)) + coeff(1)*coeff(3)*x(2)**2/(coeff(2)**3*(x(3&
    !     &)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(2)*c&
    !     &oeff(3)*x(1)**2/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x&
    !     &(1)**2/coeff(1)**2)))**(3/2)) - sqrt(coeff(1)*coeff(2)*coeff(3))/(coeff(2)*&
    !     &*2*sqrt(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)*sq&
    !     &rt(coeff(1)*coeff(2)*x(3)**2/(coeff(3)**3*(x(3)**2/coeff(3)**2 + x(2)**2/co&
    !     &eff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(1)*coeff(3)*x(2)**2/(coeff(2)**3*&
    !     &(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(&
    !     &2)*coeff(3)*x(1)**2/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2&
    !     & + x(1)**2/coeff(1)**2)))) + x(2)**2*sqrt(coeff(1)*coeff(2)*coeff(3))/(coef&
    !     &f(2)**4*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)**&
    !     &(3/2)*sqrt(coeff(1)*coeff(2)*x(3)**2/(coeff(3)**3*(x(3)**2/coeff(3)**2 + x(&
    !     &2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(1)*coeff(3)*x(2)**2/(coef&
    !     &f(2)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) &
    !     &+ coeff(2)*coeff(3)*x(1)**2/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coe&
    !     &ff(2)**2 + x(1)**2/coeff(1)**2)))) - x(1)*sqrt(coeff(1)*coeff(2)*coeff(3))*&
    !     &(coeff(2)*x(1)*x(3)**2/(coeff(1)*coeff(3)**3*(x(3)**2/coeff(3)**2 + x(2)**2&
    !     &/coeff(2)**2 + x(1)**2/coeff(1)**2)**2) + coeff(3)*x(1)*x(2)**2/(coeff(1)*c&
    !     &oeff(2)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2&
    !     &)**2) - coeff(2)*coeff(3)*x(1)/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2)**2/&
    !     &coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(2)*coeff(3)*x(1)**3/(coeff(1)**&
    !     &5*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)**2))/(c&
    !     &oeff(1)**2*sqrt(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1&
    !     &)**2)*(coeff(1)*coeff(2)*x(3)**2/(coeff(3)**3*(x(3)**2/coeff(3)**2 + x(2)**&
    !     &2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(1)*coeff(3)*x(2)**2/(coeff(2)&
    !     &**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + co&
    !     &eff(2)*coeff(3)*x(1)**2/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2&
    !     &)**2 + x(1)**2/coeff(1)**2)))**(3/2)) - sqrt(coeff(1)*coeff(2)*coeff(3))/(c&
    !     &oeff(1)**2*sqrt(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1&
    !     &)**2)*sqrt(coeff(1)*coeff(2)*x(3)**2/(coeff(3)**3*(x(3)**2/coeff(3)**2 + x(&
    !     &2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(1)*coeff(3)*x(2)**2/(coef&
    !     &f(2)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) &
    !     &+ coeff(2)*coeff(3)*x(1)**2/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coe&
    !     &ff(2)**2 + x(1)**2/coeff(1)**2)))) + x(1)**2*sqrt(coeff(1)*coeff(2)*coeff(3&
    !     &))/(coeff(1)**4*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(&
    !     &1)**2)**(3/2)*sqrt(coeff(1)*coeff(2)*x(3)**2/(coeff(3)**3*(x(3)**2/coeff(3)&
    !     &**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(1)**2)) + coeff(1)*coeff(3)*x(2)*&
    !     &*2/(coeff(2)**3*(x(3)**2/coeff(3)**2 + x(2)**2/coeff(2)**2 + x(1)**2/coeff(&
    !     &1)**2)) + coeff(2)*coeff(3)*x(1)**2/(coeff(1)**3*(x(3)**2/coeff(3)**2 + x(2&
    !     &)**2/coeff(2)**2 + x(1)**2/coeff(1)**2))))

    ellipsoid_curvature_x = abs(ellipsoid_curvature_x)

  end function ellipsoid_curvature_x

  subroutine compute_exact_curvature_field(curvature_ex, interface_centroid, int_norm, mesh, &
      gmesh, mixed_cells)

    use mesh_subset_type
    use mixed_cell_subset_constructor

    real(r8), allocatable, intent(out) :: curvature_ex(:)
    real(r8), intent(in) :: interface_centroid(:,:), int_norm(:,:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    type(mesh_subset), intent(in) :: mixed_cells

    integer :: c, i

    if (allocated(curvature_ex)) deallocate(curvature_ex)
    allocate(curvature_ex(mesh%ncell))
    curvature_ex = 0

    do c = 1,mixed_cells%ncell
      i = mixed_cells%cell_id(c)

      !curvature_ex(i) = 1 / 0.35_r8 ! cylinder
      !curvature_ex(i) = 2 / 0.35_r8 ! sphere
      curvature_ex(i) = ellipsoid_curvature_x(maxloc(abs(int_norm(:,1,i)),1), &
          !interface_centroid(:,c), &
          gmesh%xc(:,i), &
          [0.35_r8, 0.3_r8, 0.2_r8])
    end do

  end subroutine compute_exact_curvature_field

end module analytic_surface_type_test
