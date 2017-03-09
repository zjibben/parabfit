module analytic_surface_type_test

  use kinds, only: r8
  use analytic_surface_type
  use paraboloid_type
  use logging_services
  implicit none
  private

  public :: analytic_surface_test_suite

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
    ! print *

    !lnorm = mesh_2d_test (20, 'cylinder.json')
    !lnorm = mesh_2d_test (80, 'cylinder.json')
    !lnorm = mesh_2d_test (100, 'cylinder.json')

    open (98, file="ftstr_conv_dump.txt")
    write (98, '(a)') '# dx l1 l2 l3'
    do i = 1,30
      ncell = floor(10.0_r8 * 1.15_r8**i)
      lnorm = mesh_2d_test (ncell, 'cylinder.json')
      write (98, '(4es15.5)') 1.0_r8 / real(ncell, r8), lnorm
    end do
    close(98)

    print '(a)', '===================================================='
    print '(a)'

  end subroutine analytic_surface_test_suite

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
        x(1,ind) = real(i-N/2+1,r8)*dx
        x(2,ind) = real(j-N/2+1,r8)*dx
        x(3,ind) = x(1,ind) + x(2,ind)
      end do
    end do
    
    call surf%bestFit (x)

    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature([0.0_r8,0.0_r8,0.0_r8])
    
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
        x(1,ind) = real(i-N/2+1,r8)*dx
        x(2,ind) = real(j-N/2+1,r8)*dx
        x(3,ind) = x(1,ind)**2 + x(2,ind)**2 ! parabola
      end do
    end do
    
    call surf%bestFit (x)

    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature([0.0_r8,0.0_r8,0.0_r8])

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

    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature(sum(x(:,1:3), dim=2) / 3.0_r8)
    
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

    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature(sum(x(:,1:3), dim=2) / 3.0_r8)
    
  end subroutine messy_test2

  subroutine messy_test3 ()
    
    real(r8), allocatable :: x(:,:)
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

    call surf%bestFit (x)
    !call surf%bestParaboloidFit (x)
    !call surf%bestOneSheetFit (x)

    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature(sum(x(:,1:3), dim=2) / 3.0_r8)
    
  end subroutine messy_test3

  subroutine messy_test4 ()
    
    real(r8), allocatable :: x(:,:)
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

    call surf%bestFit (x)
    !call surf%bestParaboloidFit (x)
    !call surf%bestOneSheetFit (x)

    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature(sum(x(:,1:3), dim=2) / 3.0_r8)
    
  end subroutine messy_test4

  function mesh_2d_test (mesh_size, shape_filename) result(lnorm)

    use,intrinsic :: iso_c_binding, only: C_NEW_LINE
    use consts, only: cutvof
    use unstr_mesh_type
    use mesh_geom_type
    use unstr_mesh_factory
    use surface_type
    use parameter_list_type
    use parameter_list_json
    use int_norm_module
    use interface_patch_type
    use vof_init
    use multimat_cell_type
    use hex_types, only: hex_f, hex_e
    use array_utils, only: normalize, isZero
    use curvature_hf

    integer, intent(in) :: mesh_size
    character(*), intent(in) :: shape_filename

    character(:), allocatable :: errmsg
    type(unstr_mesh) :: mesh
    type(mesh_geom) :: gmesh
    type(surface) :: intrec
    type(parameter_list), pointer :: plist
    real(r8) :: curvature_exact, lnorm(3), err, &
        !vof(2,mesh_size**3), curvature(mesh_size**3)
        vof(2,mesh_size*mesh_size*3), curvature(mesh_size*mesh_size*3)
    real(r8), allocatable :: int_norm(:,:,:)
    type(multimat_cell) :: cell
    integer :: i, infile, ierr, nvofcell

    ! create a regular 2D mesh
    mesh = new_unstr_mesh ([-0.5_r8, -0.5_r8, -3.0_r8*0.5_r8/real(mesh_size,r8)], &
        [0.5_r8, 0.5_r8, 3.0_r8*0.5_r8/real(mesh_size,r8)], [mesh_size,mesh_size,3])
    ! mesh = new_unstr_mesh ([-0.5_r8, -0.5_r8, -0.5_r8], &
    !     [0.5_r8, 0.5_r8, 0.5_r8], [mesh_size,mesh_size,mesh_size])
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
    call vof_initialize (mesh, plist, vof, [1,2], 2)
    !curvature_exact = 1.0_r8/0.35_r8 + 1.0_r8/0.35_r8 ! sphere
    curvature_exact = 1.0_r8/0.25_r8 + 0.0_r8 ! cylinder
    deallocate(plist)
    
    ! get the interface reconstructions
    int_norm = interface_normal (vof, mesh, gmesh, .false.)
    do i = 1,mesh%ncell
      call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, mesh%volume(i), &
          gmesh%outnorm(:,:,i))
      if (ierr /= 0) call LS_fatal ('cell_outward_volflux failed: could not initialize cell')
      
      int_norm(:,:,i) = 0.0_r8
      int_norm(1:2,1,i) = -normalize(gmesh%xc(1:2,i))
      int_norm(:,2,i) = -int_norm(:,1,i)
      call cell%partition (vof(:,i), int_norm(:,:,i))

      call intrec%append (cell%interface_polygon(1), i)
    end do
    call intrec%write_ply ('cylsurf.ply')

    ! get the curvature
    curvature = 0.0_r8; lnorm = 0.0_r8; nvofcell = 0
    !curvature = abs(curvatureHF(vof(1,:), int_norm(:,1,:), mesh, gmesh))
    do i = 1,mesh%ncell
      ! TODO: need to handle boundaries. this might be automatic.
      if (any(gmesh%cneighbor(:,i)<1)) cycle

      ! TODO: this really should be in any cell neighboring a cell containing the interface
      if (vof(1,i) > cutvof .and. vof(1,i) < 1.0_r8-cutvof) then ! .and. .not.isZero(curvature(i))) then
      !if (vof(1,i) > 5e-2_r8 .and. vof(1,i) < 1.0_r8-5e-2_r8) then! .and. .not.isZero(curvature(i))) then
        curvature(i) = abs(curvature_from_patch (intrec%local_patch(i,gmesh)))

        ! append to norms
        err = abs((curvature(i) - curvature_exact) / curvature_exact)
        nvofcell = nvofcell + 1
        lnorm(1) = lnorm(1) + err
        lnorm(2) = lnorm(2) + err**2
        lnorm(3) = max(lnorm(3),err)
        
        ! if (err > 8.6e-2_r8) then
        !   print '(i6, 3es14.4)', i, curvature(i), curvature_exact, err !, c_new_line
        !   print '(2es15.5)', vof(1,i), 1.0_r8 - vof(1,i)
        !   call LS_fatal ("large curvature error")
        ! end if
      end if
    end do
    lnorm(1) = lnorm(1) / real(nvofcell,r8)
    lnorm(2) = sqrt(lnorm(2) / real(nvofcell,r8))

    print '(i5, a,3es15.4)', mesh_size, '  Finished. L1,L2,Linf = ',lnorm

  end function mesh_2d_test
    
end module analytic_surface_type_test
