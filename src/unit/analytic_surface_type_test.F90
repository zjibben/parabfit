module analytic_surface_type_test

  use kinds, only: r8
  use analytic_surface_type
  use logging_services
  implicit none
  private

  public :: analytic_surface_test_suite

contains

  subroutine analytic_surface_test_suite ()

    integer :: i

    print '(a)'
    print '(a)', 'ANALYTIC_SURFACE'
    print '(a)', '===================================================='

    ! call plane_test ()
    ! call parabola_test ()
    ! call messy_test ()

    !call mesh_2d_test (2**5)
    do i = 2,6
      call mesh_2d_test (2**i, 'cylinder.json')
    end do

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
    
    call surf%init (x)

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
    
    call surf%init (x)

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

    call surf%init (x)

    print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature([0.0_r8,0.0_r8,0.0_r8])
    
  end subroutine messy_test

  subroutine mesh_2d_test (mesh_size, shape_filename)

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
    use array_utils, only: normalize

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
    !curvature_exact = 0.5_r8 * (1.0_r8/0.35_r8 + 1.0_r8/0.35_r8) ! sphere
    curvature_exact = 0.5_r8 * (1.0_r8/0.35_r8 + 0.0_r8) ! cylinder
    deallocate(plist)
    
    ! get the interface reconstructions
    int_norm = interface_normal (vof, mesh, gmesh, .false.)
    do i = 1,mesh%ncell
      call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, mesh%volume(i), &
          gmesh%outnorm(:,:,i))
      if (ierr /= 0) call LS_fatal ('cell_outward_volflux failed: could not initialize cell')

      ! int_norm(:,:,i) = 0.0_r8
      ! int_norm(1:2,1,i) = -normalize(gmesh%xc(1:2,i))
      call cell%partition (vof(:,i), int_norm(:,:,i))

      call intrec%append (cell%interface_polygon(1), i)
    end do
    call intrec%write_ply ('cylsurf.ply')

    ! get the curvature
    curvature = 0.0_r8; lnorm = 0.0_r8; nvofcell = 0
    do i = 1,mesh%ncell
      ! TODO: need to handle boundaries. this might be automatic.
      if (any(gmesh%fneighbor(:,i)<1)) cycle

      ! TODO: this really should be in any cell neighboring a cell containing the interface
      if (vof(1,i) > cutvof .and. vof(1,i) < 1.0_r8-cutvof) then
        curvature(i) = abs(curvature_from_patch (intrec%local_patch(i,gmesh)))

        ! append to norms
        err = abs(curvature(i) - curvature_exact)
        nvofcell = nvofcell + 1
        lnorm(1) = lnorm(1) + err
        lnorm(2) = lnorm(2) + err**2
        lnorm(3) = max(lnorm(3),err)
        
        ! print '(i6, 3es15.4)', i, curvature(i), curvature_exact, err
        ! if (err > 1e3) stop
      end if
    end do
    lnorm(1) = lnorm(1) / real(nvofcell,r8)
    lnorm(2) = sqrt(lnorm(2) / real(nvofcell,r8))

    print '(a,3es15.4)', 'Finished. L1,L2,Linf = ',lnorm

  end subroutine mesh_2d_test
    
end module analytic_surface_type_test
