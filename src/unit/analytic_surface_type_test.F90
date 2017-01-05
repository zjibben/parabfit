module analytic_surface_type_test

  use kinds, only: r8
  use analytic_surface_type
  use logging_services
  implicit none
  private

  public :: analytic_surface_test_suite

contains

  subroutine analytic_surface_test_suite ()

    print '(a)'
    print '(a)', 'ANALYTIC_SURFACE'
    print '(a)', '===================================================='

    call plane_test ()
    call parabola_test ()
    call mesh_2d_test (16)

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

  subroutine mesh_2d_test (mesh_size)

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

    character(:), allocatable :: errmsg
    type(unstr_mesh) :: mesh
    type(mesh_geom) :: gmesh
    type(surface) :: intrec
    type(parameter_list), pointer :: plist
    real(r8) :: curvature_exact, lnorm(3)
    real(r8) :: vof(2,mesh_size*mesh_size), curvature(mesh_size*mesh_size)
    real(r8), allocatable :: int_norm(:,:,:)
    type(multimat_cell) :: cell
    integer :: i, infile, ierr, nvofcell

    ! create a regular 2D mesh
    mesh = new_unstr_mesh ([-0.5_r8, -0.5_r8, 0.0_r8], [0.5_r8, 0.5_r8, 1.0_r8], &
        [mesh_size,mesh_size,1])
    call gmesh%init (mesh)

    ! fill the mesh with volume fractions for a circle
    ! right now this relies on an input file to describe the cylinder. In the future I'd like to
    ! initialize material_geometry_types without a parameter_list_type, or initialize a
    ! parameter_list_type without a JSON file
    open(newunit=infile,file='cylinder.json',action='read',access='stream')
    call parameter_list_from_json_stream (infile, plist, errmsg)
    if (.not.associated(plist)) call LS_fatal ("error reading input file:" // C_NEW_LINE // errmsg)
    close(infile)

    plist => plist%sublist('initial-vof')
    call vof_initialize (mesh, plist, vof, [1,2], 2)
    curvature_exact = 4.0_r8
    deallocate(plist)
    
    ! get the interface reconstructions
    int_norm = interface_normal (vof, mesh, gmesh, .false.)
    do i = 1,mesh%ncell
      call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, mesh%volume(i), &
          gmesh%outnorm(:,:,i))
      if (ierr /= 0) call LS_fatal ('cell_outward_volflux failed: could not initialize cell')

      !int_norm(1:2,1,i) = normalize(gmesh%xc(1:2,i))
      call cell%partition (vof(:,i), int_norm(:,:,i))

      call intrec%append (cell%interface_polygon(1), i)
    end do

    ! get the curvature
    curvature = 0.0_r8; lnorm = 0.0_r8; nvofcell = 0
    do i = 1,mesh%ncell
      ! TODO: this really should be in any cell neighboring a cell containing the interface
      if (vof(1,i) > cutvof .and. vof(1,i) < 1.0_r8-cutvof) then
        curvature(i) = curvature_from_patch (intrec%local_patch(i,gmesh))
        nvofcell = nvofcell + 1
        lnorm(1) = lnorm(1) + abs(curvature(i) - curvature_exact)
        print *, i, curvature(i)
      end if
    end do
    lnorm(1) = lnorm(1) / real(nvofcell,r8)

    print '(a,es15.4)', 'Finished. L1 = ',lnorm(1)

  end subroutine mesh_2d_test
    
end module analytic_surface_type_test
