module interface_point_test

  use kinds, only: r8
  use unstr_mesh_type
  use mesh_geom_type
  use logging_services
  implicit none
  private

  public :: interface_point_test_suite

contains

  ! note: not sure that this should actually converge any better than 1st order.
  !       is the centroid of the plane really *supposed* to be on the interface?
  !       i think probably not. instead the plane should represent the interface
  !       sort of on average over the cell.
  subroutine interface_point_test_suite ()

    integer :: i, ncell
    real(r8) :: lnormYG(3), lnormHF(3)

    print '(a)'
    print '(a)', 'INTERFACE POINTS'
    print '(a)', '===================================================='

    open (98, file="ft_conv_dump.txt")
    open (99, file="hf_conv_dump.txt")

    write (98, '(a)') '# dx l1 l2 l3'
    write (99, '(a)') '# dx l1 l2 l3'

    do i = 1,4
    !do i = 2,2
      ncell = 10 * 2**i
    ! do i = 1,25
    !   ncell = floor(10 * 1.15_r8**i)
      call mesh_2d_test (ncell, 'cylinder.json', lnormYG, lnormHF)
      write (98, '(4es15.5)') 1.0_r8 / ncell, lnormYG
      write (99, '(4es15.5)') 1.0_r8 / ncell, lnormHF
    end do

    close(98); close(99)

    print '(a)', '===================================================='
    print '(a)'

  end subroutine interface_point_test_suite

  subroutine mesh_2d_test (mesh_size, shape_filename, lnormYG, lnormHF)

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
    use hex_types, only: hex_f, hex_e
    use array_utils, only: normalize, isZero
    use curvature_hf

    integer, intent(in) :: mesh_size
    character(*), intent(in) :: shape_filename
    real(r8), intent(out) :: lnormYG(:), lnormHF(:)

    character(:), allocatable :: errmsg
    type(unstr_mesh) :: mesh
    type(mesh_geom) :: gmesh
    type(parameter_list), pointer :: plist
    real(r8) :: &
        !vof(2,mesh_size**3), curvature(mesh_size**3)
        vof(2,mesh_size*mesh_size*3), curvature(mesh_size*mesh_size*3)
    real(r8), allocatable :: int_norm(:,:,:)
    integer :: infile

    ! create a regular 2D mesh
    mesh = new_unstr_mesh ([-0.5_r8, -0.5_r8, -3*0.5_r8/mesh_size], &
        [0.5_r8, 0.5_r8, 3*0.5_r8 / mesh_size], [mesh_size,mesh_size,3])
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
    call vof_initialize (mesh, gmesh, plist, vof, [1,2], 2)
    deallocate(plist)

    ! get the interface reconstructions
    int_norm = interface_normal (vof, mesh, gmesh, .false.)

    ! calculate errors for YG and HF curvature methods
    lnormYG = young_mesh_test (vof, int_norm, mesh, gmesh)
    lnormHF = hf_mesh_test (vof, int_norm, mesh, gmesh)

    print '(i5, 2(a,3es10.2))', mesh_size, '  YG L1,L2,Linf = ',lnormYG, &
        ',  HF L1,L2,Linf = ',lnormHF

  end subroutine mesh_2d_test

  function young_mesh_test (vof, int_norm, mesh, gmesh) result(lnorm)

    !use consts, only: cutvof
    use hex_types, only: hex_f, hex_e
    use multimat_cell_type
    use polygon_type
    use ieee_arithmetic, only: ieee_is_nan

    real(r8), intent(in) :: vof(:,:), int_norm(:,:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8) :: lnorm(3)

    integer :: i, nvofcell, imax, ierr
    real(r8) :: err
    type(multimat_cell) :: cell
    type(polygon) :: interface_polygon

    lnorm = 0; nvofcell = 0
    do i = 1,mesh%ncell
      if (any(gmesh%cneighbor(:,i)<1)) cycle ! WARN: skipping boundaries. BCs might be automatic?

      ! TODO: this really should be in any cell neighboring a cell containing the interface
      if (vof(1,i) < 1e-2_r8 .or. vof(1,i) > 1-1e-2_r8) cycle

      call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, gmesh%outnorm(:,:,i), &
          mesh%volume(i))
      if (ierr /= 0) call LS_fatal ('young_mesh_test: could not initialize cell')
      call cell%partition (vof(:,i), int_norm(:,:,i))
      interface_polygon = cell%interface_polygon(1)
      if (interface_polygon%nverts < 3) cycle
      err = abs(norm2(interface_polygon%centroid()) - 0.25_r8)

      if (err > lnorm(3)) imax = i
      nvofcell = nvofcell + 1
      lnorm(1) = lnorm(1) + err
      lnorm(2) = lnorm(2) + err**2
      lnorm(3) = max(lnorm(3),err)

      ! if (err > 3e-1_r8) then
      !   print '(i6, 3es14.4)', i, err
      !   print '(3es16.3)', int_norm(:,i)
      !   print '(3es16.3)', norm_ex
      !   call LS_fatal ("large normal error")
      ! end if
    end do
    lnorm(1) = lnorm(1) / nvofcell
    lnorm(2) = sqrt(lnorm(2) / nvofcell)

  end function young_mesh_test

  function hf_mesh_test (vof, int_norm, mesh, gmesh) result(lnorm)

    use consts, only: cutvof
    use hex_types, only: hex_f, hex_e
    use curvature_hf
    use multimat_cell_type
    use polygon_type

    real(r8), intent(in) :: vof(:,:), int_norm(:,:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8) :: lnorm(3)

    integer :: i, nvofcell, ierr, imax
    real(r8) :: err
    real(r8), allocatable :: throwaway(:), int_norm_hf(:,:), int_norm_hf2(:,:,:)
    type(multimat_cell) :: cell
    type(polygon) :: interface_polygon

    call heightFunction(throwaway, int_norm_hf, vof(1,:), int_norm(:,1,:), mesh, gmesh)

    allocate(int_norm_hf2(3,2,mesh%ncell))
    int_norm_hf2(:,1,:) =  int_norm_hf
    int_norm_hf2(:,2,:) = -int_norm_hf

    lnorm = 0; nvofcell = 0
    do i = 1,mesh%ncell
      if (any(gmesh%cneighbor(:,i)<1)) cycle ! WARN: skipping boundaries. BCs might be automatic?

      ! TODO: this really should be in any cell neighboring a cell containing the interface
      if (vof(1,i) < 1e-2_r8 .or. vof(1,i) > 1-1e-2_r8) cycle

      call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, gmesh%outnorm(:,:,i), &
          mesh%volume(i))
      if (ierr /= 0) call LS_fatal ('young_mesh_test: could not initialize cell')
      call cell%partition (vof(:,i), int_norm_hf2(:,:,i))
      interface_polygon = cell%interface_polygon(1)
      err = abs(norm2(interface_polygon%centroid()) - 0.25_r8)

      if (err > lnorm(3)) imax = i
      nvofcell = nvofcell + 1
      lnorm(1) = lnorm(1) + err
      lnorm(2) = lnorm(2) + err**2
      lnorm(3) = max(lnorm(3),err)
    end do
    lnorm(1) = lnorm(1) / nvofcell
    lnorm(2) = sqrt(lnorm(2) / nvofcell)

  end function hf_mesh_test

end module interface_point_test
