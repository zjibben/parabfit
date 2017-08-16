module normal_vector_test

  use kinds, only: r8
  use unstr_mesh_type
  use mesh_geom_type
  use logging_services
  implicit none
  private

  public :: normal_vector_test_suite

contains

  subroutine normal_vector_test_suite ()

    integer :: i, ncell, fh1, fh2, fh3
    real(r8) :: lnormYG(3), lnormHF(3), lnormLV(3)

    print '(a)'
    print '(a)', 'NORMALS'
    print '(a)', '===================================================='

    open (newunit=fh1, file="yg_norm_conv_dump.txt")
    open (newunit=fh2, file="hf_norm_conv_dump.txt")
    open (newunit=fh3, file="hf_norm_conv_dump.txt")

    write (fh1, '(a)') '# dx l1 l2 l3'
    write (fh2, '(a)') '# dx l1 l2 l3'
    write (fh3, '(a)') '# dx l1 l2 l3'

    do i = 2,2 !4
      ncell = 10 * 2**i
    ! do i = 1,25
    !   ncell = floor(10 * 1.15_r8**i)
      call mesh_2d_test (ncell, 'plane.json', lnormYG, lnormHF, lnormLV)
      write (fh1, '(4es15.5)') 1.0_r8 / ncell, lnormYG
      write (fh2, '(4es15.5)') 1.0_r8 / ncell, lnormHF
      write (fh3, '(4es15.5)') 1.0_r8 / ncell, lnormLV
    end do

    close(fh1); close(fh2); close(fh3)

    print '(a)', '===================================================='
    print '(a)'

  end subroutine normal_vector_test_suite

  ! subroutine plane_test ()

  !   real(r8), allocatable :: x(:,:)
  !   real(r8) :: dx
  !   integer :: N,ind,i,j
  !   type(analytic_surface) :: surf

  !   dx = 0.1_r8
  !   N = 9

  !   allocate(x(3,N*N))
  !   do i = 1,N
  !     do j = 1,N
  !       ind = i + (j-1)*N
  !       x(1,ind) = (i-N/2+1)*dx
  !       x(2,ind) = (j-N/2+1)*dx
  !       x(3,ind) = x(1,ind) + x(2,ind)
  !     end do
  !   end do

  !   call surf%bestFit (x)

  !   print '(dt,a,es12.4)', surf, ',     curvature: ', surf%curvature([0.0_r8,0.0_r8,0.0_r8])

  ! end subroutine plane_test

  subroutine mesh_2d_test (mesh_size, shape_filename, lnormYG, lnormHF, lnormLV)

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
    use vof_io

    integer, intent(in) :: mesh_size
    character(*), intent(in) :: shape_filename
    real(r8), intent(out) :: lnormYG(:), lnormHF(:), lnormLV(:)

    character(:), allocatable :: errmsg, filename
    character(30) :: tmp
    type(unstr_mesh) :: mesh
    type(mesh_geom) :: gmesh
    type(parameter_list), pointer :: plist
    real(r8) :: norm_ex(3), &
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

    write (tmp, '(a,i0,a)') "vof_plane_2d_", mesh_size, ".dat"
    filename = trim(adjustl(tmp))
    if (file_exists(filename)) then
      call read_vof_field(filename, vof)
    else
      call vof_initialize (mesh, plist, vof, [1,2], 2)
      call store_vof_field(filename, vof)
    end if
    deallocate(plist)

    norm_ex = -[0.8_r8, 0.6_r8, 0.0_r8]

    ! get the interface reconstructions
    int_norm = interface_normal (vof, mesh, gmesh, .false.)

    ! calculate errors for YG and HF curvature methods
    lnormYG = young_mesh_test (vof(1,:), int_norm(:,1,:), mesh, gmesh, norm_ex)
    lnormHF = hf_mesh_test (vof(1,:), int_norm(:,1,:), mesh, gmesh, norm_ex)
    lnormLV = lvira_mesh_test (vof, int_norm(:,1,:), mesh, gmesh, norm_ex)

    print '(i5, 3(a,3es10.2))', mesh_size, '  YG L1,L2,Linf = ',lnormYG, &
        ',  HF L1,L2,Linf = ',lnormHF, &
        ',  LV L1,L2,Linf = ',lnormLV

  end subroutine mesh_2d_test

  function young_mesh_test (vof, int_norm, mesh, gmesh, norm_ex) result(lnorm)

    !use consts, only: cutvof

    real(r8), intent(in) :: vof(:), int_norm(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8), intent(in) :: norm_ex(:)
    real(r8) :: lnorm(3)

    integer :: i, nvofcell, imax
    real(r8) :: err

    lnorm = 0; nvofcell = 0
    do i = 1,mesh%ncell
      if (any(gmesh%cneighbor(:,i)<1)) cycle ! WARN: skipping boundaries. BCs might be automatic?

      ! TODO: this really should be in any cell neighboring a cell containing the interface
      if (vof(i) < 1e-2_r8 .or. vof(i) > 1-1e-2_r8) cycle

      err = norm2(int_norm(:,i) - norm_ex)

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

  function hf_mesh_test (vof, int_norm, mesh, gmesh, norm_ex) result(lnorm)

    use consts, only: cutvof
    use curvature_hf

    real(r8), intent(in) :: vof(:), int_norm(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8), intent(in) :: norm_ex(3)
    real(r8) :: lnorm(3)

    integer :: i, nvofcell, ierr, imax
    real(r8) :: err
    real(r8), allocatable :: throwaway(:), int_norm_hf(:,:)

    call heightFunction(throwaway, int_norm_hf, vof, int_norm, mesh, gmesh)

    lnorm = 0; nvofcell = 0
    do i = 1,mesh%ncell
      if (any(gmesh%cneighbor(:,i)<1)) cycle ! WARN: skipping boundaries. BCs might be automatic?

      ! TODO: this really should be in any cell neighboring a cell containing the interface
      if (vof(i) < 1e-2_r8 .or. vof(i) > 1-1e-2_r8) cycle

      err = norm2(int_norm_hf(:,i) - norm_ex)

      if (err > lnorm(3)) imax = i
      nvofcell = nvofcell + 1
      lnorm(1) = lnorm(1) + err
      lnorm(2) = lnorm(2) + err**2
      lnorm(3) = max(lnorm(3),err)

      ! if (err > 1e-5_r8) then
      !   print '(i6, 3es14.4)', i, err
      !   print '(3es16.5)', int_norm_hf(:,i)
      !   print '(3es16.5)', norm_ex
      !   print '(3es16.5)', int_norm_hf(:,i) - norm_ex
      !   print '(3es16.5)', norm2(int_norm_hf(:,i) - norm_ex)
      !   call LS_fatal ("large normal error")
      ! end if
    end do
    lnorm(1) = lnorm(1) / nvofcell
    lnorm(2) = sqrt(lnorm(2) / nvofcell)

  end function hf_mesh_test

  function lvira_mesh_test (vof, int_norm, mesh, gmesh, norm_ex) result(lnorm)

    use consts, only: cutvof
    use lvira_normals

    real(r8), intent(in) :: vof(:,:), int_norm(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8), intent(in) :: norm_ex(3)
    real(r8) :: lnorm(3)

    integer :: i, nvofcell, ierr, imax
    real(r8) :: err
    real(r8), allocatable :: int_norm_lvira(:,:,:)

    call interface_normals_lvira(int_norm_lvira, vof, mesh, gmesh)

    lnorm = 0; nvofcell = 0
    do i = 1,mesh%ncell
      if (any(gmesh%cneighbor(:,i)<1)) cycle ! WARN: skipping boundaries. BCs might be automatic?

      ! TODO: this really should be in any cell neighboring a cell containing the interface
      if (vof(1,i) < 1e-2_r8 .or. vof(1,i) > 1-1e-2_r8) cycle

      err = norm2(int_norm_lvira(:,1,i) - norm_ex)

      if (err > lnorm(3)) imax = i
      nvofcell = nvofcell + 1
      lnorm(1) = lnorm(1) + err
      lnorm(2) = lnorm(2) + err**2
      lnorm(3) = max(lnorm(3),err)

      if (err > 1.7e-5_r8) then
        print '(i6, 3es14.4)', i, err
        print '(3es16.5)', int_norm_lvira(:,1,i)
        print '(3es16.5)', norm_ex
        print '(3es16.5)', int_norm_lvira(:,1,i) - norm_ex
        print '(3es16.5)', norm2(int_norm_lvira(:,1,i) - norm_ex)
        call LS_fatal ("large normal error")
      end if
    end do
    lnorm(1) = lnorm(1) / nvofcell
    lnorm(2) = sqrt(lnorm(2) / nvofcell)

  end function lvira_mesh_test

end module normal_vector_test
