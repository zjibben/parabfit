module curvature_tests

  use kinds, only: r8
  use logging_services
  use unstr_mesh_type
  use mesh_geom_type
  implicit none
  private

  public :: curvature_grid_refinement_study

  real(r8), parameter :: pi = 4*atan(1.0_r8)

  integer, parameter :: &
      SPHERE = 1, &
      ELLIPSOID = 2, &
      SINUSOID = 3

  integer, parameter :: &
      REG = 1, &
      RND = 2, &
      TET = 3

  type :: curvature_test
    integer :: mesh_type, shape_type
    character(:), allocatable :: mesh_name, shape_name, shape_filename
    real(r8) :: error_cutoff
  contains
    procedure :: init => init_curvature_test
  end type curvature_test

contains

  subroutine curvature_grid_refinement_study (mesh_type, shape_type, error_cutoff)

    integer :: mesh_type, shape_type
    real(r8), intent(in) :: error_cutoff

    type(curvature_test) :: test_params
    integer :: i, ncell, fh_fit, fh_hf
    real(r8) :: lnormFT(3), lnormHF(3)
    character(256) :: filename

    call test_params%init(mesh_type, shape_type, error_cutoff)

    print '(a)'
    print '(3a)', test_params%shape_name, ' ', test_params%mesh_name
    print '(a)', '===================================================='

    ! TODO: create conv directory
    write(filename, '(5a)') "conv/ft_conv_", test_params%shape_name, &
        "_", test_params%mesh_name, ".txt"
    open(newunit=fh_fit, file=trim(adjustl(filename)))
    write(fh_fit, '(a)') '# ncell l1 l2 linf'

    fh_hf = 0
    if (test_params%mesh_type == REG) then
      write(filename, '(5a)') "conv/hf_conv_", &
          test_params%shape_name, "_", test_params%mesh_name, ".txt"
      open(newunit=fh_hf,  file=trim(adjustl(filename)))

      write(fh_hf,  '(a)') '# ncell l1 l2 linf'
    end if

    do i = 1,4
      ncell = 10 * 2**i
      call mesh_test (test_params, ncell, lnormFT, lnormHF, fh_fit, fh_hf)
    end do

    close(fh_fit)
    if (test_params%mesh_type == REG) close(fh_hf)

    print '(a)', '===================================================='
    print '(a)'

  end subroutine curvature_grid_refinement_study

  subroutine init_curvature_test (this, mesh_type, shape_type, error_cutoff)

    class(curvature_test), intent(out) :: this
    integer, intent(in) :: mesh_type, shape_type
    real(r8), intent(in) :: error_cutoff

    this%mesh_type = mesh_type
    this%shape_type = shape_type
    this%error_cutoff = error_cutoff

    select case (mesh_type)
    case(REG)
      this%mesh_name = 'reg'
    case(RND)
      this%mesh_name = 'rnd'
    case(TET)
      this%mesh_name = 'tet'
    end select

    select case (shape_type)
    case(SPHERE)
      this%shape_name = 'sphere'
      this%shape_filename = 'sphere.json'
    case(ELLIPSOID)
      this%shape_name = 'ellipsoid'
      this%shape_filename = 'ellipsoid.json'
    case(SINUSOID)
      this%shape_name = 'sinusoid'
      this%shape_filename = 'sinusoid.json'
    end select

  end subroutine init_curvature_test

  subroutine mesh_test(test_params, mesh_size, lnormFT, lnormHF, fh_fit, fh_hf)

    use,intrinsic :: iso_c_binding, only: C_NEW_LINE
    use unstr_mesh_factory
    use parameter_list_type
    use parameter_list_json
    use vof_init
    use array_utils, only: normalize, isZero
    use vof_io
    use mixed_cell_subset_constructor
    use mesh_subset_type

    type(curvature_test), intent(in) :: test_params
    integer, intent(in) :: mesh_size
    real(r8), intent(out) :: lnormFT(:), lnormHF(:)
    integer, intent(in) :: fh_fit, fh_hf

    character(:), allocatable :: errmsg, filename
    character(256) :: tmp
    type(unstr_mesh) :: mesh
    type(mesh_geom) :: gmesh
    type(mesh_subset) :: mixed_cells
    type(parameter_list), pointer :: plist
    real(r8), allocatable, target :: vof(:,:)
    real(r8), allocatable :: curvature_ex(:)
    integer :: infile

    ! get mesh
    if (test_params%mesh_type == REG) then
      ! create a regular 3D mesh
      mesh = new_unstr_mesh ([-0.5_r8, -0.5_r8, -0.5_r8], [0.5_r8, 0.5_r8, 0.5_r8], &
          [mesh_size,mesh_size,mesh_size])
    else
      ! read in exodus mesh
      write (tmp, '(a,i0,3a)') "data/meshes/cube_", mesh_size, "_", test_params%mesh_name, ".exo"
      filename = trim(adjustl(tmp))
      mesh = new_unstr_mesh (filename)
    end if
    call gmesh%init (mesh)
    if (test_params%mesh_type == RND) call recalculate_mesh_volumes(mesh,gmesh) ! TODO: remove this?
    allocate(vof(2,mesh%ncell))
    print '(a)', 'mesh initialized'

    ! initialize shape type
    ! TODO: initialize without parameter_list_type or without JSON file
    open(newunit=infile,file=test_params%shape_filename,action='read',access='stream')
    call parameter_list_from_json_stream (infile, plist, errmsg)
    if (.not.associated(plist)) call LS_fatal ("error reading input file:" // C_NEW_LINE // errmsg)
    close(infile)

    plist => plist%sublist('initial-vof')

    write (tmp, '(2a,i0,5a)') "data/fields/", &
        "vof_", mesh%ncell, "_", test_params%mesh_name, "_", test_params%shape_name, ".dat"
    filename = trim(adjustl(tmp))
    if (file_exists(filename)) then
      call read_vof_field(filename, vof)
    else
      print '(a)', 'initializing volume fractions...'
      call vof_initialize (mesh, gmesh, plist, vof, [1,2], 2)
      call store_vof_field(filename, vof)
      print '(a)', 'done'
    end if
    deallocate(plist)

    ! get mixed cells subset
    mixed_cells = mixed_cell_subset(vof, mesh)

    ! calculate errors for FT and HF curvature methods
    lnormFT = ft_mesh_test(test_params, vof, mesh, gmesh, mixed_cells, curvature_ex)
    write (fh_fit, '(i12,3es15.5)') mesh%ncell, lnormFT
    if (test_params%mesh_type == REG) then
      lnormHF = hf_mesh_test(test_params, vof, mesh, gmesh, curvature_ex)

      print '(i8, 2(a,3es10.2))', mesh%ncell, '  FT L1,L2,Linf = ',lnormFT, &
          ',  HF L1,L2,Linf = ',lnormHF
      write (fh_hf, '(i12,3es15.5)') mesh%ncell, lnormHF
    else
      lnormHF = 0

      print '(i8,a,3es10.2)', mesh%ncell, '  FT L1,L2,Linf = ',lnormFT
    end if

  end subroutine mesh_test

  function ft_mesh_test (test_params, vof, mesh, gmesh, mixed_cells, curvature_ex) result(lnorm)

    use interface_patch_type
    use multimat_cell_type
    use lvira_normals
    use mesh_subset_type
    use surface_type
    use vof_io

    type(curvature_test), intent(in) :: test_params
    real(r8), intent(in) :: vof(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    type(mesh_subset), intent(in) :: mixed_cells
    real(r8), allocatable, intent(inout) :: curvature_ex(:)
    real(r8) :: lnorm(3)

    character(:), allocatable :: filename
    character(256) :: tmp
    integer :: i, nvofcell, ierr, imax, c
    type(multimat_cell) :: cell
    type(surface) :: intrec
    real(r8) :: err, curvature(mesh%ncell), totvolume, xc(3)
    real(r8), allocatable :: int_norm(:,:,:), interface_centroid(:,:), reconstruction_area(:)

    write (tmp, '(2a,i0,5a)') "data/fields/", &
        "normals_", mesh%ncell, "_", test_params%mesh_name, "_", test_params%shape_name, ".dat"
    filename = trim(adjustl(tmp))
    if (file_exists(filename)) then
      allocate(int_norm(3,2,mesh%ncell))
      call read_normals(filename, int_norm)
    else
      print '(a)', 'calculating normals ...'
      call interface_normals_lvira(int_norm, vof, mesh, gmesh, mixed_cells)
      call store_normals(filename, int_norm)
      print '(a)', 'done'
    end if

    allocate(interface_centroid(3,mixed_cells%ncell), reconstruction_area(mixed_cells%ncell))

    ! get the interface reconstructions
    lnorm = 0; nvofcell = 0; totvolume = 0
    do c = 1,mixed_cells%ncell
      i = mixed_cells%cell_id(c)

      call cell%init (ierr, i, mesh, gmesh, tesselate=(test_params%mesh_type/=REG))
      if (ierr /= 0) call LS_fatal ('cell_outward_volflux failed: could not initialize cell')

      call cell%partition (vof(:,i), int_norm(:,:,i))

      call intrec%append (cell%interface_polygons(1), i)

      interface_centroid(:,c) = cell%interface_polygons(1)%centroid()
      reconstruction_area(c) = cell%interface_polygons(1)%area()
    end do
    print '(a)', 'interfaces located'

    call compute_exact_curvature_field(test_params, curvature_ex, interface_centroid, int_norm, &
        mesh, gmesh, mixed_cells)

    ! get the curvature
    curvature = 0; lnorm = 0; totvolume = 0
    do i = 1,mesh%ncell
      ! skipping boundaries
      if (any(gmesh%cneighbor(:, &
          pack(gmesh%cneighbor(:,i), mask=gmesh%cneighbor(:,i)>0)) < 1)) cycle

      if (vof(1,i) > test_params%error_cutoff .and. vof(1,i) < 1-test_params%error_cutoff) then
        curvature(i) = curvature_from_patch (intrec%local_patch(i,gmesh, vof(1,:)), &
            0.0_r8, int_norm(:,1,i), vof(1,:), mesh, gmesh, i, centroid=xc)

        ! append error to norms
        if (test_params%shape_type == SINUSOID) then
          err = abs(curvature(i) - curvature_ex(i))
        else
          err = abs((curvature(i) - curvature_ex(i)) / curvature_ex(i))
        end if

        if (err > lnorm(3)) imax = i
        totvolume = totvolume + mesh%volume(i)
        lnorm(1) = lnorm(1) + err * mesh%volume(i)
        lnorm(2) = lnorm(2) + err**2 * mesh%volume(i)
        lnorm(3) = max(lnorm(3),err)
      end if
    end do
    lnorm(1) = lnorm(1) / totvolume
    lnorm(2) = sqrt(lnorm(2) / totvolume)

  end function ft_mesh_test

  function hf_mesh_test (test_params, vof, mesh, gmesh, curvature_ex) result(lnorm)

    use array_utils, only: isZero
    use int_norm_module
    use curvature_hf

    type(curvature_test), intent(in) :: test_params
    real(r8), intent(in) :: vof(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8), intent(in) :: curvature_ex(:)
    real(r8) :: lnorm(3)

    integer :: i, nvofcell, ierr, imax
    real(r8) :: err, curvature(mesh%ncell)
    real(r8), allocatable :: int_norm(:,:,:)

    int_norm = interface_normal (vof, mesh, gmesh, .false.) ! get the initial normal estimate
    curvature = curvatureHF(vof(1,:), int_norm(:,1,:), mesh, gmesh)

    ! calculate error
    lnorm = 0; nvofcell = 0
    do i = 1,mesh%ncell
      if (any(gmesh%cneighbor(:,i)<1)) cycle ! skipping boundaries

      if (vof(1,i) > test_params%error_cutoff .and. vof(1,i) < 1-test_params%error_cutoff &
          .and. .not.isZero(curvature(i))) then

        if (test_params%shape_type == SINUSOID) then
          err = abs(curvature(i) - curvature_ex(i))
        else
          err = abs((curvature(i) - curvature_ex(i)) / curvature_ex(i))
        end if

        if (err > lnorm(3)) imax = i
        nvofcell = nvofcell + 1
        lnorm(1) = lnorm(1) + err
        lnorm(2) = lnorm(2) + err**2
        lnorm(3) = max(lnorm(3),err)
      end if
    end do
    lnorm(1) = lnorm(1) / nvofcell
    lnorm(2) = sqrt(lnorm(2) / nvofcell)

  end function hf_mesh_test

  subroutine recalculate_mesh_volumes(mesh,gmesh)

    use polyhedron_type

    type(unstr_mesh), intent(inout) :: mesh
    type(mesh_geom), intent(in) :: gmesh

    integer :: i, ierr
    type(polyhedron) :: cell

    print '(a)', 'recalculating mesh volumes ... '
    do i = 1,mesh%ncell
      call cell%init (ierr, i, mesh, gmesh)
      mesh%volume(i) = cell%volume()
    end do
    print '(a)', 'done'

  end subroutine recalculate_mesh_volumes

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

  end function ellipsoid_curvature_x

  real(r8) function sinusoid_curvature_x(dir, x, coeff)
    integer, intent(in) :: dir
    real(r8), intent(in) :: x(:), coeff(:)
    sinusoid_curvature_x = -pi**2*(coeff(2)*coeff(3)**2*(pi**2*coeff(5)**2*coeff(6)*&
        &*2*sin(pi*coeff(6)*(coeff(7) - x(2)))**2 + 1)*cos(pi*coeff(3)*(coeff(4) - x&
        &(1))) + coeff(5)*coeff(6)**2*(pi**2*coeff(2)**2*coeff(3)**2*sin(pi*coeff(3)&
        &*(coeff(4) - x(1)))**2 + 1)*cos(pi*coeff(6)*(coeff(7) - x(2))))*(pi**2*coef&
        &f(2)**2*coeff(3)**2*sin(pi*coeff(3)*(coeff(4) - x(1)))**2 + pi**2*coeff(5)*&
        &*2*coeff(6)**2*sin(pi*coeff(6)*(coeff(7) - x(2)))**2 + 1)**(-1.5_r8)
  end function sinusoid_curvature_x

  subroutine compute_exact_curvature_field(test_params, curvature_ex, interface_centroid, int_norm, &
      mesh, gmesh, mixed_cells)

    use mesh_subset_type
    use mixed_cell_subset_constructor

    type(curvature_test), intent(in) :: test_params
    real(r8), allocatable, intent(out) :: curvature_ex(:)
    real(r8), intent(in) :: interface_centroid(:,:), int_norm(:,:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    type(mesh_subset), intent(in) :: mixed_cells

    integer :: c, i
    real(r8), allocatable :: coeff(:)

    if (allocated(curvature_ex)) deallocate(curvature_ex)
    allocate(curvature_ex(mesh%ncell))
    curvature_ex = 0

    select case (test_params%shape_type)
    case(SPHERE)
      do c = 1,mixed_cells%ncell
        i = mixed_cells%cell_id(c)
        curvature_ex(i) = -2 / 0.35_r8 ! sphere
      end do
    case(ELLIPSOID)
      coeff = [0.35_r8, 0.3_r8, 0.2_r8]
      do c = 1,mixed_cells%ncell
        i = mixed_cells%cell_id(c)
        curvature_ex(i) = ellipsoid_curvature_x(maxloc(abs(int_norm(:,1,i)),1), &
            !nearest_ellipsoid_point(interface_centroid(:,c), int_norm(:,1,i), coeff), &
            !nearest_ellipsoid_point(gmesh%xc(:,i), int_norm(:,1,i), coeff), &
            !interface_centroid(:,c), &
            gmesh%xc(:,i), &
            coeff)
      end do
    case(SINUSOID)
      coeff = [0.0_r8, 0.125_r8, 2.5_r8, 0.2_r8, 0.125_r8, 2.5_r8, 0.2_r8]
      do c = 1,mixed_cells%ncell
        i = mixed_cells%cell_id(c)
        curvature_ex(i) = sinusoid_curvature_x(maxloc(abs(int_norm(:,1,i)),1), &
            !interface_centroid(:,c), &
            gmesh%xc(:,i), &
            coeff)
      end do
    end select

  end subroutine compute_exact_curvature_field

  ! find the closest on-ellipsoid point to x, along the direction n
  function nearest_ellipsoid_point(x, n, coeff) result(xr)

    real(r8), intent(in) :: x(:), n(:), coeff(:)
    real(r8) :: xr(3)

    real(r8) :: l, t, a, b, c

    ! coefficients in quadratic formula for the parameter l
    a = sum(n**2 / coeff**2)
    b = 2 * sum(x * n / coeff**2)
    c = sum(x**2 / coeff**2) - 1

    ! parameter l is the smallest (absolute value) solution to the quadratic
    ! unless we are seriously underresolved, one solution will be much smaller
    ! than the others and there will be no imaginary component
    l = -b + sqrt(b**2 - 4*a*c)
    t = -b - sqrt(b**2 - 4*a*c)
    if (abs(t) < abs(l)) l = t
    l = l / (2*a)

    xr = x + l*n

  end function nearest_ellipsoid_point

end module curvature_tests
