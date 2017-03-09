!!
!! vof_init_test
!!

module vof_init_test

  use kinds, only: r8
  use vof_init
  implicit none
  private

  public :: vof_init_test_suite

contains

  subroutine vof_init_test_suite ()

    print '(a)'
    print '(a)', 'VOF_INIT'
    print '(a)', '===================================================='

    call cylinder_test ([-5.5e-2_r8, 0.245_r8, 0.0_r8], 1e-2_r8)

    print '(a)', '===================================================='
    print '(a)'

  end subroutine vof_init_test_suite

  ! check that the volume of a cut of a cylinder is correct
  subroutine cylinder_test (xc, dx)

    use region_class
    use material_geometry_type
    use region_factories, only: alloc_cylinder_region, alloc_fill_region
    use scalar_func_class
    use scalar_func_containers, only: scalar_func_box
    use scalar_func_factories, only: alloc_const_scalar_func, alloc_piecewise_scalar_func

    real(r8), intent(in) :: xc(:), dx

    type(dnc_hex) :: hex
    type(material_geometry) :: matl_geom
    type(region_box) :: rgn(2)
    type(scalar_func_box) :: subfunc(2)
    class(scalar_func), allocatable :: matl_index
    real(r8) :: vof(2), vof_ex

    ! set up the cylinder geometry
    call alloc_cylinder_region (rgn(1)%r, [0.0_r8, 0.0_r8, 0.0_r8], [0.0_r8, 0.0_r8, 1.0_r8], &
        0.25_r8, 1.0_r8)
    call alloc_fill_region (rgn(2)%r)

    call alloc_const_scalar_func (subfunc(1)%f, 1.0_r8)
    call alloc_const_scalar_func (subfunc(2)%f, 2.0_r8)

    call alloc_piecewise_scalar_func (matl_index, subfunc, rgn)
    call matl_geom%init (matl_index)

    ! set up cell
    hex%node(:,1) = xc + 0.5_r8 * dx * [-1.0_r8, -1.0_r8, -1.0_r8]
    hex%node(:,2) = xc + 0.5_r8 * dx * [ 1.0_r8, -1.0_r8, -1.0_r8]
    hex%node(:,3) = xc + 0.5_r8 * dx * [ 1.0_r8,  1.0_r8, -1.0_r8]
    hex%node(:,4) = xc + 0.5_r8 * dx * [-1.0_r8,  1.0_r8, -1.0_r8]
    hex%node(:,5) = xc + 0.5_r8 * dx * [-1.0_r8, -1.0_r8,  1.0_r8]
    hex%node(:,6) = xc + 0.5_r8 * dx * [ 1.0_r8, -1.0_r8,  1.0_r8]
    hex%node(:,7) = xc + 0.5_r8 * dx * [ 1.0_r8,  1.0_r8,  1.0_r8]
    hex%node(:,8) = xc + 0.5_r8 * dx * [-1.0_r8,  1.0_r8,  1.0_r8]

    ! calculate the vof
    vof = hex%vof (matl_geom, 2, 0, tiny(1.0_r8))

    ! print results
    vof_ex = 0.614298769812239_r8
    print '(a,3es15.5)', "vof, vofex, err: ", vof(2), vof_ex, abs(vof(2) - vof_ex)
    
  end subroutine cylinder_test

end module vof_init_test
