#include "f90_assert.fpp"

module region_factories

  use kinds, only: r8
  use region_class
  implicit none
  private

  !! higher-level procedure which takes a parameter list as input
  public :: alloc_region

contains

  subroutine alloc_plane_region (r, n, p)
    use plane_region_type
    class(region), allocatable, intent(out) :: r
    real(r8), intent(in) :: n(:), p
    allocate(r, source=plane_region(n, p))
  end subroutine alloc_plane_region

  subroutine alloc_sphere_region (r, xc, radius)
    use sphere_region_type
    class(region), allocatable, intent(out) :: r
    real(r8), intent(in) :: xc(:), radius
    allocate(r, source=sphere_region(xc, radius))
  end subroutine alloc_sphere_region

  subroutine alloc_halfsphere_region (r, xc, radius, n)
    use halfsphere_region_type
    class(region), allocatable, intent(out) :: r
    real(r8), intent(in) :: xc(:), radius, n(:)
    allocate(r, source=halfsphere_region(xc, radius, n))
  end subroutine alloc_halfsphere_region

  subroutine alloc_fill_region (r)
    use fill_region_type
    class(region), allocatable, intent(out) :: r
    allocate(fill_region :: r)
  end subroutine alloc_fill_region

  subroutine alloc_region (r, params)

    use parameter_list_type
    use logging_services
    use array_utils, only: normalize

    class(region), allocatable, intent(out) :: r
    type(parameter_list), intent(inout) :: params

    real(r8), allocatable :: normal(:), x(:)
    real(r8) :: p
    character(:), allocatable :: rtype, context, errmsg
    integer :: stat

    context = 'processing ' // params%name() // ': '
    call params%get ('type', rtype, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    select case (rtype)
    case ('sphere')
      call params%get ('center', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('radius', p, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      ASSERT(size(x)==3)
      call alloc_sphere_region (r, x, p)
    case ('plane')
      call params%get ('normal', normal, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('plane-const', p, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      ASSERT(size(normal)==3)
      ! make sure the normal is normalized
      normal = normalize(normal)
      call alloc_plane_region (r, normal, p)
    case ('halfsphere', 'half-sphere')
      call params%get ('normal', normal, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('center', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('radius', p, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      ASSERT(size(x)==3)
      ASSERT(size(normal)==3)
      normal = normalize(normal)
      call alloc_halfsphere_region (r, x, p, normal)
    case ('fill', 'all')
      call alloc_fill_region (r)
    case default
      call LS_fatal (context//'unknown "type" value: '//rtype)
    end select
    ASSERT(allocated(f))

  end subroutine alloc_region

end module region_factories
