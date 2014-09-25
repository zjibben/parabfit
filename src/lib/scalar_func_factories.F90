!!
!! SCALAR_FUNC_FACTORIES
!!
!! Procedures for instantiating new SCALAR_FUNC class objects.  All return
!! return a new CLASS(SCALAR_FUNC) object, but the dynamic type is determined
!! by the particular function and arguments.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! PROGRAMMING INTERFACE
!!
!! Two flavors of procedures are provided.  In one, the functions return a
!! pointer to a newly allocated CLASS(SCALAR_FUNC) object, and in the second,
!! subroutines allocate an allocatable CLASS(SCALAR_FUNC) argument.  It is
!! typical for a CLASS(SCALAR_FUNC) object to be passed from one procedure
!! to another before ultimately being stored as a data structure component.
!! Pointers have been used to do this to avoid copying of the objects, however
!! the new MOVE_ALLOC instrinsic allows allocatable variables to be used for
!! the same purpose, bringing the advantages of allocatable over pointer.
!!

#include "f90_assert.fpp"

module scalar_func_factories

  use kinds, only: r8
  use scalar_func_class
  implicit none
  private

  public :: scalar_func ! re-export

  !! These return CLASS(SCALAR_FUNC) pointers
  public :: new_const_scalar_func
  public :: new_poly_scalar_func
  public :: new_mpoly_scalar_func
  public :: new_tabular_scalar_func
  public :: new_fptr_scalar_func

  !! These subroutines allocate an allocatable CLASS(SCALAR_FUNC) argument
  public :: alloc_const_scalar_func
  public :: alloc_poly_scalar_func
  public :: alloc_mpoly_scalar_func
  public :: alloc_tabular_scalar_func
  public :: alloc_fptr_scalar_func

  !! These higher-level procedures take a parameter list as input.
  public :: new_scalar_func
  public :: alloc_scalar_func

contains

  subroutine alloc_const_scalar_func (f, const)
    use const_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: const
    allocate(f, source=const_scalar_func(const))
  end subroutine alloc_const_scalar_func

  subroutine alloc_poly_scalar_func (f, c, e, x0)
    use poly_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:)
    real(r8), intent(in), optional :: x0
    allocate(f, source=poly_scalar_func(c, e, x0))
  end subroutine alloc_poly_scalar_func

  subroutine alloc_mpoly_scalar_func (f, c, e, x0)
    use mpoly_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:,:)
    real(r8), intent(in), optional :: x0(:)
    allocate(f, source=mpoly_scalar_func(c, e, x0))
  end subroutine alloc_mpoly_scalar_func

  subroutine alloc_tabular_scalar_func (f, x, y)
    use tabular_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    real(r8), intent(in) :: x(:), y(:)
    allocate(f, source=tabular_scalar_func(x, y))
  end subroutine alloc_tabular_scalar_func

  subroutine alloc_fptr_scalar_func (f, fptr, p)
    use fptr_scalar_func_type
    class(scalar_func), allocatable, intent(out) :: f
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    allocate(f, source=fptr_scalar_func(fptr, p))
  end subroutine alloc_fptr_scalar_func


  function new_const_scalar_func (const) result (f)
    use const_scalar_func_type
    real(r8), intent(in) :: const
    class(scalar_func), pointer :: f
    allocate(f, source=const_scalar_func(const))
  end function new_const_scalar_func

  function new_poly_scalar_func (c, e, x0) result (f)
    use poly_scalar_func_type
    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:)
    real(r8), intent(in), optional :: x0
    class(scalar_func), pointer :: f
    allocate(f, source=poly_scalar_func(c, e, x0))
  end function new_poly_scalar_func

  function new_mpoly_scalar_func (c, e, x0) result (f)
    use mpoly_scalar_func_type
    real(r8), intent(in) :: c(:)
    integer,  intent(in) :: e(:,:)
    real(r8), intent(in), optional :: x0(:)
    class(scalar_func), pointer :: f
    allocate(f, source=mpoly_scalar_func(c, e, x0))
  end function new_mpoly_scalar_func

  function new_tabular_scalar_func (x, y) result (f)
    use tabular_scalar_func_type
    real(r8), intent(in) :: x(:), y(:)
    class(scalar_func), pointer :: f
    allocate(f, source=tabular_scalar_func(x, y))
  end function new_tabular_scalar_func

  function new_fptr_scalar_func (fptr, p) result (f)
    use fptr_scalar_func_type
    procedure(fptr_func) :: fptr
    real(r8), intent(in), optional :: p(:)
    class(scalar_func), pointer :: f
    allocate(f, source=fptr_scalar_func(fptr, p))
  end function new_fptr_scalar_func

  !! A higher-level constructors that takes a parameter list
  !! that describes the function to be instantiated.

  function new_scalar_func (params) result (f)

    use parameter_list_type
    use logging_services

    type(parameter_list), intent(inout) :: params
    class(scalar_func), pointer :: f

    real(r8) :: x0, t0, v0, t1, v1
    integer, allocatable :: expon(:), expon2(:,:)
    real(r8), allocatable :: coef(:), x(:), y(:), x0_def(:), x02(:)
    character(:), allocatable :: ftype, context, errmsg
    integer :: stat

    context = 'processing ' // params%name() // ': '
    call params%get ('type', ftype, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    select case (ftype)
    case ('constant')
      call params%get ('value', v0, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      f => new_const_scalar_func(v0)
    case ('polynomial')
      call params%get ('poly-coef', coef, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (params%is_vector('poly-powers')) then
        call params%get ('poly-powers', expon, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        call params%get ('poly-center', x0, default=0.0_r8, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        if (size(coef) /= size(expon)) &
            call LS_fatal (context//'inconsistent sizes for "poly-coef" and "poly-powers"')
        f => new_poly_scalar_func(coef, expon, x0)
      else if (params%is_matrix('poly-powers')) then
        call params%get ('poly-powers', expon2, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        if (size(coef) /= size(expon2,2)) &
            call LS_fatal (context//'inconsistent sizes for "poly-coef" and "poly-powers"')
        allocate(x0_def(size(expon2,1)))
        x0_def = 0.0_r8
        call params%get ('poly-center', x02, default=x0_def, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        f => new_mpoly_scalar_func(coef, expon2, x02)
      else
        call LS_fatal (context//'"poly-powers" must be vector or matrix-valued')
      end if
    case ('tabular')
      call params%get ('tabular-x', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('tabular-y', y, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (size(x) /= size(y)) &
          call LS_fatal (context//'inconsistent sizes for "tabular-x" and "tabular-y"')
      f => new_tabular_scalar_func(x, y)
    case ('smooth-ramp')
      call params%get ('begin-time', t0, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('begin-value', v0, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('end-time', t1, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('end-value', v1, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      f => new_fptr_scalar_func(smooth_ramp,[t0,v0,t1,v1])
    case ('sinusoidal')
      call params%get ('mean-value', v0, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('amplitude',  v1, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('frequency',  t1, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('time-shift', t0, default=0.0_r8, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      f => new_fptr_scalar_func(sinusoid,[t0,t1,v0,v1])
    case default
      call LS_fatal (context//'unknown "type" value: '//ftype)
    end select
    ASSERT(associated(f))

  end function new_scalar_func


  subroutine alloc_scalar_func (f, params)

    use parameter_list_type
    use logging_services

    class(scalar_func), allocatable, intent(out) :: f
    type(parameter_list), intent(inout) :: params

    real(r8) :: x0, t0, v0, t1, v1
    integer, allocatable :: expon(:), expon2(:,:)
    real(r8), allocatable :: coef(:), x(:), y(:), x0_def(:), x02(:)
    character(:), allocatable :: ftype, context, errmsg
    integer :: stat

    context = 'processing ' // params%name() // ': '
    call params%get ('type', ftype, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    select case (ftype)
    case ('constant')
      call params%get ('value', v0, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call alloc_const_scalar_func (f, v0)
    case ('polynomial')
      call params%get ('poly-coef', coef, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (params%is_vector('poly-powers')) then
        call params%get ('poly-powers', expon, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        call params%get ('poly-center', x0, default=0.0_r8, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        if (size(coef) /= size(expon)) &
            call LS_fatal (context//'inconsistent sizes for "poly-coef" and "poly-powers"')
        call alloc_poly_scalar_func (f, coef, expon, x0)
      else if (params%is_matrix('poly-powers')) then
        call params%get ('poly-powers', expon2, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        if (size(coef) /= size(expon2,2)) &
            call LS_fatal (context//'inconsistent sizes for "poly-coef" and "poly-powers"')
        allocate(x0_def(size(expon2,1)))
        x0_def = 0.0_r8
        call params%get ('poly-center', x02, default=x0_def, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        call alloc_mpoly_scalar_func (f, coef, expon2, x02)
      else
        call LS_fatal (context//'"poly-powers" must be vector or matrix-valued')
      end if
    case ('tabular')
      call params%get ('tabular-x', x, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('tabular-y', y, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (size(x) /= size(y)) &
          call LS_fatal (context//'inconsistent sizes for "tabular-x" and "tabular-y"')
      call alloc_tabular_scalar_func (f, x, y)
    case ('smooth-ramp')
      call params%get ('begin-time', t0, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('begin-value', v0, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('end-time', t1, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('end-value', v1, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call alloc_fptr_scalar_func (f, smooth_ramp,[t0,v0,t1,v1])
    case ('sinusoidal')
      call params%get ('mean-value', v0, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('amplitude',  v1, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('frequency',  t1, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call params%get ('time-shift', t0, default=0.0_r8, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call alloc_fptr_scalar_func(f, sinusoid, [t0,t1,v0,v1])
    case default
      call LS_fatal (context//'unknown "type" value: '//ftype)
    end select
    ASSERT(allocated(f))

  end subroutine alloc_scalar_func


  real(r8) function smooth_ramp (x, p)
    real(r8), intent(in) :: x(*), p(*)
    real(r8) :: s
    associate (t => x(1), t0 => p(1), v0 => p(2), t1 => p(3), v1 => p(4))
      if (t <= t0) then
        smooth_ramp = v0
      else if (t >= t1) then
        smooth_ramp = v1
      else
        s = (t - t0) / (t1 - t0)
        smooth_ramp = v0 + (v1-v0) * (s**2) * (3 - 2*s)
      end if
    end associate
  end function smooth_ramp

  real(r8) function sinusoid (x, p)
    real(r8), intent(in) :: x(*), p(*)
    real(r8), parameter :: TWOPI = 6.2831853071795864769_r8
    associate (t => x(1), t0 => p(1), w => p(2), m => p(3), a => p(4))
      sinusoid = m + a*sin(TWOPI*w*(t-t0))
    end associate
  end function sinusoid

end module scalar_func_factories
