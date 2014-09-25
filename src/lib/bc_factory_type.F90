!!
!! BC_FACTORY_TYPE
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2014
!!

#include "f90_assert.fpp"

module bc_factory_type

  use unstr_mesh_type
  use parameter_list_type
  implicit none
  private

  type, public :: bc_factory
    private
    type(unstr_mesh),     pointer :: mesh   => null() ! reference only - do not own
    type(parameter_list), pointer :: params => null() ! reference only - do not own
  contains
    procedure :: init
    procedure :: alloc_bc
  end type

contains

  subroutine init (this, mesh, params)
    class(bc_factory), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list), intent(in), target :: params
    this%mesh => mesh
    this%params => params
  end subroutine init

  subroutine alloc_bc (this, condition, bd)

    use kinds, only: r8
    use scalar_func_class
    use scalar_func_factories
    use bndry_func_class
    use bndry_face_func_type
    use logging_services
    use parameter_list_type !BUG: not necessary but Intel needs it (should be reported)

    class(bc_factory), intent(in) :: this
    character(*), intent(in) :: condition
    class(bndry_func), allocatable, intent(out) :: bd

    type(bndry_face_func), allocatable :: bc
    type(parameter_list_iterator) :: piter
    type(parameter_list), pointer :: bc_params, func_params
    character(:), allocatable :: bc_type
    class(scalar_func), allocatable :: f
    integer, allocatable :: setids(:)
    integer :: stat
    character(:), allocatable :: context, errmsg
    real(r8) :: value

    call LS_info ('Generating "' // condition // '" boundary condition', LS_VERB_NOISY)

    allocate(bc)
    call bc%init (this%mesh)

    piter = parameter_list_iterator(this%params, sublists_only=.true.)
    do while (.not.piter%at_end())
      bc_params => piter%sublist()
      ASSERT(associated(bc_params))
      context = 'processing ' // bc_params%name() // ': '
      call bc_params%get('condition', bc_type, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (bc_type == condition) then
        call LS_info ('  using "' // piter%name() // '"', LS_VERB_NOISY)
        if (bc_params%is_parameter('data-constant')) then
          call bc_params%get('data-constant', value, stat=stat, errmsg=errmsg)
          if (stat /= 0) call LS_fatal (context//errmsg)
          call alloc_const_scalar_func (f, value)
        else if (bc_params%is_sublist('data-function')) then
          func_params => bc_params%sublist('data-function', stat=stat, errmsg=errmsg)
          if (stat /= 0) call LS_fatal (context//errmsg)
          call alloc_scalar_func (f, func_params)
        else
          call LS_fatal (context//'missing "data-constant" or "data-function" parameter')
        end if
        INSIST(allocated(f))
        call bc_params%get ('face-sets', setids, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        call bc%add (f, setids, stat, errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
      end if
      call piter%next
    end do
    call bc%add_complete

    call move_alloc (bc, bd)

  end subroutine alloc_bc

end module bc_factory_type
