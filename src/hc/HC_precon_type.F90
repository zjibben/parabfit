!!
!! HC_PRECON_TYPE
!!
!! This module defines a derived type that encapsulates the preconditioner for
!! the heat conduction model.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, April 2014
!!
!! PROGRAMMING INTERFACE
!!
!! The module defines the derived type HC_PRECON_TYPE.  It has the following
!! type bound procedures.
!!
!!  INIT(MODEL, PARAMS) initializes the object.  MODEL is of type HC_MODEL.
!!    The object will hold a reference to the model, and so the actual argument
!!    must have the target attribute and persist for the lifetime of the object.
!!    The PARAMETER_LIST type argument PARAMS gives the parameters for the
!!    preconditioner.  For this model there is only the single heat equation
!!    and the expected parameters are those described for MFD_DIFF_PRECON.
!!
!!  COMPUTE(T, U, DT) computes the preconditioner for the model at time T,
!!    unknown vector U, and time step DT.  It must be called before calling
!!    the APPLY procedure.
!!
!!  APPLY(F) applies the preconditioner for the model to the vector F, which is
!!    overwritten with the result.
!!

#include "f90_assert.fpp"

module HC_precon_type

  use kinds, only: r8
  use HC_model_type
  use unstr_mesh_type
  use mfd_diff_precon_type
  use mfd_diff_matrix_type
  use parameter_list_type
  use timer_tree_type
  implicit none
  private

  type, public :: HC_precon
    type(HC_model),   pointer :: model => null()  ! reference only -- do not own
    type(unstr_mesh), pointer :: mesh  => null()  ! reference only -- do not own
    real(r8) :: dt  ! time step
    real(r8), allocatable :: dHdT(:)  ! derivative of the enthalpy/temperature relation
    type(mfd_diff_precon) :: hcprecon ! heat equation preconditioner
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type HC_precon

contains

  subroutine init (this, model, params)

    class(HC_precon), intent(out) :: this
    type(HC_model), intent(in), target :: model
    type(parameter_list) :: params

    type(mfd_diff_matrix), allocatable :: dm

    this%model => model
    this%mesh  => model%mesh

    !! Heat density/temperature relation derivative.
    allocate(this%dHdT(this%mesh%ncell))

    !! Create the preconditioner for the heat equation.
    !! The preconditioner assumes ownership of the matrix.
    allocate (dm)
    call dm%init (model%disc)
    call this%hcprecon%init (dm, params)

  end subroutine init

  subroutine compute (this, t, u, dt)

    class(HC_precon), intent(inout) :: this
    real(r8), intent(in) :: t, dt
    real(r8), intent(in), target :: u(:)

    real(r8), pointer :: temp(:)
    real(r8) :: coef(this%mesh%ncell)
    type(mfd_diff_matrix), pointer :: dm

    call start_timer ('hc-precon-compute')
    call start_timer ('matrix-compute')

    ASSERT(size(u) == this%model%num_dof())
    ASSERT(dt > 0.0_r8)

    this%dt = dt  ! the time step size

    !! Grab the matrix; we will update its values.
    dm => this%hcprecon%matrix()

    !! Jacobian of the heat equation diffusion operator ignoring nonlinearities.
    call this%model%get_cell_temp_view (u, temp)
    call this%model%conductivity (temp, coef)
    call dm%compute (coef)

    !! Contribution from the time derivative (H/T relation eliminated).
    call this%model%dHdT (temp, this%dHdT)
    call dm%incr_cell_diag (this%mesh%volume * this%dHdT / dt)

    !! Dirichlet boundary condition fixups.
    call this%model%temp_bc%compute (t)
    call dm%set_dir_faces (this%model%temp_bc%index)
    call stop_timer ('matrix-compute')

    !! The matrix is now complete; re-compute the preconditioner.
    call this%hcprecon%compute

    call stop_timer ('hc-precon-compute')

  end subroutine compute

  subroutine apply (this, f)

    class(HC_precon), intent(in) :: this
    real(r8), intent(inout), target :: f(:)

    real(r8), pointer :: f0(:), f1(:), f2(:)

    call start_timer ('hc-precon-apply')

    !! Heat equation cell residual with the H/T relation residual eliminated.
    call this%model%get_cell_heat_view (f, f0)
    call this%model%get_cell_temp_view (f, f1)
    f1 = f1 - (this%mesh%volume/this%dt)*f0

    !! Heat equation face residual.
    call this%model%get_face_temp_view (f, f2)

    !! Precondition the heat equation.
    call this%hcprecon%apply (f1, f2)

    !! Backsubstitute to obtain the preconditioned H/T-relation residual.
    f0 = f0 + this%dHdT*f1

    call stop_timer ('hc-precon-apply')

  end subroutine apply

end module HC_precon_type
