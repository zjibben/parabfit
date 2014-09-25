!!
!! HC_IDAESOL_MODEL_TYPE
!!
!! This module defines an extension of the IDAESOL_MODEL abstract class that
!! implements the methods required by the ODE integrator.  It bundles up
!! several different computational pieces and adapts them to the required
!! interface.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! August 2014
!!

#include "f90_assert.fpp"

module HC_idaesol_model_type

  use kinds, only: r8
  use HC_model_type
  use HC_precon_type
  use HC_norm_type
  use idaesol_type, only: idaesol_model
  implicit none
  private

  type, extends(idaesol_model), public :: HC_idaesol_model
    type(HC_model),  pointer :: model  => null()  ! reference only -- do not own
    type(HC_precon), pointer :: precon => null()  ! reference only -- do not own
    type(HC_norm),   pointer :: norm   => null()  ! reference only -- do not own
  contains
    procedure :: init
    !! Deferred procedures from IDAESOL_MODEL
    procedure :: size => model_size
    procedure :: compute_f
    procedure :: apply_precon
    procedure :: compute_precon
    procedure :: du_norm
    procedure :: schk
  end type HC_idaesol_model

contains

  subroutine init (this, model, precon, norm)
    class(HC_idaesol_model), intent(out) :: this
    class(HC_model),  intent(in), target :: model
    class(HC_precon), intent(in), target :: precon
    class(HC_norm),   intent(in), target :: norm
    this%model => model
    this%precon => precon
    this%norm => norm
    ASSERT(associated(this%model, precon%model))
    !ASSERT(associated(this%model, norm%model))
  end subroutine init

  integer function model_size (this)
    class(HC_idaesol_model), intent(in) :: this
    model_size = this%model%num_dof()
  end function model_size

  subroutine compute_f (this, t, u, udot, f)
    class(HC_idaesol_model) :: this
    real(r8), intent(in)  :: t, u(:), udot(:)
    real(r8), intent(out) :: f(:)
    call this%model%residual (t, u, udot, f)
  end subroutine compute_f

  subroutine apply_precon (this, t, u, f)
    class(HC_idaesol_model) :: this
    real(r8), intent(in) :: t, u(:)
    real(r8), intent(inout) :: f(:)
    call this%precon%apply (f)
  end subroutine apply_precon

  subroutine compute_precon (this, t, u, dt)
    class(HC_idaesol_model) :: this
    real(r8), intent(in) :: t, u(:), dt
    call this%precon%compute (t, u, dt)
  end subroutine compute_precon

  subroutine du_norm (this, u, du, error)
    class(HC_idaesol_model) :: this
    real(r8), intent(in) :: u(:), du(:)
    real(r8), intent(out) :: error
    call this%norm%compute (u, du, error)
  end subroutine du_norm

  subroutine schk (this, u, stage, errc)
    class(HC_idaesol_model) :: this
    real(r8), intent(in) :: u(:)
    integer, intent(in)  :: stage
    integer, intent(out) :: errc
    errc = 0  ! solution is fine
  end subroutine schk

end module HC_idaesol_model_type
