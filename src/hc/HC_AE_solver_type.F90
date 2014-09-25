!!
!! HC_AE_SOLVER_TYPE
!!
!! A solver for the algebraic equation portion of the heat conduction model
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! September 2014
!!
!! NOTES
!!
!!  This is a rough first pass and much work remains to clean this up.
!!  In particular the following points should be addressed.
!!
!!  1. Only the (2,2) (face-face) block of the diffusion matrix is needed
!!  but we are hijacking the machinery of the full diffusion matrix to
!!  get it.  Significant wasted memory here.  We are also hijacking the
!!  full HC residual when only the algebraic portion is needed.
!!
!!  2. This potentially has several use cases beyond just solving for
!!  the initial face temperatures; e.g., making them exactly consistent
!!  as part of the nonlinear time step solves (the residual in the
!!  algebraic equation is exactly the local face flux mismatch).  For
!!  this, some further attention needs to be paid to the different
!!  phases of the calculation, and exposing them individually to avoid
!!  unnecessary repeated work when doing multiple solves, as opposed to
!!  a single one-off solve.
!!
!!  4. The calculation of the matrix probably should be a HC_model
!!  method.  This is because it depends on BC, and otherwise code that
!!  depends on BC will be scattered all over the place, where ever the
!!  diffusion matrix is used (preconditioners, for example).
!!
!!  5. The HC_model residual method is being used to get the rhs for the
!!  algebraic system.  A significant portion of the residual is ignored
!!  and it also requires udot (we pass 0) -- more wasted memory here.
!!  It would be nice to avoid this if possible, but not if that means
!!  replicated residual code.
!!

#include "f90_assert.fpp"

module HC_AE_solver_type

  use kinds, only: r8
  use HC_model_type
  use mfd_diff_matrix_type
  use hypre_hybrid_type
  use parameter_list_type
  use logging_services
  implicit none
  private

  type, public :: HC_AE_solver
    private
    type(HC_model), pointer :: model => null()  ! reference only -- do not own
    type(parameter_list), pointer :: params => null() ! reference only - do not own
  contains
    procedure :: init
    procedure :: solve
  end type HC_AE_solver

contains

  subroutine init (this, model, params)
    class(HC_AE_solver), intent(out) :: this
    type(HC_model), intent(in), target :: model
    type(parameter_list), pointer :: params
    this%model => model
    this%params => params
  end subroutine init

  subroutine solve (this, t, tcell, tface, stat, errmsg)

    class(HC_AE_solver), intent(inout) :: this
    real(r8), intent(in) :: t, tcell(:)
    real(r8), intent(inout) :: tface(:)
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    type(mfd_diff_matrix), target :: dm
    type(hypre_hybrid) :: solver
    real(r8), allocatable :: coef(:), z(:), u(:), udot(:)
    real(r8), allocatable, target :: r(:)
    real(r8), pointer :: rface(:)
    real(r8) :: norm
    integer :: num_itr, num_dscg_itr, num_pcg_itr
    character(80) :: string

    ASSERT(size(tcell) == this%model%mesh%ncell)
    ASSERT(size(tface) == this%model%mesh%nface)

    !! Compute the RHS.
    allocate(u(this%model%num_dof()), udot(this%model%num_dof()), r(this%model%num_dof()))
    udot = 0.0_r8
    u = 0.0_r8
    call this%model%set_cell_temp (tcell, u)
    call this%model%set_face_temp (tface, u)
    call this%model%residual (t, u, udot, r)
    call this%model%get_face_temp_view (r, rface)
    norm = sqrt(sum(rface**2))
    if (LS_VERBOSITY >= LS_VERB_NOISY) then
      write(string,'(a,es9.2)') 'HC_AE_solver%solve: initial ||rface||_2 = ', norm
      call LS_info (string)
    end if
    if (norm == 0.0_r8) return

    !! Compute the matrix (the A22 submatrix)
    allocate(coef(size(tcell)))
    call this%model%conductivity (tcell, coef)
    call dm%init (this%model%disc)
    call dm%compute (coef)
    call this%model%temp_bc%compute (t)
    call dm%set_dir_faces (this%model%temp_bc%index)
    deallocate(coef)

    !! Setup the linear solver.
    call this%params%set ('krylov-method', 'cg')
    call this%params%get ('max-iter', num_itr)
    call this%params%set ('max-ds-iter', num_itr)
    call this%params%set ('max-amg-iter', num_itr)
    if (LS_VERBOSITY >= LS_VERB_NOISY) then
      call this%params%set ('print-level', 1)
      call this%params%set ('logging-level', 1)
    else
      call this%params%set ('print-level', 0)
      call this%params%set ('logging-level', 0)
    end if
    call solver%init (dm%a22, this%params)
    call solver%setup

    !! Solve
    allocate(z(size(tface)))
    z = 0.0_r8
    call solver%solve (rface, z, stat)
    tface = tface - z

    if (LS_VERBOSITY >= LS_VERB_NOISY) then
      call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(string,'(3(a,i0),a,es9.2)') 'HC_AE_solver%solve: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      call LS_info (string)

      !! Check the residual.
      call this%model%set_face_temp (tface, u)
      call this%model%residual (t, u, udot, r)
      write(string,'(a,es9.2)') 'HC_AE_solver%solve: ||rface||_2 = ', sqrt(sum(rface**2))
      call LS_info (string)
    end if

    if (stat /= 0) then
      call solver%get_metrics (num_itr, num_dscg_itr, num_pcg_itr, norm)
      write(string,'(3(a,i0),a,es9.2)') 'failed to converge: num_itr = ', num_itr, &
          ' (', num_dscg_itr, ', ', num_pcg_itr, '), ||r||/||b|| = ', norm
      errmsg = trim(string)
    end if

  end subroutine solve

end module HC_AE_solver_type
