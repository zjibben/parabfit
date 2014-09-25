!!
!! IDAESOL_TYPE
!!
!! A solver for index-1 DAE in implicit form using the BDF2 method.
!!
!! This F2008 version is adapted from much earlier F95 and F77 implementations
!! and was obtained from http://sourceforge.net/projects/calliope.mfeproject.p/
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) 2013  Neil N. Carlson <neil.n.carlson@gmail.com>
!!
!! Permission is hereby granted, free of charge, to any person obtaining a
!! copy of this software and associated documentation files (the "Software"),
!! to deal in the Software without restriction, including without limitation
!! the rights to use, copy, modify, merge, publish, distribute, sublicense,
!! and/or sell copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following conditions:
!!
!! The above copyright notice and this permission notice shall be included
!! in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!! THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!! DEALINGS IN THE SOFTWARE.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "f90_assert.fpp"

module idaesol_type

  use,intrinsic :: iso_fortran_env, only: r8 => real64
  use state_history_type
  use nka_type
  implicit none
  private

  type, public :: idaesol
    private
    class(idaesol_model), pointer :: model => null()
    integer  :: n                   ! number of unknowns
    integer  :: seq = -1            ! number of steps taken
    real(r8) :: hlast               ! last step size
    real(r8) :: hpc                 ! step size built into the current preconditioner
    logical  :: usable_pc = .false. ! whether the current preconditioner is usable
    integer  :: freeze_count = 0    ! don't increase step size for this number of steps
    integer  :: mitr = 5            ! maximum number of nonlinear iterations
    real(r8) :: ntol = 0.1_r8       ! nonlinear solver error tolerance (relative to 1)
    type(nka) :: nka          ! nonlinear solver (AIN) accelerator structure
    type(state_history)   :: uhist        ! solution history structure

    !! Perfomance counters
    integer :: pcfun_calls = 0      ! number of calls to PCFUN
    integer :: updpc_calls = 0      ! number of calls to UPDPC
    integer :: updpc_failed = 0     ! number of UPDPC calls returning an error
    integer :: retried_bce = 0      ! number of retried BCE steps
    integer :: failed_bce = 0       ! number of completely failed BCE steps
    integer :: rejected_steps = 0   ! number of steps rejected on error tolerance
    real(r8) :: hmin = huge(1.0_r8) ! minimum step size used on a successful step
    real(r8) :: hmax = tiny(1.0_r8) ! maximum step size used on a successful step

    !! Diagnostics
    integer :: unit = 0
    logical :: verbose = .false.
  contains
    procedure :: init
    procedure :: set_initial_state
    procedure :: integrate => bdf2_step_driver
    procedure :: step
    procedure :: commit_state
    procedure :: get_interpolated_state
    procedure :: get_last_state_copy
    procedure :: get_last_state_view
    procedure :: last_time
    procedure :: last_step_size
    procedure :: set_verbose_stepping
    procedure :: set_quiet_stepping
    generic :: write_metrics => write_metrics_unit, write_metrics_string
    procedure, private :: write_metrics_unit
    procedure, private :: write_metrics_string
  end type idaesol

  type, abstract, public :: idaesol_model
  contains
    procedure(model_size), deferred :: size
    procedure(compute_f), deferred :: compute_f
    procedure(apply_precon), deferred :: apply_precon
    procedure(compute_precon), deferred :: compute_precon
    procedure(du_norm), deferred :: du_norm
    procedure(schk),  deferred :: schk
  end type

  abstract interface
    integer function model_size (this)
      import idaesol_model
      class(idaesol_model), intent(in) :: this
    end function model_size
    subroutine compute_f (this, t, u, udot, f)
      import idaesol_model, r8
      class(idaesol_model) :: this
      real(r8), intent(in)  :: t, u(:), udot(:)
      real(r8), intent(out) :: f(:)
    end subroutine compute_f
    subroutine apply_precon (this, t, u, f)
      import idaesol_model, r8
      class(idaesol_model) :: this
      real(r8), intent(in) :: t, u(:)
      real(r8), intent(inout) :: f(:)
    end subroutine apply_precon
    subroutine compute_precon (this, t, u, dt)
      import idaesol_model, r8
      class(idaesol_model) :: this
      real(r8), intent(in)  :: t, u(:), dt
    end subroutine compute_precon
    subroutine du_norm (this, u, du, error)
      import :: idaesol_model, r8
      class(idaesol_model) :: this
      real(r8), intent(in) :: u(:), du(:)
      real(r8), intent(out) :: error
    end subroutine du_norm
    subroutine schk (this, u, stage, errc)
      import :: idaesol_model, r8
      class(idaesol_model) :: this
      real(r8), intent(in)  :: u(:)
      integer,  intent(in)  :: stage
      integer,  intent(out) :: errc
    end subroutine schk
  end interface

  real(r8), parameter, private :: RMIN = 0.25_r8
  real(r8), parameter, private :: RMAX = 4.0_r8
  real(r8), parameter, private :: MARGIN = 3.0_r8

  !! Successful STATUS return codes:
  integer, parameter, public :: SOLVED_TO_TOUT = 1
  integer, parameter, public :: SOLVED_TO_NSTEP = 2

  !! Unsuccessful STATUS return codes:
  integer, parameter, public :: BAD_INPUT = -1
  integer, parameter, public :: STEP_FAILED = -2
  integer, parameter, public :: STEP_SIZE_TOO_SMALL = -3

contains

  subroutine get_last_state_view (this, view)
    class(idaesol), intent(in) :: this
    real(r8), pointer, intent(out) :: view(:)
    call this%uhist%get_last_state_view (view)
  end subroutine get_last_state_view

  subroutine get_last_state_copy (this, copy)
    class(idaesol), intent(in) :: this
    real(r8), intent(out) :: copy(:)
    call this%uhist%get_last_state_copy (copy)
  end subroutine get_last_state_copy

  function last_time (this) result (t)
    class(idaesol), intent(in) :: this
    real(r8) :: t
    t = this%uhist%last_time()
  end function last_time

  function last_step_size (this) result (h)
    class(idaesol), intent(in) :: this
    real(r8) :: h
    h = this%hlast
  end function last_step_size

  subroutine get_interpolated_state (this, t, u, first)
    class(idaesol), intent(in) :: this
    real(r8), intent(in) :: t
    real(r8), intent(out) :: u(:)
    integer, intent(in), optional :: first
    call this%uhist%interp_state (t, u, first)
  end subroutine get_interpolated_state

  subroutine write_metrics_unit (this, unit)
    class(idaesol), intent(in) :: this
    integer, intent(in) :: unit
    character(80) :: string(2)
    call write_metrics_string (this, string)
    write(unit,'(/,(a))') string
  end subroutine write_metrics_unit

  subroutine write_metrics_string (this, string)
    class(idaesol), intent(in) :: this
    character(*), intent(out) :: string(:)
    ASSERT(size(string) == 2)
    write(string(1),fmt='(a,i6,a,es11.5,a,es9.3)') &
      'STEP=', this%seq, ', T=', this%uhist%last_time(), ', H=', this%hlast
    write(string(2),fmt='(a,i7.7,":",i5.5,a,5(i4.4,:,":"))') &
      'NFUN:NPC=', this%pcfun_calls, this%updpc_calls, &
      ', NPCF:NNR:NNF:NSR=', this%updpc_failed, &
      this%retried_bce, this%failed_bce, this%rejected_steps
  end subroutine write_metrics_string

  subroutine set_verbose_stepping (this, unit)
    class(idaesol), intent(inout) :: this
    integer, intent(in) :: unit
    this%unit = unit
    this%verbose = .true.
  end subroutine set_verbose_stepping

  subroutine set_quiet_stepping (this)
    class(idaesol), intent(inout) :: this
    this%verbose = .false.
  end subroutine set_quiet_stepping

  subroutine init (this, model, params)

    use parameter_list_type
    use logging_services

    class(idaesol), intent(out) :: this
    !class(idaesol_model), pointer, intent(in) :: model
    class(idaesol_model), intent(in), target :: model
    type(parameter_list) :: params

    integer :: stat
    character(:), allocatable :: context, errmsg
    integer :: maxv
    real(r8) :: vtol

    this%model => model

    this%n = model%size()
    INSIST(this%n > 0)

    context = 'processing ' // params%name() // ': '
    call params%get ('nlk-max-iter', this%mitr, default=5, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (this%mitr < 2) call LS_fatal (context//'"nlk-max-iter" must be > 1')

    call params%get ('nlk-tol', this%ntol, default=0.1_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (this%ntol <= 0.0_r8 .or. this%ntol > 1.0_r8) &
        call LS_fatal (context//'"nlk-tol" must be > 0.0 and <= 1.0')

    call params%get ('nlk-max-vec', maxv, default=this%mitr-1, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (maxv <= 0) call LS_fatal (context//'"nlk-max-vec" must be > 0')
    maxv = min(maxv, this%mitr-1)

    call params%get ('nlk-vec-tol', vtol, default=0.01_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (vtol <= 0.0_r8) call LS_fatal (context//'"nlk-vec-tol" must be > 0.0')

    !! Initialize the NKA structure.
    call this%nka%init (this%n, maxv)
    call this%nka%set_vec_tol (vtol)

    !call this%nka%set_dot_prod (pardp)

    !! We need to maintain 3 solution vectors for quadratic extrapolation.
    call this%uhist%init (3, this%n)

  end subroutine init

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SET_INITIAL_STATE
 !!

  subroutine set_initial_state (this, t, u, udot)
    class(idaesol), intent(inout) :: this
    real(r8), intent(in) :: t, u(:), udot(:)
    ASSERT(size(u) == this%n)
    ASSERT(size(udot) == this%n)
    call this%uhist%flush (t, u, udot)
    this%seq  = 0
  end subroutine set_initial_state

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  BDF2_STEP_DRIVER
 !!

  subroutine bdf2_step_driver (this, hnext, status, &
                               nstep, tout, hmin, hmax, mtry)

    class(idaesol), intent(inout) :: this
    real(r8), intent(inout) :: hnext
    integer,  intent(out)   :: status
    integer,  intent(in) :: nstep, mtry
    real(r8), intent(in) :: tout, hmin, hmax
    optional :: nstep, tout, hmin, hmax, mtry

    integer  :: max_step, max_try, step, errc
    real(r8) :: tlast, h, t_out, h_min, h_max, u(this%n)

    ASSERT( this%seq >= 0 )

   !!!
   !!! PROCESS THE INPUT AND CHECK IT FOR CORRECTNESS

    status = 0

    !! Set the maximum number of time steps; default is unlimited.
    max_step = huge(1)
    if (present(nstep)) max_step = nstep
    if (max_step < 1) status = BAD_INPUT

    !! Set the target integration time; default is +infinity.
    t_out = huge(1.0_r8)
    if (present(tout)) t_out = tout
    if (t_out <= this%uhist%last_time()) status = BAD_INPUT

    !! Verify that at least one of NSTEP and TOUT were specified.
    if (.not.present(nstep) .and. .not.present(tout)) status = BAD_INPUT

    !! Set the bounds on the step size; default is none.
    h_min = tiny(1.0_r8)
    h_max = huge(1.0_r8)
    if (present(hmin)) h_min = hmin
    if (present(hmax)) h_max = hmax
    if (h_min < 0.0_r8 .or. hnext < h_min .or. hnext > h_max) status = BAD_INPUT

    !! Set the maximum number of attempts allowed for a step.
    max_try = 10
    if (present(mtry)) max_try = mtry
    if (max_try < 1) status = BAD_INPUT

    if (status == BAD_INPUT) return

   !!!
   !!! BEGIN TIME STEPPING

    step = 0
    STEP_LOOP: do

      h = hnext
      tlast = this%uhist%last_time()

      !! Check for a normal return before proceeding.
      if (t_out <= tlast) then
        status = SOLVED_TO_TOUT
        exit STEP_LOOP
      end if
      if (step >= max_step) then
        status = SOLVED_TO_NSTEP
        exit STEP_LOOP
      end if

      step = step + 1

      call bdf2_step (this, h, h_min, max_try, u, hnext, errc)
      if (errc /= 0) then
        status = STEP_FAILED
        exit STEP_LOOP   ! step failed
      end if

      !! BDF2 step was successful; commit the solution.
      call commit_state (this, tlast+h, u)

      !! Set the next step size.
      hnext = min(h_max, hnext)
      if (this%verbose) write(this%unit,fmt=1) hnext/h

    end do STEP_LOOP

    1 format(2x,'Changing H by a factor of ',f6.3)

  end subroutine bdf2_step_driver

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! COMMIT_SOLUTION
 !!

  subroutine commit_state (this, t, u)

    class(idaesol), intent(inout) :: this
    real(r8), intent(in) :: t, u(:)

    real(r8) :: h

    h = t - this%uhist%last_time()
    ASSERT(h > 0.0_r8)
    call this%uhist%record_state (t, u)

    this%hlast = h
    this%seq = this%seq + 1
    this%freeze_count = max(0, this%freeze_count - 1)

    this%hmin = min(h, this%hmin)
    this%hmax = max(h, this%hmax)

  end subroutine commit_state

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! BDF2_STEP
 !!

  subroutine bdf2_step (this, h, hmin, mtry, u, hnext, errc)

    type(idaesol), intent(inout) :: this
    real(r8), intent(inout) :: h
    real(r8), intent(in)    :: hmin
    integer,  intent(in)    :: mtry
    real(r8), intent(out)   :: u(:), hnext
    integer,  intent(out)   :: errc

    integer :: try

    errc = 0
    if (hmin < 0.0_r8 .or. h < hmin) errc = BAD_INPUT
    if (mtry < 1) errc = BAD_INPUT
    if (errc == BAD_INPUT) return

    try = 0
    do
      try = try + 1

      !! Check for too many attempts at a single step.
      if (try > mtry) then
        errc  = STEP_FAILED
        exit
      end if

      !! Check for a too-small step size.
      if (h < hmin) then
        errc  = STEP_SIZE_TOO_SMALL
        exit
      end if

      !! Attempt a BDF2 step.
      call step (this, this%uhist%last_time()+h, u, hnext, errc)
      !call step (this, h, u, hnext, errc)
      if (errc == 0) exit

      !! Step failed; try again with the suggested step size.
      if (this%verbose) write(this%unit,fmt=1) hnext/h
      h = hnext

    end do

    1 format(2x,'Changing H by a factor of ',f6.3)

  end subroutine bdf2_step

  subroutine step (this, t, u, hnext, errc)

    class(idaesol), intent(inout) :: this
    real(r8), intent(in)  :: t
    real(r8), intent(out) :: u(:)
    real(r8), intent(out) :: hnext
    integer,  intent(out) :: errc

    real(r8) :: eta, etah, h, t0, tlast, perr, dt(3)
    real(r8) :: u0(size(u)), up(size(u))
    logical  :: fresh_pc, predictor_error

    ASSERT(this%seq >= 0)
    ASSERT(size(u) == this%n)

    tlast = this%uhist%last_time()
    h = t - tlast
    INSIST(h > 0)

    if (this%verbose) write(this%unit,fmt=1) this%seq+1, tlast, h

    !! Predicted solution and base point for BCE step.
    if (this%uhist%depth() == 2) then ! BDF1
      etah = h
      t0 = tlast
      call this%uhist%interp_state (t,  up, order=1)
      call this%uhist%interp_state (t0, u0, order=0)
    else  ! BDF2
      eta = (this%hlast + h) / (this%hlast + 2.0_r8 * h)
      etah = eta * h
      t0 = tlast + (1.0_r8 - eta)*h
      call this%uhist%interp_state (t,  up, order=2)
      call this%uhist%interp_state (t0, u0, order=1)
    end if

    !! Check the predicted solution for admissibility.
    call this%model%schk (up, 0, errc)
    if (errc /= 0) then ! it's bad; cut h and retry.
      this%rejected_steps = this%rejected_steps + 1
      if (this%verbose) write(this%unit,fmt=7)
      hnext = 0.25_r8 * h
      this%freeze_count = 1
      this%usable_pc = .false.
      errc = 3
      return
    end if

    fresh_pc = .false.

    !! If the PC step size is too different than the current step size we tag
    !! it as unusable in order to preempt a possible nonlinear solve failure.
    if (this%usable_pc) then
      if (this%hpc/etah > 1.0_r8 + MARGIN) this%usable_pc = .false.
      if (etah/this%hpc > 1.0_r8 + MARGIN) this%usable_pc = .false.
    end if

    BCE: do

      !! Update the preconditioner if necessary.
      if (.not.this%usable_pc) then
        this%updpc_calls = this%updpc_calls + 1
        call this%model%compute_precon (t, up, etah)
        if (this%verbose) write(this%unit,fmt=3) t
        this%hpc = etah
        this%usable_pc = .true.
        fresh_pc = .true.
      end if

      !! Solve the nonlinear BCE system.
      u = up ! Initial solution guess is the predictor.
      call bce_step (this, t, etah, u0, u, errc)
      if (errc == 0) exit BCE ! the BCE step was successful.

      if (fresh_pc) then ! preconditioner was fresh; cut h and return error condition.
        if (errc == 2) then ! inadmissible iterate generated
          this%rejected_steps = this%rejected_steps + 1
          hnext = 0.25_r8 * h
          this%freeze_count = 1
          this%usable_pc = .false.
        else
          this%failed_bce = this%failed_bce + 1
          hnext = 0.5_r8 * h
          this%freeze_count = 2
        end if
        errc = 1
        return
      else ! update the preconditioner and retry the nonlinear solve.
        this%retried_bce = this%retried_bce + 1
        this%usable_pc = .false.
        cycle BCE
      end if

    end do BCE

    predictor_error = (this%seq >= 3)

    if (predictor_error) then

      !! Predictor error control.
      u0 = u - up
      call this%model%du_norm (u, u0, perr)
      if (perr < 4.0_r8) then ! accept the step.
        if (this%verbose) write(this%unit,fmt=4) perr
        errc = 0
      else ! reject the step; cut h and return error condition.
        this%rejected_steps = this%rejected_steps + 1
        if (this%verbose) write(this%unit,fmt=5) perr
        hnext = 0.5_r8 * h
        this%freeze_count = 1
        errc = 2
        return
      end if

      !! Select the next step size based on the predictor error and past step
      !! size history, but don't change the step size by too great a factor.
      dt(1) = h
      dt(2:) = h + this%uhist%time_deltas()
      call select_step_size (dt, perr, hnext)
      hnext = max(RMIN*h, min(RMAX*h, hnext))
      if (this%freeze_count /= 0) hnext = min(h, hnext)

    else

      if (this%verbose) write(this%unit,fmt=6)
      hnext = h
      errc = 0

    end if

    1 format(/,'BDF2 step ',i6,': T=',es12.5,', H=',es12.5)
    3 format(2x,'Preconditioner updated at T=',es12.5)
    4 format(2x,'Step accepted: perr=',es12.5)
    5 format(2x,'Step REJECTED: perr=',es12.5)
    6 format(2x,'Step accepted: no local error control')
    7 format(2x,'Step REJECTED: inadmissible predicted solution')

  end subroutine step

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !!  SELECT_STEP_SIZE -- Choose a new time step.
 !!

  subroutine select_step_size (dt, perr, h)

    real(r8), intent(in)  :: dt(:), perr
    real(r8), intent(out) :: h

    real(r8), parameter :: tol = 0.001_r8
    real(r8) :: a, dh, phi, dphi

    ASSERT( size(dt) == 3 )

    a = 0.5_r8*dt(1)*dt(2)*dt(3)/max(perr,0.001_r8)
    h = dt(1)

    do ! until converged -- DANGEROUS!

      phi  = h*(h + dt(1))*(h + dt(2)) - a
      dphi = (2.0_r8*h + dt(1))*(h + dt(2)) + h*(h + dt(1))

      dh = phi / dphi
      h = h - dh
      if (abs(dh) / h < tol) exit

    end do

  end subroutine select_step_size

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!
 !! SOLVE_BCE_AIN -- Solve the Backward Cauchy-Euler system using AIN.
 !!
 !! The backward Cauchy-Euler (BCE) method applied to the implicit DAE
 !!
 !!     f(t,u,u') = 0
 !!
 !! yields a nonlinear system of equations for advancing the solution from a
 !! given state u0 at time t - h to the unknown solution u at time t,
 !!
 !!     f(t,u,(u-u0)/h) = 0.
 !!
 !! This subroutine solves this nonlinear system using an accelerated fixed
 !! point iteration [1] for the preconditioned system
 !! g(u) = pc(f(t,u,(u-u0)/h)) = 0:
 !!
 !!    u given
 !!    Do until converged:
 !!      du <-- g(u)
 !!      du <-- NKA(du)
 !!      u  <-- u - du
 !!    End do
 !!
 !! The procedure NKA uses information about g' gleaned from the unaccelerated
 !! correction du=g(u) and previous g values to compute an improved correction.
 !! The preconditioning function pc() is typically an approximate solution
 !! of the Newton correction equation  J*du = f(t,u,(u-u0)/h) where J is an
 !! approximation to the Jacobian of f(t,u,(u-u0)/h) as a function of u.  Thus
 !! this method can be regarded as an accelerated inexact Newton (AIN) method.
 !!
 !! The dummy procedure PCFUN evaluates the preconditioned function g.
 !!
 !! [1] N.N.Carlson and K.Miller, "Design and application of a gradient-
 !!     weighted moving finite element code I: in one dimension", SIAM J.
 !!     Sci. Comput;, 19 (1998), pp. 728-765..
 !!

  subroutine bce_step (this, t, h, u0, u, errc)

    type(idaesol), intent(inout) :: this
    real(r8), intent(in)    :: t, h, u0(:)
    real(r8), intent(inout) :: u(:)
    integer,  intent(out)   :: errc

    integer  :: itr
    real(r8) :: error, du(size(u))

    call this%nka%restart

    itr = 0
    do

      if (itr >= this%mitr) then  ! too many nonlinear iterations
        if (this%verbose) write(this%unit,fmt=1) itr, error
        errc = 1
        exit
      end if

      itr = itr + 1

      !! Evaluate the nonlinear function and precondition it.
      this%pcfun_calls = this%pcfun_calls + 1
      call this%model%compute_f (t, u, (u-u0)/h, du)
      call this%model%apply_precon (t, u, du)

      !! NKA accelerated correction.
      call this%nka%accel_update (du)

      !! Next solution iterate.
      u  = u - du

      !! Check the solution iterate for admissibility.
      call this%model%schk (u, 1, errc)
      if (errc /= 0) then ! iterate is bad; bail.
        if (this%verbose) write(this%unit,fmt=4) itr
        errc = 2
        exit
      end if

      !! Error estimate.
      call this%model%du_norm(u, du, error)
      if (this%verbose) write(this%unit,fmt=3) itr, error

      !! Check for convergence.
      if (((error < this%ntol) .and. (itr > 1)) .or. (error < 0.01_r8 * this%ntol)) then
        if (this%verbose) write(this%unit,fmt=2) itr, error
        errc = 0
        exit
      end if

    end do

    1 format(2x,'NLK BCE solve FAILED: ',i3,' iterations (max), error=',es12.5)
    2 format(2x,'NLK BCE solve succeeded: ',i3,' iterations, error=',es12.5)
    3 format(2x,i3,': error=',es12.5)
    4 format(2x,'NLK BCE solve FAILED: inadmissible solution iterate: itr=',i3)

  end subroutine bce_step

end module idaesol_type
