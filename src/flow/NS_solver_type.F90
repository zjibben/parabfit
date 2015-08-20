!!
!! FLOW_SOLVER_TYPE
!!
!! This module defines a class that encapsulates a (time) solver for the
!! Navier-Stokes equations. 
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

#include "f90_assert.fpp"
!#define DEBUG

module NS_solver_type

  use kinds, only: r8
  use unstr_mesh_type
  use logging_services
  implicit none
  private
  
  type, public :: NS_solver
     private
     !type(HC_model), pointer :: model => null() ! reference only -- do not own
     type(unstr_mesh), pointer :: mesh => null() ! reference only -- do not own
     
     !! Pending/current state
     real(r8) :: t, dt
     real(r8), public, allocatable :: velocity(:,:) ! potentially a target
     real(r8), public, allocatable :: fluidRho(:)
     real(r8), pointer :: vof(:,:) ! reference only -- do not own
     logical, public   :: use_prescribed_velocity
     integer, public   :: prescribed_velocity_case
   contains
     procedure :: init
     procedure :: init_matls
     procedure :: set_initial_state
     ! procedure :: test_initial_state
     ! procedure :: time
     ! procedure :: get_interpolated_solution
     ! procedure :: get_solution_view
     ! procedure :: get_solution_copy
     ! procedure :: write_metrics
     !procedure :: advance_state
     !procedure :: commit_pending_state
     final :: NS_solver_delete
  end type NS_solver

contains
  
  subroutine NS_solver_delete (this)
    type(NS_solver), intent(inout) :: this
    ! if (associated(this%precon)) deallocate(this%precon)
    ! if (associated(this%norm)) deallocate(this%norm)
    ! if (associated(this%integ_model)) deallocate(this%integ_model)
    !if (associated(this%velocity)) deallocate(this%velocity)
  end subroutine NS_solver_delete
  
  subroutine init (this, mesh, params)
    
    use parameter_list_type
    
    class(NS_solver), intent(out) :: this
    !type(HC_model), intent(in), target :: model
    type(unstr_mesh), intent(in), target :: mesh
    type(parameter_list) :: params
    
    type(parameter_list), pointer :: plist
    character(:), allocatable :: context,errmsg
    integer :: stat
    ! this%model => model

    this%mesh => mesh
    allocate(this%velocity(3,this%mesh%ncell), this%fluidRho(this%mesh%ncell))
    ! allocate(this%u(this%model%num_dof()))

    
    !! check for prescribed velocity case
    context = 'processing ' // params%name() // ': '
    this%use_prescribed_velocity = params%is_scalar('prescribed-velocity')
    if (this%use_prescribed_velocity) then
       call params%get ('prescribed-velocity', this%prescribed_velocity_case, stat=stat, errmsg=errmsg)
       if (stat /= 0) call LS_fatal (context//errmsg)
    end if
    
    ! !! Create the preconditioner
    ! context = 'processing ' // params%name() // ': '
    ! if (params%is_sublist('preconditioner')) then
    !   plist => params%sublist('preconditioner')
    !   allocate(this%precon)
    !   call this%precon%init (this%model, plist)
    ! else
    !   call LS_fatal (context//'missing "preconditioner" sublist parameter')
    ! end if

    ! !! Create the error norm
    ! if (params%is_sublist('error-norm')) then
    !   allocate(this%norm)
    !   plist => params%sublist('error-norm')
    !   call this%norm%init (this%model, plist)
    ! else
    !   call LS_fatal (context//'missing "error-norm" sublist parameter')
    ! end if
    
    ! !! Create the IDAESOL model
    ! allocate(this%integ_model)
    ! call this%integ_model%init (this%model, this%precon, this%norm)

    ! !! Create the IDAESOL integrator
    ! if (params%is_sublist('integrator')) then
    !   plist => params%sublist('integrator')
    !   call this%integ%init (this%integ_model, plist)
    ! else
    !   call LS_fatal (context//'missing "integrator" sublist parameter')
    ! end if
    
  end subroutine init

  subroutine init_matls(this, vof)
    class(NS_solver), intent(inout) :: this
    real(r8), dimension(:,:), intent(in), target :: vof
    
    this%vof => vof
  end subroutine init_matls
  
  subroutine set_initial_state (this) !, t, temp, dt)
    use prescribed_velocity_fields, only: prescribed_velocity
    
    class(NS_solver), intent(inout) :: this
    !real(r8), intent(in) :: t, temp(:), dt

    integer :: i
    ! integer :: stat
    ! character(:), allocatable :: errmsg
    !real(r8), allocatable :: udot(:)

    ! allocate(udot(size(this%u)))
    ! call compute_initial_state (this%model, t, temp, dt, this%u, udot, stat, errmsg)
    ! if (stat /= 0) call LS_fatal ('HC_SOLVER%SET_INITIAL_STATE: ' // errmsg)
    ! call this%integ%set_initial_state (t, this%u, udot)

    do i = 1,this%mesh%ncell
       this%velocity(:,i) = prescribed_velocity (this%mesh%x(:,i), 0.0_r8, this%prescribed_velocity_case)
    end do
    
  end subroutine set_initial_state
  
  ! !=======================================================================
  ! ! Purpose(s):
  ! !
  ! !   Navier-Stokes (NS) driver: increment NS equations by one time step.
  ! !
  ! !======================================================================
  ! subroutine step (t)
  !   use fluid_data_module,      only: fluid_flow, fluidRho, Solid_Face, &
  !        isPureImmobile, fluidDeltaRho,      &
  !        fluid_to_move, Fluxing_Velocity
  !   use parameter_module,       only: ncells, nfc, ndim
  !   use predictor_module,       only: PREDICTOR
  !   use projection_data_module, only: Face_Density, mac_projection_iterations, &
  !        prelim_projection_iterations
  !   use projection_module,      only: PROJECTION
  !   use property_module,        only: FLUID_PROPERTIES
  !   use time_step_module,       only: cycle_number
  !   use viscous_data_module,    only: prelim_viscous_iterations, &
  !        viscous_iterations
  !   use zone_module,            only: Zone

  !   real(r8), intent(in) :: t
    
  !   integer :: status

  !   if (.not.fluid_flow) return

  !   allocate (Solid_Face(nfc,ncells), isPureImmobile(ncells), &
  !        fluidDeltaRho(ncells), Face_Density(nfc,ncells), STAT = status)
  !   if (status /= 0) call LS_fatal ('Navier-Stokes step: allocation failed')
    
  !   ! Evaluate cell properties excluding immobile materials, and
  !   ! check that there are at least some flow equations to solve
  !   fluidRho = 0.0_r8
  !   call fluid_properties (fluid_to_move, t)

  !   if (fluid_to_move) then
  !     call predictor ()  ! Predictor Step
  !     call projection () ! Projection Step
      
  !     if (cycle_number == 0) then
  !       ! Special operations required during the prepass
  !       prelim_projection_iterations = mac_projection_iterations
  !       prelim_viscous_iterations = viscous_iterations
  !       Zone%Vc = Zone%Vc_Old
  !     end if
  !   else
  !     ! Everything solid; set velocities to zero and check again in the next timestep.
  !     Fluxing_Velocity = 0
  !     Zone%Vc = 0
  !     Zone%Vc_Old = 0

  !     if (cycle_number == 0) then
  !       prelim_projection_iterations = 0
  !       prelim_viscous_iterations = 0
  !     end if
  !   end if
    
  !   deallocate (Solid_Face,isPureImmobile,fluidDeltaRho,Face_Density)

  ! end subroutine step
  
  ! !===============================================================================
  ! ! Purpose(s):
  ! !
  ! !   Predictor phase of the incompressible Navier-Stokes solution
  ! !   algorithm. This routine is the main driver for advancing the
  ! !   cell-centered velocity from time "n" to time "*". The "*" time
  ! !   level velocity field is an estimate of the "n+1" velocity field
  ! !   without regard for the solenoidal constraint.
  ! !
  ! !   This routine basically completes the entire RHS for the momentum equations,
  ! !   then calls SOLVE_FOR_VELOCITY to get the predicited velocity, i.e., the *
  ! !   state.
  ! !
  ! !===============================================================================
  ! subroutine predictor (dt)
  !   use advection_module,        only: ADVECT_MOMENTUM
  !   use body_data_module,        only: body_force_face_method
  !   use body_force_module,       only: add_cell_body_force
  !   use fluid_data_module,       only: fluidRho_n,           &
  !        fluidVof, fluidVof_n,           &
  !        Drag_Coefficient,               &
  !        Mom_Delta,                      &
  !        momentum_solidify_implicitness
  !   use flow_phase_change,       only: have_solidifying_flow, solidified_rho
  !   use parameter_module,        only: ncells, ndim
  !   use porous_drag_data,        only: porous_flow
  !   use porous_drag_module,      only: POROUS_DRAG
  !   use timing_tree
  !   use turbulence_module,       only: turbulence_model, TURBULENCE
  !   use viscous_data_module,     only: inviscid, stokes, viscous_implicitness
  !   use viscous_module,          only: viscousExplicit,viscousCleanup
  !   use zone_module,             only: Zone
  !   use surface_tension_module,  only: CSF, csf_tangential

  !   real(r8), intent(in) :: dt

  !   integer :: i, n
  !   real(r8) :: tweight
  !   logical :: solidifying_flow

  !   call start_timer ("Predictor")

  !   ! Initialize
  !   Mom_Delta        = 0
  !   Drag_Coefficient = 0

  !   call predictorSetup (dt, Mom_Delta)

  !   ! Explicit Pressure Gradient and Buoyant Force
  !   call PressureExplicit (dt, Mom_Delta)

  !   ! Explicit Viscous Stress.
  !   if (.not.inviscid) then
  !     call turbulence (turbulence_model)
      
  !     ! Add in the explicit evaluation of the Standard Newtonian stress tensor.
  !     if (viscous_implicitness < 1) call viscousExplicit (dt, Mom_Delta)
  !   end if
    
  !   ! MAC note: Here, there can be problems with the explicit term of the
  !   ! momentum  sink due to solidification.  If the fluid velocity at time
  !   ! level n is zero, there is no contribution to the sink term.  For this
  !   ! reason, the recommended practice, for now, is to use
  !   ! momentum_solidify_implicitness = 1.0.  This is the default.
  !   tweight = 1.0 - momentum_solidify_implicitness
  !   solidifying_flow = have_solidifying_flow ()
    
  !   do i = 1, mesh%ncell
  !     ! Multiply Pressure Gradient and Viscous Stress by FluidVof
  !     ! (This is a crude way to account for solid material within the cell.
  !     ! In SOLVE_FOR_VELOCITY we also divide by FluidVof to account for
  !     ! the mass of fluid in the cell, more or less canceling out this
  !     ! term.  Momentum advection is specifically excluded because the VOF
  !     ! is already accounting for the solid material.)
  !     Mom_Delta(:,i) = FluidVof(i)*Mom_Delta(:,i)
      
  !     ! Advect Momentum... No real advection is done here, instead, the increment from
  !     ! momentum advection is just accumulated into Mom_Delta.
  !     if (.not.stokes) mom_delta(:,i) = mom_delta(:,i) + momentum_delta(:,i)
      
  !     ! calculate the RHS momentum change due to the solidification of fluid
  !     if (solidifying_flow) Mom_Delta(:,i) = Mom_Delta(:,i) - tweight*solidified_rho(i)*Zone(i)%Vc(:)
  !   end do

  !   if (csf_tangential) call CSF (dt, Mom_Delta)      ! Tangential surface tension 
  !   if (porous_flow) call porous_drag (dt, Mom_Delta) ! Porous drag
  !   if (.not.body_force_face_method) call add_cell_body_force (dt, Mom_Delta) ! Add in the simple cell-centered body force

  !   ! Complete the RHS terms in Mom_Delta by including the momentum at time 'n'
  !   do i = 1, ncells
  !     Mom_Delta(:,i) = Mom_Delta(:,i) + fluidRho_n(i)*fluidVof_n(i)*Zone(i)%Vc_old(:)
  !   end do

  !   ! Solve for Star time velocity 
  !   call solve_for_velocity(Mom_Delta)

  !   if (.not.inviscid) call viscousCleanup()

  !   call stop_timer("Predictor")

  ! end subroutine predictor

  ! subroutine test_initial_state (this, t, temp, dt, u, udot)

  !   class(NS_solver), intent(inout) :: this
  !   real(r8), intent(in) :: t, temp(:), dt
  !   real(r8), intent(out) :: u(:), udot(:)

  !   integer :: stat
  !   character(:), allocatable :: errmsg

  !   ASSERT(size(u) == this%model%num_dof())
  !   ASSERT(size(u) == size(udot))

  !   call compute_initial_state (this%model, t, temp, dt, u, udot, stat, errmsg)
  !   if (stat /= 0) call LS_fatal ('HC_SOLVER%SET_INITIAL_STATE: ' // errmsg)

  ! end subroutine test_initial_state

  ! !! Returns the current integration time.
  ! real(r8) function time (this)
  !   class(NS_solver), intent(in) :: this
  !   time = this%integ%last_time()
  ! end function time
  
!   !! Returns the solution U interpolated to time T.  This should only
!   !! be called when the integrator has first stepped across time T, so
!   !! that T lies within an interval of very recent time steps where
!   !! solution data is currently available.
!   subroutine get_interpolated_solution (this, t, u)
!     class(NS_solver), intent(in) :: this
!     real(r8), intent(in)  :: t
!     real(r8), intent(out) :: u(:)
!     ASSERT(size(u) == this%model%num_dof())
!     call this%integ%get_interpolated_state (t, u)
!   end subroutine get_interpolated_solution

!   subroutine get_solution_view (this, u)
!     class(NS_solver), intent(in) :: this
!     real(r8), pointer, intent(out) :: u(:)
!     call this%integ%get_last_state_view (u)
!   end subroutine get_solution_view

!   subroutine get_solution_copy (this, u)
!     class(NS_solver), intent(in) :: this
!     real(r8), intent(out) :: u(:)
!     call this%integ%get_last_state_copy (u)
!   end subroutine get_solution_copy

!   subroutine write_metrics (this, string)
!     class(NS_solver), intent(in) :: this
!     character(*), intent(out) :: string(:)
!     ASSERT(size(string) == 2)
!     call this%integ%write_metrics (string)
!   end subroutine write_metrics
  
!   !! This auxiliary procedure computes the consistent initial state (u, du/dt)
!   !! given the initial cell temperatures.  For a typical explicit ODE system
!   !! du/dt = F(t,u) this is trivial; u is given and F evaluated to get du/dt.
!   !! However for our implicit index-1 DAE system F(t,u,du/dt) = 0 this is much
!   !! more involved.  We are only given part of u; the remaining part must
!   !! obtained by solving the algebraic equation portion of the DAE system.
!   !! Furthermore F=0 only defines du/dt for the cell enthalpies; the remaining
!   !! time derivatives must be solved for (by differentiating F=0 with respect
!   !! to time) or approximated (which we do here).
!   !!
!   !! NB: the time step DT is used to approximate the time derivatives of the
!   !! cell and face temperatures.  The current integration algorithm uses FE
!   !! to get an initial guess for either a BDF1 or trapezoid starting step.
!   !! Consequently, the best choice of DT would be the initial time step, as
!   !! this will give a predicted state that is exactly consistent.
  
!   subroutine compute_initial_state (model, t, temp, dt, u, udot, stat, errmsg)

!     use HC_AE_solver_type
!     use parameter_list_type

!     class(HC_model), intent(inout), target :: model
!     real(r8), intent(in) :: t, temp(:), dt
!     real(r8), intent(out), target :: u(:), udot(:)
!     integer, intent(out) :: stat
!     character(:), allocatable, intent(out) :: errmsg

!     type(HC_AE_solver) :: solver
!     real(r8), allocatable, target :: f(:)
!     real(r8), pointer :: u1(:), u2(:), u3(:), f1(:), f2(:), f3(:), hdot(:)
!     type(parameter_list), pointer :: params

!     ASSERT(size(temp) == model%mesh%ncell)
!     ASSERT(size(u) == model%num_dof())
!     ASSERT(size(udot) == size(u))

!     call model%get_cell_heat_view (u, u1)  ! enthalpy
!     call model%get_cell_temp_view (u, u2)  ! cell temp
!     call model%get_face_temp_view (u, u3)  ! face temp

!     u2 = temp ! set the cell temperatures from the input
!     call model%H_of_T (u2, u1)  ! compute the cell enthalpy

!     !! Solve for the face temperatures.
!     allocate(params)
!     call params%set ('max-iter', 100) !TODO: expose as input
!     call params%set ('rel-tol', 1.0d-6) !TODO: expose as input
!     call solver%init (model, params)
!     u3 = 0.0_r8 ! initial guess (we could do much better)
!     call solver%solve (t, u2, u3, stat, errmsg)
!     if (stat /= 0) then
!       errmsg = 'face temp solve 1: ' // errmsg
!       return
!     end if

!     !! The DAE system F(t,u,udot) = 0 gives the time derivative of the cell
!     !! enthalpy as a a function of the cell and face temperatures.  We back
!     !! out what it is by computing F with udot set equal 0.  The info is
!     !! contained in the cell temperature section of F.  By construction, the
!     !! the remaining sections should be zero (the cell enthalpy section to
!     !! round-off, and the face temperature section to the solver tolerance).

!     allocate(f(size(u)))
!     call model%get_cell_heat_view (f, f1)  ! enthalpy / enthalpy-temp AE
!     call model%get_cell_temp_view (f, f2)  ! cell temp / heat conduction DE
!     call model%get_face_temp_view (f, f3)  ! face temp / face-cell temp AE

!     udot = 0.0_r8
!     call model%residual (t, u, udot, f)

!     call model%get_cell_heat_view (udot, hdot)
!     hdot = -f2 / model%mesh%volume

!     !! The time derivative of the cell and face temperatures are approximated
!     !! by a finite difference.  The enthalpy is advanced by a small time step
!     !! using its time derivative (forward Euler), and then associated advanced
!     !! cell and face temperatures are solved for using the algebraic relations.

!     f1 = u1 + dt*hdot ! advance the enthalpy
!     call model%T_of_H (f1, f2) ! compute cell temperature
!     f3 = u3 ! initial guess (probably not half bad)
!     call solver%solve (t+dt, f2, f3, stat, errmsg) ! compute face temperature
!     if (stat /= 0) then
!       errmsg = 'face temp solve 2: ' // errmsg
!       return
!     end if

!     f2 = (f2 - u2) / dt
!     call model%set_cell_temp (f2, udot)

!     f3 = (f3 - u3) / dt
!     call model%set_face_temp (f3, udot)

!     deallocate(params)

! #ifdef DEBUG
!     call model%residual (t, u, udot, f)
!     print *, '||f1||_max =', maxval(abs(f1))
!     print *, '||f2||_max =', maxval(abs(f2))
!     print *, '||f3||_max =', maxval(abs(f3))
! #endif
!   end subroutine compute_initial_state

end module NS_solver_type
