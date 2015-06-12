!!
!! FLOW_SOLVER_TYPE
!!
!! This module defines a class that encapsulates a (time) solver for the
!! discrete heat conduction model.  It is based on IDAESOL which is an
!! integrator for index-1 DAE systems that uses BDF2 time discretization.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

#include "f90_assert.fpp"
!#define DEBUG

module vof_solver_type

  use kinds, only: r8
  use unstr_mesh_type
  use logging_services
  use vof_tools, only: cell_materials
  implicit none
  private

  
  type, public :: vof_solver
    private
    type(unstr_mesh), pointer :: mesh => null() ! reference only -- do not own
     
    !! Pending/current state
    real(r8)                      :: t, dt
    integer, public               :: nmat                       ! number of materials present globally
    real(r8), pointer             :: velocity_cc(:,:) => null() ! cell centered velocity. reference only -- do not own
    integer, dimension(:), allocatable, public :: matl_id
    type(cell_materials), public, pointer :: cell_matls(:) => null()    ! potentially a target
  contains
    procedure :: init
    procedure :: set_initial_state
    ! procedure :: test_initial_state
    ! procedure :: advect_mass
    ! procedure, private :: advect_volume
    ! procedure :: time
    ! procedure :: get_interpolated_solution
    ! procedure :: get_solution_view
    ! procedure :: get_solution_copy
    ! procedure :: write_metrics
    !procedure :: advance_state
    !procedure :: commit_pending_state
    final :: vof_solver_delete
 end type vof_solver
 
contains

  subroutine vof_solver_delete (this)
    type(vof_solver), intent(inout) :: this
    integer :: i
    if (associated(this%cell_matls)) then
       do i = 1,size(this%cell_matls)
          if (allocated(this%cell_matls(i)%matl)) deallocate(this%cell_matls(i)%matl)
       end do
       deallocate(this%cell_matls)
    end if
  end subroutine vof_solver_delete

  subroutine init (this, mesh, velocity_cc, params)
    
    use parameter_list_type
    use vof_init

    class(vof_solver), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    real(r8), dimension(:,:), intent(in), target :: velocity_cc
    type(parameter_list) :: params

    type(parameter_list), pointer :: plist
    character(:), allocatable :: context

    this%velocity_cc => velocity_cc
    this%mesh        => mesh

    allocate(this%cell_matls(this%mesh%ncell))

  end subroutine init

  subroutine set_initial_state (this, plist)
    use parameter_list_type
    use vof_init
    use vof_tools, only: distinct_matls, volume_of_matl
    
    class(vof_solver), intent(inout) :: this
    type(parameter_list), intent(in) :: plist
    ! real(r8), intent(in) :: t, temp(:), dt

    ! integer :: stat
    ! character(:), allocatable :: errmsg
    ! real(r8), allocatable :: udot(:)
    
    !! Initialize the Vof
    call vof_initialize (this%mesh, plist, this%cell_matls)
    call distinct_matls(this%nmat, this%matl_id, this%cell_matls)

    write(*,*) volume_of_matl(1, this%cell_matls, this%mesh%volume)
    write(*,*) volume_of_matl(2, this%cell_matls, this%mesh%volume)
    
    ! allocate(udot(size(this%u)))
    ! call compute_initial_state (this%model, t, temp, dt, this%u, udot, stat, errmsg)
    ! if (stat /= 0) call LS_fatal ('HC_SOLVER%SET_INITIAL_STATE: ' // errmsg)
    ! call this%integ%set_initial_state (t, this%u, udot)

  end subroutine set_initial_state

  ! subroutine test_initial_state (this, t, temp, dt, u, udot)

  !   class(HC_solver), intent(inout) :: this
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
  !   class(HC_solver), intent(in) :: this
  !   time = this%integ%last_time()
  ! end function time

  ! !! Returns the solution U interpolated to time T.  This should only
  ! !! be called when the integrator has first stepped across time T, so
  ! !! that T lies within an interval of very recent time steps where
  ! !! solution data is currently available.
  ! subroutine get_interpolated_solution (this, t, u)
  !   class(HC_solver), intent(in) :: this
  !   real(r8), intent(in)  :: t
  !   real(r8), intent(out) :: u(:)
  !   ASSERT(size(u) == this%model%num_dof())
  !   call this%integ%get_interpolated_state (t, u)
  ! end subroutine get_interpolated_solution

  ! subroutine get_solution_view (this, u)
  !   class(HC_solver), intent(in) :: this
  !   real(r8), pointer, intent(out) :: u(:)
  !   call this%integ%get_last_state_view (u)
  ! end subroutine get_solution_view

  ! subroutine get_solution_copy (this, u)
  !   class(HC_solver), intent(in) :: this
  !   real(r8), intent(out) :: u(:)
  !   call this%integ%get_last_state_copy (u)
  ! end subroutine get_solution_copy

  ! subroutine write_metrics (this, string)
  !   class(HC_solver), intent(in) :: this
  !   character(*), intent(out) :: string(:)
  !   ASSERT(size(string) == 2)
  !   call this%integ%write_metrics (string)
  ! end subroutine write_metrics

  
  
  !   subroutine advect_mass(this)
  !   integer :: status
  !   real(r8), dimension(:,:), allocatable :: Vof, Vof_n

    
  !   ! Allocate working arrays
  !   allocate (Mask(ncells),                     &
  !             Tmp(ncells),                      &
  !             Vof  (this%nmat,this%mesh%ncell), &
  !             Vof_n(this%nmat,this%mesh%ncell), STAT = status)
  !   if (status /= 0) call TLS_panic ('ADVECT_MASS: allocation failed')
    
  !   if (volume_track_interfaces .and. .not. ASSOCIATED(VT_Interface_Mask)) then
  !      ALLOCATE (VT_Interface_Mask(nfc,ncells), STAT = status)
  !      if (status /= 0) call TLS_panic ('ADVECT_MASS: VT_Interface_Mask(nfc,ncells) allocation failed')
  !   end if
    
  !   if (.not. ALLOCATED(Volume_Flux)) then
  !      ALLOCATE (Volume_Flux(nmat,nfc,ncells), STAT=status)
  !      if (status /= 0) call TLS_panic ('ADVECT_MASS: Volume_Flux(nmat,nfc,ncells) allocation failed')
  !   end if
    
  !   ! Grab a copy of the volume fractions from Matl; will return time n+1 values to Matl at the end
  !   ! of this subroutine, via MATL_SET_VOF
  !   call matl_get_vof (vof, this%cell_matls)

  !   ! Will also need a copy of the time n Vof values.
  !   Vof_n = Vof

  !   ! Advect material volumes.
  !   call advect_volume (Vof, Vof_n, Fluxing_Velocity, Volume_Flux)

  !   ! Update the mass & concentration distributions.
  !   if (volume_track_interfaces) then 
  !      call UPDATE_MASS (Fluxing_Velocity, Vof, Vof_n, VT_Interface_Mask)
  !   else
  !      call UPDATE_MASS (Fluxing_Velocity, Vof, Vof_n)
  !   end if

  !   ! Return Vof values at this point back into the Matl structure.
  !   call MATL_SET_VOF (Vof)

  !   ! Accumulate Momentum Advection Array
  !   call ADVECT_MOMENTUM_ACCUMULATION ()

  !   ! Find and store the accumulated inflow and outflow mass.
  !   do f = 1, nfc
  !      Mask = IN_FLOW (f, Fluxing_Velocity) .OR. OUT_FLOW(f, Fluxing_Velocity)

  !      if (PGSLib_Global_ANY(Mask)) then
  !         ! Accumulate inflow volume - Fluxing_Velocity < zero
  !         Tmp = MIN(dt*Cell%Face_Area(f)*Fluxing_Velocity(f,:),0.0_r8)
  !         qin = qin + ABS(PGSLib_Global_SUM(Tmp, MASK = Mask))

  !         ! Accumulate outflow volume - Fluxing_Velocity > zero
  !         Tmp = MAX(dt*Cell%Face_Area(f)*Fluxing_Velocity(f,:),0.0_r8)
  !         qout = qout + PGSLib_Global_SUM(Tmp, MASK = Mask)

  !      end if
  !   end do

  !   deallocate (Mask)
  !   deallocate (Tmp)
  !   deallocate (Vof)
  !   deallocate (Vof_n)

  !   ! Stop the Volume Advection Timer
  !   call stop_timer("Mass Advection")
    
  ! end subroutine advect_mass

!   subroutine advect_volume (this, vof, vof_n, fluxing_velocity, volume_flux)
!     real(r8), dimension(this%nmat,ncells),     intent(INOUT) :: Vof
!     real(r8), dimension(this%nmat,ncells),     intent(IN)    :: Vof_n
!     real(r8), dimension(nfc,ncells),      intent(IN)    :: Fluxing_Velocity
!     real(r8), dimension(nmat,nfc,ncells), intent(OUT)   :: Volume_Flux_Tot

!     ! Local Variables
!     integer :: status
!     integer :: p, vps
!     real(r8), dimension(:,:,:), allocatable :: Volume_Flux_Sub
    
    
!     ! Start the volume advection timer.
!     call start_timer("Volume Tracking")

!     ! Zero the total flux array.
!     Volume_Flux_Tot = 0.0_r8

!     ! No subcycling if we're not volume tracking.
!     vps = volume_track_subcycles
    
!     ! Set the advection timestep.
!     adv_dt = dt/vps

!     ! If we're volume tracking, we'll need an array to keep track of volume
!     ! changes in every subcycle.
!     ALLOCATE (Volume_Flux_Sub(nmat,nfc,ncells), STAT = status)
!     if (status /= 0) call TLS_panic ('ADVECT_VOLUME: Volume_Flux_Sub(nmat,nfc,ncells) allocation failed')
    
!     FLUXING_PASSES: do p = 1,vps
!        ! Initialize the array that'll keep track of volume fluxes in this subcycle.
!        Volume_Flux_Sub = 0.0_r8

!        ! Get the donor fluxes.
!        call VOLUME_TRACK (Vof, Fluxing_Velocity, Volume_Flux_Sub)

!        ! Normalize the donor fluxes.
!        call FLUX_RENORM (Fluxing_Velocity, Vof_n, Volume_Flux_Tot, Volume_Flux_Sub)

!        ! Compute the acceptor fluxes.
!        call FLUX_ACCEPTOR (Volume_Flux_Sub)

!        ! Compute BC (inflow) fluxes.
!        call FLUX_BC (Fluxing_Velocity, Vof_n, Volume_Flux_Sub)

!        ! Add the volume fluxes from this subcycle (Volume_Flux_Sub) to the
!        ! total flux array (Volume_Flux_Tot), and update the volume fraction
!        ! array (Vof).
!        call VOLUME_ADVANCE (Volume_Flux_Sub, Volume_Flux_Tot, Vof)

!        ! Make sure volume fractions of a particular material are within
!        ! the allowed range (0 <= Vof <= 1) and that all materials sum to one.
!        call VOF_BOUNDS (Vof, Volume_Flux_Tot)
!     end do FLUXING_PASSES

!     DEALLOCATE (Volume_Flux_Sub)

!     ! Stop the volume advection timer.
!     call stop_timer("Volume Tracking")
    
!   end subroutine advect_volume

  
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
  
end module vof_solver_type
