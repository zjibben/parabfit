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

    this%fluidRho = 1.0_r8 ! just set this to 1 everywhere for now
    
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
  !   use fluid_data_module,      only: fluid_flow, fluidRho, Solid_Face, fluid_to_move, Fluxing_Velocity
  !   use parameter_module,       only: ncells, nfc, ndim
  !   use predictor_module,       only: predictor
  !   use projection_data_module, only: Face_Density, mac_projection_iterations, prelim_projection_iterations
  !   use projection_module,      only: PROJECTION
  !   use property_module,        only: FLUID_PROPERTIES
  !   use time_step_module,       only: cycle_number
  !   use viscous_data_module,    only: prelim_viscous_iterations, viscous_iterations
  !   use zone_module,            only: Zone

  !   real(r8), intent(in) :: t
    
  !   real(r8) :: fluidDeltaRho(ncells)    
  !   logical :: solid_face(nfc,ncells), isPureImmobile(ncells)
  !   integer :: status


  !   if (.not.fluid_flow) return

  !   allocate (Face_Density(nfc,ncells), STAT = status)
    
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
    
  !   deallocate (Face_Density)

  ! end subroutine step
  
end module NS_solver_type
