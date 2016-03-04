!!
!! FLOW_SOLVER_TYPE
!!
!! This module defines a class that encapsulates a (time) solver for the
!! Navier-Stokes equations.
!!
!! TODO: * Figure out what the prepass is and how to do it. Basically,
!!         Truchas executes
!!           use projection_data_module, only: mac_projection_iterations, prelim_projection_iterations
!!           use viscous_data_module,    only: prelim_viscous_iterations, viscous_iterations
!!           if (cycle_number == 0) then
!!             ! Special operations required during the prepass
!!             prelim_projection_iterations = mac_projection_iterations !or 0 if all solid
!!             prelim_viscous_iterations = viscous_iterations           !or 0 if all solid
!!             velocity_cc = velocity_cc_n
!!           end if
!!         on the first cycle.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

#include "f90_assert.fpp"
!#define DEBUG

module NS_solver_type

  use kinds, only: r8
  use parameter_list_type
  use unstr_mesh_type
  use mesh_geom_type
  use matl_props_type
  use logging_services
  use bndry_func_class
  implicit none
  private
  
  type, public :: NS_solver
    private

    type(unstr_mesh), pointer :: mesh => null()  ! reference only -- do not own
    type(mesh_geom),  pointer :: gmesh => null() ! reference only -- do not own
    type(matl_props), pointer :: mprop => null() ! reference only -- do not own

    real(r8), pointer :: vof(:,:)           => null() ! reference only -- do not own
    real(r8), pointer :: volume_flux(:,:,:) => null() ! reference only -- do not own

    !! Pending/current state
    real(r8) :: t, dt
    type(parameter_list), pointer :: projection_solver_params => null()
    
    real(r8), public, allocatable :: velocity_cc(:,:), pressure_cc(:), fluxing_velocity(:,:), &
        fluidRho(:)
    real(r8), allocatable :: gradP_dynamic_cc(:,:), fluidVof(:), body_force(:)
    class(bndry_func), allocatable :: pressure_bc, velocity_bc
    logical         :: boussinesq_approximation
    logical, public :: use_prescribed_velocity
    integer, public :: prescribed_velocity_case
  contains
    procedure :: init
    procedure :: init_matls
    procedure :: init_mprop
    procedure :: set_initial_state
    procedure :: step
    procedure, private :: update_prescribed_velocity
    procedure, private :: fluid_to_move
    procedure, private :: fluid_properties
  end type NS_solver

contains

  subroutine init (this, mesh, gmesh, params)

    use consts, only: ndim,nfc
    use bc_factory_type

    class(NS_solver), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(mesh_geom),  intent(in), target :: gmesh
    type(parameter_list) :: params

    type(parameter_list), pointer :: plist
    character(:), allocatable :: context,errmsg
    integer :: stat
    type(bc_factory) :: bcfact

    this%mesh => mesh
    this%gmesh => gmesh
    allocate(this%velocity_cc(ndim,this%mesh%ncell), this%pressure_cc(this%mesh%ncell), &
        this%fluidRho(this%mesh%ncell), this%fluidVof(this%mesh%ncell), &
        this%gradP_dynamic_cc(ndim,this%mesh%ncell), this%fluxing_velocity(nfc,this%mesh%ncell))

    context = 'processing ' // params%name() // ': '

    !! check for prescribed velocity case
    this%use_prescribed_velocity = params%is_scalar('prescribed-velocity')
    if (this%use_prescribed_velocity) then
      call params%get ('prescribed-velocity', this%prescribed_velocity_case, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      return
    end if
    
    !! for now, ignore buoyancy
    this%boussinesq_approximation = .false.

    !! check for a body force
    if (params%is_scalar('body-force')) then
      call params%get ('body-force', this%body_force, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
    else
      allocate(this%body_force(ndim))
      this%body_force = 0.0_r8
    end if

    !! store the projection solver's hypre parameters
    if (params%is_sublist('projection-solver')) then
      this%projection_solver_params => params%sublist('projection-solver')
    else
      call LS_fatal (context//'missing "projection-solver" sublist parameter')
    end if

    !! initialize boundary conditions
    if (params%is_sublist('bc')) then
      call bcfact%init (this%mesh, params%sublist('bc'))
      call bcfact%alloc_bc ('dirichlet-pressure', this%pressure_bc)
      call bcfact%alloc_bc ('dirichlet-velocity', this%velocity_bc)
    else
      call LS_fatal (context//'missing "bc" sublist parameter')
    end if
    
  end subroutine init

  ! point vof-related quantities to the arrays owned by the vof solver
  subroutine init_matls (this, vof, volume_flux)

    class(NS_solver), intent(inout) :: this
    real(r8), target, intent(in)    :: vof(:,:), volume_flux(:,:,:)

    this%vof => vof
    this%volume_flux => volume_flux

  end subroutine init_matls

  subroutine init_mprop (this, mprop)

    class(NS_solver), intent(inout) :: this
    type(matl_props), target        :: mprop

    this%mprop => mprop

  end subroutine init_mprop

  ! initialize state variables
  subroutine set_initial_state (this)

    class(NS_solver), intent(inout) :: this

    integer :: i

    ! set initial time
    this%t = 0.0_r8

    ! TODO: * set initial values based on user input, especially fluidRho and fluidVof
    !       * maybe thread this?
    
    if (this%use_prescribed_velocity) then
      call this%update_prescribed_velocity ()
    else
      this%velocity_cc = 0.0_r8
      this%fluxing_velocity = 0.0_r8
    end if

    write(*,*) 'WARNING: hard-coded initial condition for single-phase channel flow'
    do i = 1,this%mesh%ncell
      this%pressure_cc(i) = 1.0_r8 - this%gmesh%xc(2,i) / 4.0_r8
      this%gradP_dynamic_cc(:,i) = [0.0_r8, -0.25_r8, 0.0_r8]
      
      ! this%pressure_cc(i) = 0.0_r8
      ! this%gradP_dynamic_cc(:,i) = 0.0_r8
      this%fluidRho(i) = 1.0_r8
      this%fluidVof(i) = 1.0_r8
    end do
    
  end subroutine set_initial_state

  ! sets the velocity to a prescribed value across the entire domain
  subroutine update_prescribed_velocity (this)

    use consts, only: nfc
    use prescribed_velocity_fields, only: prescribed_velocity

    class(NS_solver), intent(inout) :: this

    integer :: i,f
    
    do i = 1,this%mesh%ncell
      this%velocity_cc(:,i) = prescribed_velocity (this%gmesh%xc(:,i), this%t, &
          this%prescribed_velocity_case)

      do f = 1,nfc
        this%fluxing_velocity(f,i) = dot_product( &
            prescribed_velocity (this%gmesh%fc(:,this%mesh%cface(f,i)), this%t, this%prescribed_velocity_case), &
            this%gmesh%outnorm(:,f,i))
      end do
    end do

  end subroutine update_prescribed_velocity

  ! Navier-Stokes (NS) driver: increment NS equations by one time step.
  ! TODO: figure out a smarter way of sending the hypre_hybrid parameters down the pipeline
  subroutine step (this, dt, t)

    use consts, only: ndim
    use predictor_module
    use projection_module

    class(NS_solver), intent(inout) :: this
    real(r8),         intent(in)    :: dt, t

    real(r8) :: fluidRho_n(this%mesh%ncell), fluidVof_n(this%mesh%ncell), rho(this%mesh%ncell), &
        velocity_cc_n(ndim,this%mesh%ncell), min_fluidRho
    logical  :: solid_face(this%mesh%nface), is_pure_immobile(this%mesh%ncell)

    this%t = t

    ! check if we are using a prescribed velocity
    if (this%use_prescribed_velocity) then
      call this%update_prescribed_velocity ()
    else
      call this%pressure_bc%compute (this%t)
      call this%velocity_bc%compute (this%t)

      ! evaluate cell properties excluding immobile materials, and
      ! check that there are at least some flow equations to solve
      call this%fluid_properties (rho, fluidRho_n, fluidVof_n, velocity_cc_n, &
          min_fluidRho, solid_face, is_pure_immobile)
      
      ! step the incompressible Navier-Stokes equations
      if (this%fluid_to_move(is_pure_immobile, solid_face)) then
        ! WARNING: hardcoding inviscid flow for now
        call predictor (this%velocity_cc, this%gradP_dynamic_cc, dt, this%mprop%density, &
            this%volume_flux, this%fluidRho, fluidRho_n, this%fluidVof, fluidVof_n, velocity_cc_n, &
            this%fluxing_velocity, .true., this%mesh, this%gmesh)
        call projection (this%velocity_cc, this%fluxing_velocity, this%pressure_cc, this%gradP_dynamic_cc, &
            dt, this%fluidRho, fluidRho_n, this%fluidVof, this%vof, this%mprop%sound_speed, this%body_force, &
            solid_face, is_pure_immobile, this%mprop%is_immobile, min_fluidRho, &
            this%velocity_bc, this%pressure_bc, this%mesh, this%gmesh, this%projection_solver_params)
      else ! Everything solid; set velocities to zero and check again in the next timestep.
        this%fluxing_velocity = 0.0_r8
        this%velocity_cc = 0.0_r8
      end if
    end if

    this%t = t+dt
    
  end subroutine step

  ! calculate fluid properties such as various density quantities,
  ! whether cells are immobile, whether faces are solid,
  ! and zero out velocities on immobile cells and solid faces
  subroutine fluid_properties (this, cellRho, fluidRho_n, fluidVof_n, velocity_cc_n, &
      min_fluidRho, solid_face, is_pure_immobile)

    use consts, only: nfc, cutvof, fluid_cutoff

    class(NS_solver), intent(inout) :: this
    real(r8),         intent(out)   :: cellRho(:), fluidRho_n(:), fluidVof_n(:), &
        velocity_cc_n(:,:), min_fluidRho
    logical,          intent(out)   :: solid_face(:), is_pure_immobile(:)

    logical :: real_fluidVof(size(fluidVof_n))
    integer :: i, f
    
    do i = 1,this%mesh%ncell
      ! save previous timestep values
      fluidRho_n(i) = this%fluidRho(i) 
      fluidVof_n(i) = this%fluidVof(i) 
      
      ! update current/next timestep values from the recent interface advection update
      cellRho(i)  = sum(this%vof(:,i)*this%mprop%density) ! TODO: is this used anywhere?
      
      ! TODO: modify fluidRho with fluidDeltaRho, if following the Boussinesq approximation
      this%fluidRho(i) = sum(this%vof(:,i)*this%mprop%density, mask=.not.this%mprop%is_immobile) &
          / merge(this%fluidVof(i), 1.0_r8, this%fluidVof(i) > 0.0_r8)

      this%fluidVof(i) = sum(this%vof(:,i),         mask=.not.this%mprop%is_immobile)
      
      ! vof of mobile and non-void materials
      real_fluidVof(i) = sum(this%vof(:,i), &
          mask=.not.this%mprop%is_immobile .and. .not.this%mprop%is_void)
      
      ! cutRho(i) = sum(cutVof*this%mprop%density, &
      !     mask=.not.this%mprop%is_immobile .and. .not.this%mprop%is_void .and. this%vof(:,i)>0.0_r8)
      
      is_pure_immobile(i) = this%fluidVof(i) < fluid_cutoff !&
          ! ! subroutine turn_off_flow
          ! .or. .not.any(this%gmesh%xc(:,i) < region(:)%x1(:) .or. this%gmesh%xc(:,i) > region(:)%x2(:))

      ! zero out velocities in cells which have completely
      ! solidified between the last flow call and this one
      ! TODO: don't modify velocity here--enforce this condition in projection where it is set
      if (is_pure_immobile(i) .or. this%fluidRho(i) == 0.0_r8) this%velocity_cc(:,i) = 0.0_r8

      ! save previous timestep velocity
      velocity_cc_n(:,i) = this%velocity_cc(:,i)
    end do
    
    min_fluidRho = minval(this%fluidRho*this%fluidVof/real_FluidVof, mask=real_FluidVof > 0.0_r8)
    
    ! a solid face is one where either connected cell is pure immobile
    do f = 1,this%mesh%nface
      solid_face(f) = is_pure_immobile(this%gmesh%fcell(1,f))
      if (this%gmesh%fcell(2,f) > 0) &
          solid_face(f) = solid_face(f) .or. is_pure_immobile(this%gmesh%fcell(2,f))
    end do
    
    ! with some way to grab the local face id of a face given its global id,
    ! this could be merged into the above face loop
    ! for it to work in parallel, would also need to be sure threads are not
    ! trying to write to a cell simultaneously
    ! could also make fluxing_velocity a 1D array over faces, but that is much more work
    ! TODO: don't modify fluxing_velocity here
    do i = 1,this%mesh%ncell
      do f = 1,nfc
        if (solid_face(this%mesh%cface(f,i))) this%fluxing_velocity(f,i) = 0.0_r8
      end do
    end do
    
  end subroutine fluid_properties

  ! decide whether or not fluid flow is needed for this step
  logical function fluid_to_move (this, is_pure_immobile, solid_face)

    use consts, only: nfc

    class(NS_solver), intent(in) :: this
    logical,          intent(in) :: is_pure_immobile(:), solid_face(:)

    integer :: i, f

    ! check if we can skip flow in this timestep
    ! we can skip the flow if there is no inflow of real fluid
    ! and there is no flow within the mesh
    
    fluid_to_move = .false.
    do i = 1,this%mesh%ncell
      if (is_pure_immobile(i)) cycle

      ! check if there is any flow inside the mesh
      do f = 1,nfc
        fluid_to_move = .not.solid_face(f) .and. this%gmesh%cneighbor(f,i) > 0
        
        ! if we have found any fluid motion, we don't need to check the rest of the mesh
        if (fluid_to_move) return
      end do
    end do

    ! TODO: check boundaries
    ! if (boundary_cond(f)==DIRICHLET_PRESSURE) then
    !   if (BC_mat(f) /= NULL_I) fluid_to_move = material_density(BC_mat(f)) > 0.0_r8
    ! else if (boundary_cond(f)==DIRICHLET_VELOCITY) then
    !   ! unfortunately, can't check both these if statements at once
    !   ! because fortran doesn't short-circuit, which can lead to segfaults
    !   if (BC_mat(f) /= NULL_I) &
    !       if (material_density(BC_mat(f)) > 0.0_r8) &
    !       fluid_to_move = dot_product(BC_vel(:,f), gmesh%outnorm(:,f,c)) < 0.0_r8
    ! end if

  end function fluid_to_move

end module NS_solver_type
