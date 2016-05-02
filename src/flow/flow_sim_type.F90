!!
!! FLOW_SIM_TYPE
!!
!! This module defines a class that encapsulates a flow simulation.
!! This drives the time integration and generates the output; it is very basic.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

#include "f90_assert.fpp"

module flow_sim_type

  use kinds, only: r8
  use unstr_mesh_type
  use mesh_geom_type
  use matl_props_type
  use NS_solver_type
  use vof_solver_type
  use logging_services
  use timer_tree_type
  implicit none
  private

  type, public :: flow_sim
    type(unstr_mesh), pointer :: mesh       => null()
    type(mesh_geom),  pointer :: gmesh      => null()
    type(NS_solver),  pointer :: ns_solver  => null()
    type(vof_solver), pointer :: vof_solver => null()
    type(matl_props), pointer :: mprop      => null()

    !! Integration control
    real(r8) :: dt_init
    real(r8) :: dt_min
    real(r8), allocatable :: tout(:)
    !!
    integer :: nfile = 0  ! output file counter
    logical :: dump_intrec
  contains
    procedure :: init
    procedure :: run
    procedure, private :: step
    procedure, private :: write_solution
    final :: flow_sim_delete
  end type flow_sim

contains

  !! Final subroutine for HC_SIM objects.
  subroutine flow_sim_delete (this)
    type(flow_sim), intent(inout) :: this
    if (associated(this%mesh)) deallocate(this%mesh)
    if (associated(this%gmesh)) deallocate(this%gmesh)
    if (associated(this%ns_solver)) deallocate(this%ns_solver)
    if (associated(this%vof_solver)) deallocate(this%vof_solver)
  end subroutine flow_sim_delete

  subroutine init (this, params)

    use parameter_list_type
    use unstr_mesh_factory
    use unstr_mesh_func
    use vof_init
#ifdef _OPENMP
    use omp_lib, only: omp_get_max_threads
#endif

    class(flow_sim), intent(out) :: this
    type(parameter_list) :: params

    integer :: stat
    type(parameter_list), pointer :: plist
    character(:), allocatable :: errmsg, context
    integer, allocatable :: tmp(:)
    real(r8) :: t_init

    call start_timer ('initialization')
    call LS_info ('Initializing the simulation', LS_VERB_NOISY)

    !! Create the mesh and discretization object.
    write(*,*) 'initializing mesh...'
    call start_timer ('mesh')
    if (params%is_sublist('mesh')) then
      plist => params%sublist('mesh')
      this%mesh => new_unstr_mesh(plist)
    else
      call LS_fatal ('missing "mesh" sublist parameter')
    end if
    allocate(this%gmesh)
    call this%gmesh%init (this%mesh)
    call stop_timer ('mesh')

    !! Create the navier stokes solver.
    write(*,*) 'initializing Navier-Stokes solver...'
    call start_timer ('ns-solver')
    if (params%is_sublist('ns-solver')) then
      plist => params%sublist('ns-solver')
      allocate(this%ns_solver)
      call this%ns_solver%init (this%mesh, this%gmesh, plist)
    else
      call LS_fatal ('missing "ns-solver" sublist parameter')
    end if
    call stop_timer ('ns-solver')

    !! Load material properties
    write(*,*) 'Reading material properties...'
    if (params%is_sublist('material-properties')) then
      plist => params%sublist('material-properties')
      allocate(this%mprop)
      call this%mprop%init (plist)
      call this%ns_solver%init_mprop (this%mprop)
    else if (this%ns_solver%use_prescribed_velocity) then
      allocate(this%mprop)
      ! get the number of materials from the material id array
      if (params%is_sublist('vof-solver')) then
        plist => params%sublist('vof-solver')
        call plist%get ('materials', tmp, stat=stat, errmsg=errmsg)
        if (stat /= 0) call LS_fatal (context//errmsg)
        ! set all materials to non-void
        allocate(this%mprop%is_void(size(tmp)))
        this%mprop%is_void = .false.
      else
        call LS_fatal ('missing "vof-solver" sublist parameter')
      end if
    else
      call LS_fatal ('missing "material-properties" sublist parameter')
    end if
    
    !! Create the volume of fluid solver.
    write(*,*) 'initializing vof solver...'
    call start_timer ('vof-solver')
    if (params%is_sublist('vof-solver')) then
      plist => params%sublist('vof-solver')
      allocate(this%vof_solver)
      call this%vof_solver%init (this%mesh, this%gmesh, this%mprop%nmat, &
          this%ns_solver%fluxing_velocity, this%ns_solver%fluidRho, plist, &
          this%ns_solver%velocity_bc)
      call this%ns_solver%init_matls(this%vof_solver%vof, this%vof_solver%volume_flux_tot)
    else
      call LS_fatal ('missing "vof-solver" sublist parameter')
    end if

    ! print plane reconstructions 
    this%dump_intrec = .true. ! TODO: grab from input file
#ifdef _OPENMP
    ! only dump interface reconstruction if running in serial
    this%dump_intrec = this%dump_intrec .and. omp_get_max_threads() == 1
#endif
    call stop_timer ('vof-solver')

    !! Simulation control parameters
    write(*,*) 'loading simulation control parameters...'
    if (params%is_sublist('sim-control')) then
      plist => params%sublist('sim-control')
      context = 'processing ' // plist%name() // ': '
      call plist%get ('initial-time', t_init, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      call plist%get ('initial-time-step', this%dt_init, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (this%dt_init <= 0.0_r8) call LS_fatal (context//'"initial-time-step" must be > 0.0')
      call plist%get ('min-time-step', this%dt_min, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      if (this%dt_min > this%dt_init) &
          call LS_fatal (context//'require "min-time-step" <= "initial-time-step"')
      call plist%get ('output-times', this%tout, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      !TODO: check for strictly increasing values in TOUT, TOUT > t_init, or sort
      !and cull those < t_init.
    else
      call LS_fatal ('missing "sim-control" sublist parameter')
    end if

    !! set initial state
    call start_timer ('initial-state')

    !! Generate the initial material configuration
    write(*,*) 'initializing material layout...'
    plist => params%sublist('vof-solver')
    if (plist%is_sublist('initial-vof')) then
      plist => plist%sublist('initial-vof')
      call this%vof_solver%set_initial_state(plist)
    else
      call LS_fatal ('missing "initial-vof" sublist parameter')
    end if

    !! Set the initial flow, pressure, density, fluidvof fields
    write(*,*) 'initializing flow...'
    call this%ns_solver%set_initial_state ()

    call stop_timer ('initial-state')

    call stop_timer ('initialization')

  end subroutine init

  ! this is the main driver, only separating dumps
  ! it calls a separate routine for flow subcycles
  subroutine run (this, stat, errmsg)
    
    class(flow_sim),           intent(inout) :: this
    integer,                   intent(out)   :: stat
    character(:), allocatable, intent(out)   :: errmsg

    integer       :: n
    real(r8)      :: t
    character(80) :: string(2)
    
    call start_timer ('integration')
    write(*,*) 'running simulation...'

    !! Write the initial solution.
    t = 0.0_r8
    if (this%dump_intrec) call this%vof_solver%update_intrec_surf ()
    call this%write_solution (t)

    call LS_info ('')
    write(string(1),'(a,es12.5)') 'Beginning integration at T = ', t
    call LS_info (string(1))

    do n = 1, size(this%tout)
      call this%step (this%tout(n)-t, t)

      t = this%tout(n)
      call this%write_solution (t)

      write(string(1),'(a,es12.5,a)') 'Completed integration to T = ', this%tout(n)
      call LS_info (string(1))
      !call LS_info (string(2))
    end do

    call LS_info ('')
    write(string(1),'(a,es12.5,a)') 'Completed integration to T = ', this%tout(size(this%tout))
    call LS_info (string(1))
    
    stat = 0
    call stop_timer ('integration')

    ! print L1, L2, and Linf error norms (TODO: optionally)
    call this%vof_solver%print_error_norms ()
  end subroutine run

  ! update the flow from t to t+dt, subcycling as necessary
  subroutine step (this, dt, t)

    class(flow_sim), intent(inout) :: this
    real(r8),        intent(in)    :: dt, t

    real(r8)           :: flow_dt, tlocal
    integer, save      :: iter = 0

    ! ! DEBUGGING/SCALING ##################
    ! integer, parameter :: nsteps = 1
    ! write(*,*) 'WARNING - number of timesteps forced'
    ! ! DEBUGGING/SCALING ##################

    tlocal = 0.0_r8
    do while (tlocal<dt) ! .and. iter<nsteps)
      ! update the timestep size
      flow_dt = this%ns_solver%timestep_size (dt-tlocal)
      
      ! advect material
      call this%vof_solver%advect_mass (flow_dt, this%mprop%is_void, &
          this%dump_intrec .and. flow_dt==dt-tlocal)

      ! update the pressure and velocity fields
      call this%ns_solver%step (flow_dt, t+tlocal)

      ! increment the local time
      tlocal = tlocal + flow_dt
      iter = iter+1
    end do

    write(*,'(a,i6,es14.4,a,2es20.10)') '  cumulative iterations, dt: ', iter, flow_dt, &
        ',   minmaxvel: ', &
        minval(this%ns_solver%velocity_cc(2,:)), maxval(this%ns_solver%velocity_cc(2,:))
    call this%vof_solver%print_vofs ()
    
  end subroutine step


  !! This auxiliary subroutine writes the solution to a GMV format viz file.
  !! There are lots of better things we could do here (write back to an Exodus
  !! file, VTK, or an HDF5 file), but procedures for writing GMV output were
  !! immediately available.
  subroutine write_solution (this, t)

    use unstr_mesh_gmv
    use array_utils, only: int2str
    
    class(flow_sim), intent(inout) :: this
    real(r8), intent(in) :: t

    integer       :: m
    character(21) :: nstr

    call start_timer ('output')

    write(nstr,'(i4.4)') this%nfile
    
    call gmv_open (trim('flow_out.gmv.'//trim(nstr)))
    call gmv_write_unstr_mesh (this%mesh)
    call gmv_begin_variables (time=t)

    ! write velocities (TODO: how do I write a vector to each cell?)
    call gmv_write_cell_var (this%mesh, this%ns_solver%velocity_cc(1,:), 'u1')
    call gmv_write_cell_var (this%mesh, this%ns_solver%velocity_cc(2,:), 'u2')
    call gmv_write_cell_var (this%mesh, this%ns_solver%velocity_cc(3,:), 'u3')
    
    call gmv_write_cell_var (this%mesh, this%ns_solver%pressure_cc, 'pressure')

    ! write all material vofs
    do m = 1,this%vof_solver%nmat
      call gmv_write_cell_var (this%mesh, this%vof_solver%vof(m,:), 'vof_'//trim(int2str(m)) )
    end do

    call gmv_end_variables
    call gmv_write_unstr_mesh_surf (this%mesh)
    ! call gmv_begin_surfvars
    ! call gmv_write_surf_var (this%mesh, tface, 'T')
    ! call gmv_end_surfvars
    call gmv_close

    ! dump the interface reconstruction
    if (this%dump_intrec) then
      do m = 1,size(this%vof_solver%intrec)
        call this%vof_solver%intrec(m)%write_ply ('surf'//trim(int2str(m))//'_'//trim(nstr)//'.ply')
      end do
    end if
    this%nfile = this%nfile + 1
    
    call stop_timer ('output')

  end subroutine write_solution

end module flow_sim_type
