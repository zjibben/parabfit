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
  use NS_solver_type
  use vof_solver_type
  use logging_services
  use timer_tree_type
  implicit none
  private

  type, public :: flow_sim
    type(unstr_mesh), pointer :: mesh       => null()
    type(mesh_geom) :: gmesh
    type(NS_solver), pointer :: ns_solver  => null()
    type(vof_solver), pointer :: vof_solver => null()
    !! Integration control
    real(r8) :: dt_init
    real(r8) :: dt_min
    real(r8), allocatable :: tout(:)
    !!
    integer :: nfile = 0  ! output file counter
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
    ! if (associated(this%disc)) deallocate(this%disc)
    ! if (associated(this%model)) deallocate(this%model)
    if (associated(this%ns_solver)) deallocate(this%ns_solver)
    if (associated(this%vof_solver)) deallocate(this%vof_solver)
  end subroutine flow_sim_delete

  subroutine init (this, params)

    use parameter_list_type
    use unstr_mesh_factory
    use unstr_mesh_func
    use vof_init

    class(flow_sim), intent(out) :: this
    type(parameter_list) :: params

    integer :: stat
    type(parameter_list), pointer :: plist
    character(:), allocatable :: errmsg, context
    !real(r8), allocatable :: temp(:)
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
    call this%gmesh%init (this%mesh)
    call stop_timer ('mesh')
    ! call start_timer ('mfd-discretization')
    ! allocate(this%disc)
    ! call this%disc%init (this%mesh)
    ! call stop_timer ('mfd-discretization')

    ! !! Create the heat conduction model.
    ! call start_timer ('hc-model')
    ! if (params%is_sublist('hc-model')) then
    !   plist => params%sublist('hc-model')
    !   allocate(this%model)
    !   call this%model%init (this%disc, plist)
    ! else
    !   call LS_fatal ('missing "hc-model" sublist parameter')
    ! end if
    ! call stop_timer ('hc-model')

    !! Create the navier stokes solver.
    write(*,*) 'initializing Navier-Stokes solver...'
    call start_timer ('ns-solver')
    if (params%is_sublist('ns-solver')) then
      plist => params%sublist('ns-solver')
      allocate(this%ns_solver)
      call this%ns_solver%init (this%mesh, plist)
    else
      call LS_fatal ('missing "ns-solver" sublist parameter')
    end if
    call stop_timer ('ns-solver')

    !! Create the volume of fluid solver.
    write(*,*) 'initializing vof solver...'
    call start_timer ('vof-solver')
    if (params%is_sublist('vof-solver')) then
      plist => params%sublist('vof-solver')
      allocate(this%vof_solver)
      call this%vof_solver%init (this%mesh, this%gmesh, this%ns_solver%velocity, this%ns_solver%fluidRho, plist)
      call this%ns_solver%init_matls(this%vof_solver%vof)
    else
      call LS_fatal ('missing "vof-solver" sublist parameter')
    end if
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
      if (this%dt_min > this%dt_init) call LS_fatal (context//'require "min-time-step" <= "initial-time-step"')
      call plist%get ('output-times', this%tout, stat=stat, errmsg=errmsg)
      if (stat /= 0) call LS_fatal (context//errmsg)
      !TODO: check for strictly increasing values in TOUT, TOUT > t_init, or sort
      !and cull those < t_init.
    else
      call LS_fatal ('missing "sim-control" sublist parameter')
    end if

    !! Set the initial flow field
    write(*,*) 'initializing velocity field...'
    call start_timer ('initial-state')
    call this%ns_solver%set_initial_state()
    ! plist => params%sublist('vof-solver')
    ! if (plist%is_sublist('initial-vof')) then
    !    plist => plist%sublist('initial-vof')
    !    call this%vof_solver%set_initial_state(plist)
    ! else
    !    call LS_fatal ('missing "initial-vof" sublist parameter')
    ! end if

    !! Generate the initial material configuration
    write(*,*) 'initializing material layout...'
    plist => params%sublist('vof-solver')
    if (plist%is_sublist('initial-vof')) then
      plist => plist%sublist('initial-vof')
      call this%vof_solver%set_initial_state(plist)
    else
      call LS_fatal ('missing "initial-vof" sublist parameter')
    end if

    call stop_timer ('initial-state')

    call stop_timer ('initialization')

  end subroutine init

  ! this is the main driver, only separating dumps
  ! it calls a separate routine for flow subcycles
  subroutine run (this, stat, errmsg)
    !use velocity_to_faces_func

    class(flow_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n, status
    real(r8) :: t, hnext
    real(r8), allocatable :: u(:)
    character(80) :: string(2)

    call start_timer ('integration')
    write(*,*) 'running simulation...'

    !! Write the initial solution.
    t = 0.0_r8
    call this%write_solution (t)

    call LS_info ('')
    write(string(1),'(a,es12.5)') 'Beginning integration at T = ', t
    call LS_info (string(1))

    do n = 1, size(this%tout)
      ! !call this%ns_solver%

      ! ! call velocity_to_faces(this%vof_solver%velocity, this%ns_solver%velocity, this%mesh, t, &
      ! !      this%ns_solver%use_prescribed_velocity, this%ns_solver%prescribed_velocity_case)
      ! call velocity_to_faces(this%vof_solver%fluxing_velocity, this%ns_solver%velocity, this%mesh, t, &
      !      this%ns_solver%use_prescribed_velocity, this%ns_solver%prescribed_velocity_case)

      ! call this%vof_solver%advect_mass(this%tout(n)-t)

      call this%step (this%tout(n)-t, t)

      t = this%tout(n)
      call this%write_solution (t)

      !call this%solver%write_metrics (string)
      !call LS_info ('')
      write(string(1),'(a,es12.5,a)') 'Completed integration to T = ', this%tout(n)
      call LS_info (string(1))
      !call LS_info (string(2))
    end do

    call LS_info ('')
    write(string(1),'(a,es12.5,a)') 'Completed integration to T = ', this%tout(size(this%tout))
    call LS_info (string(1))

    stat = 0
    call stop_timer ('integration')

  end subroutine run

  subroutine step (this,dt,t)
    use velocity_to_faces_func

    class(flow_sim), intent(inout) :: this
    real(r8), intent(in) :: dt,t

    real(r8) :: flow_dt,tlocal,CFL
    integer :: ns_subcycles
    
    tlocal = 0.0_r8
    
    ! ! set the flow timestep. this can be improved, probably needs to be corrected for non-cubic cell
    ! call velocity_to_faces (this%vof_solver%fluxing_velocity, this%ns_solver%velocity, this%mesh, this%gmesh, t, &
    !      this%ns_solver%use_prescribed_velocity, this%ns_solver%prescribed_velocity_case)
    ! ns_subcycles = ceiling(dt / (0.25_r8*minval(this%mesh%volume**(1.0_r8/3.0_r8))/maxval(this%vof_solver%fluxing_velocity))) 
    ! flow_dt = dt / real(ns_subcycles,r8)
    
    do while (tlocal<dt)
      !call this%ns_solver%update_flowfield (flow_dt)

      ! project the cell centered velocity from the flowsolver to the faces
      call velocity_to_faces (this%vof_solver%fluxing_velocity, this%ns_solver%velocity, this%mesh, this%gmesh, t+tlocal, &
           this%ns_solver%use_prescribed_velocity, this%ns_solver%prescribed_velocity_case)

      ! update the timestep
      CFL = 0.25_r8
      ns_subcycles = max(&
           ceiling((dt-tlocal) / (CFL*minval(this%mesh%volume**(1.0_r8/3.0_r8))/maxval(this%vof_solver%fluxing_velocity))),&
           1)
      flow_dt = (dt-tlocal) / real(ns_subcycles,r8)

      ! advect material
      call this%vof_solver%advect_mass (flow_dt)

      ! increment the local time
      tlocal = tlocal + flow_dt
    end do
    
  end subroutine step

  !! This auxiliary subroutine writes the solution to a GMV format viz file.
  !! There are lots of better things we could do here (write back to an Exodus
  !! file, VTK, or an HDF5 file), but procedures for writing GMV output were
  !! immediately available.

  subroutine write_solution (this, t)
    use unstr_mesh_gmv
    
    class(flow_sim), intent(inout) :: this
    real(r8), intent(in) :: t

    integer :: m
    character(21) :: filename,mstr

    call start_timer ('output')

    write(filename,'(a,i4.4)') 'flow_out.gmv.', this%nfile
    this%nfile = this%nfile + 1

    call gmv_open (trim(filename))
    call gmv_write_unstr_mesh (this%mesh) ! segfault in here
    call gmv_begin_variables (time=t)

    ! write velocities (TODO: how do I write a vector to each cell?)
    call gmv_write_cell_var (this%mesh, this%ns_solver%velocity(1,:), 'u1')
    call gmv_write_cell_var (this%mesh, this%ns_solver%velocity(2,:), 'u2')
    call gmv_write_cell_var (this%mesh, this%ns_solver%velocity(3,:), 'u3')

    ! get the vof
    !call matl_get_vof (vof, this%vof_solver%cell_matls, this%vof_solver%nmat, this%vof_solver%matl_id)
    do m = 1,this%vof_solver%nmat ! write all material vofs (TODO should I somehow be writing the interface reconstruction?)
      write(mstr, '(i3)') this%vof_solver%matl_id(m)
      call gmv_write_cell_var (this%mesh, this%vof_solver%vof(m,:), 'Vof_'//trim(adjustl(mstr)) )
    end do

    call gmv_end_variables
    call gmv_write_unstr_mesh_surf (this%mesh)
    ! call gmv_begin_surfvars
    ! call gmv_write_surf_var (this%mesh, tface, 'T')
    ! call gmv_end_surfvars
    call gmv_close

    call stop_timer ('output')

  end subroutine write_solution

end module flow_sim_type
