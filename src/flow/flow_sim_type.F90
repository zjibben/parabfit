!!
!! FLOW_SIM_TYPE
!!
!! This module defines a class that encapsulates a flow simulation.
!! This drives the time integration and generates the output; it is very basic.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! September 2014
!!

#include "f90_assert.fpp"

module flow_sim_type

  use kinds, only: r8
  use unstr_mesh_type
  ! use mfd_disc_type
  ! use flow_model_type
  use NS_solver_type
  use vof_solver_type
  use logging_services
  use timer_tree_type
  implicit none
  private

  type, public :: flow_sim
     type(unstr_mesh)         , pointer :: mesh       => null()
     ! type(mfd_disc)           , pointer :: disc       => null()
     ! type(HC_model)           , pointer :: model      => null()
     type(NS_solver_t), pointer :: ns_solver  => null()
     type(vof_solver)         , pointer :: vof_solver => null()
     !! Integration control
     real(r8) :: dt_init
     real(r8) :: dt_min
     real(r8), allocatable :: tout(:)
     !!
     integer :: nfile = 0  ! output file counter
   contains
     procedure :: init
     procedure :: run
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
    write(*,*) 'initializing Navier Stokes solver...'
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
      call this%vof_solver%init (this%mesh, this%ns_solver%velocity, plist)
      call this%ns_solver%init_matls(this%vof_solver%cell_matls)
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

    !! Generate the initial material configuration
    write(*,*) 'initializing material layout...'
    call start_timer ('initial-state')
    plist => params%sublist('vof-solver')
    if (plist%is_sublist('initial-vof')) then
      plist => plist%sublist('initial-vof')
      call this%vof_solver%set_initial_state(plist)
   else
      call LS_fatal ('missing "initial-vof" sublist parameter')
    end if
    
    !! Define the initial heat conduction state
    !call this%solver%set_initial_state (t_init, temp, this%dt_init)
    call stop_timer ('initial-state')

    call stop_timer ('initialization')
    
  end subroutine init

  subroutine run (this, stat, errmsg)

    class(flow_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n, status
    real(r8) :: t, hnext
    real(r8), allocatable :: u(:)
    character(80) :: string(2)

    call start_timer ('integration')
    write(*,*) 'running simulation...'
    ! allocate(u(this%model%num_dof()))

    !! Write the initial solution.
    t = 0.0_r8 !this%ns_solver%time()
    call this%write_solution (t)

    ! call LS_info ('')
    ! write(string(1),'(a,es12.5)') 'Beginning integration at T = ', t
    ! call LS_info (string(1))

    ! ! hnext = this%dt_init
    ! ! do n = 1, size(this%tout)
    ! !    call this%ns_solver%integrate  (hnext, status, tout=this%tout(n), hmin=this%dt_min)
    ! !    call this%vof_solver%integrate (hnext, status, tout=this%tout(n), hmin=this%dt_min)
    ! !    t = this%ns_solver%time()
    ! ! end do

    ! ! do n = 1, size(this%tout)
    ! !   if (this%tout(n) <= this%solver%time()) then
    ! !     t = this%tout(n)
        
    ! !     call this%vof_solver%get_interpolated_solution (t, cells)
    
    ! !     call this%write_solution (t, u)
    ! !   else
    ! !     call this%ns_solver%integrate (hnext, status, tout=this%tout(n), hmin=this%dt_min)
    ! !     call this%ns_solver%integrate (hnext, status, tout=this%tout(n), hmin=this%dt_min)
    ! !     t = this%ns_solver%time()
    ! !     select case (status)
    ! !     case (SOLVED_TO_TOUT)
    ! !       call this%solver%get_interpolated_solution (this%tout(n), u)
    ! !       call this%write_solution (this%tout(n), u)
    ! !     case (STEP_FAILED)
    ! !       call this%solver%get_solution_copy (u)
    ! !       call this%write_solution (t, u)
    ! !       stat = -1
    ! !       errmsg = 'failed to take a step'
    ! !       return
    ! !     case (STEP_SIZE_TOO_SMALL)
    ! !       call this%solver%get_solution_copy (u)
    ! !       call this%write_solution (t, u)
    ! !       stat = -1
    ! !       errmsg = 'next time step is too small'
    ! !       return
    ! !     case (BAD_INPUT)
    ! !       stat = -1
    ! !       errmsg = 'bad integrator input parameters'
    ! !       return
    ! !     case default
    ! !       stat = -1
    ! !       errmsg = 'unknown integrator return status'
    ! !       return
    ! !     end select
    ! !     call this%solver%write_metrics (string)
    ! !     call LS_info ('')
    ! !     call LS_info (string(1))
    ! !     call LS_info (string(2))
    ! !   end if
    ! ! end do

    ! call LS_info ('')
    ! write(string(1),'(a,es12.5,a)') 'Completed integration to T = ', this%tout(size(this%tout))
    ! call LS_info (string(1))
    
    stat = 0
    call stop_timer ('integration')

  end subroutine run
  
  !! This auxiliary subroutine writes the solution to a GMV format viz file.
  !! There are lots of better things we could do here (write back to an Exodus
  !! file, VTK, or an HDF5 file), but procedures for writing GMV output were
  !! immediately available.

  subroutine write_solution (this, t)

    use unstr_mesh_gmv
    use vof_tools
    
    class(flow_sim), intent(inout) :: this
    real(r8), intent(in) :: t
    
    integer :: m
    character(21) :: filename,mstr
    !real(r8), dimension(this%vof_solver%nmat,this%mesh%ncell) :: Vof
    real(r8), dimension(this%mesh%ncell,this%vof_solver%nmat) :: Vof

    call start_timer ('output')
    write(*,*) 'dumping solution...'
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
    call matl_get_vof (vof, this%vof_solver%cell_matls, this%vof_solver%nmat, this%vof_solver%matl_id)
    do m = 1,this%vof_solver%nmat ! write all material vofs (TODO should I somehow be writing the interface reconstruction?)
       write(mstr, '(i3)') this%vof_solver%matl_id(m)
       call gmv_write_cell_var (this%mesh, vof(:,m), 'Vof_'//trim(adjustl(mstr)) )
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
