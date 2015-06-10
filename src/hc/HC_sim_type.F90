!!
!! HC_SIM_TYPE
!!
!! This module defines a class that encapsulates a heat conduction simulation.
!! This drives the time integration and generates the output; it is very basic.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! September 2014
!!

#include "f90_assert.fpp"

module HC_sim_type

  use kinds, only: r8
  use unstr_mesh_type
  use unstr_mesh_partition_type
  use mfd_disc_type
  use HC_model_type
  use HC_solver_type
  use logging_services
  use timer_tree_type
  implicit none
  private

  type, public :: HC_sim
    type(unstr_mesh), pointer :: mesh => null()
    type(unstr_mesh_partition), pointer :: partition => null()
    type(mfd_disc), pointer :: disc => null()
    type(HC_model), pointer :: model => null()
    type(HC_solver), pointer :: solver => null()
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
    final :: HC_sim_delete
  end type HC_sim

contains

  !! Final subroutine for HC_SIM objects.
  subroutine HC_sim_delete (this)
    type(HC_sim), intent(inout) :: this
    if (associated(this%mesh)) deallocate(this%mesh)
    if (associated(this%partition)) deallocate(this%partition)
    if (associated(this%disc)) deallocate(this%disc)
    if (associated(this%model)) deallocate(this%model)
    if (associated(this%solver)) deallocate(this%solver)
  end subroutine HC_sim_delete

  subroutine init (this, params)

    use parameter_list_type
    use unstr_mesh_factory
    use unstr_mesh_func

    class(HC_sim), intent(out) :: this
    type(parameter_list) :: params

    integer :: stat
    type(parameter_list), pointer :: plist
    character(:), allocatable :: errmsg, context
    real(r8), allocatable :: temp(:)
    real(r8) :: t_init

    ! Enable wall-clock timing.
    call set_timer_type(realtime_timing)

    call start_timer ('initialization')
    call LS_info ('Initializing the simulation', LS_VERB_NOISY)

    !! Create the mesh
    call start_timer ('mesh')
    if (params%is_sublist('mesh')) then
      plist => params%sublist('mesh')
      this%mesh => new_unstr_mesh(plist)
    else
      call LS_fatal ('missing "mesh" sublist parameter')
    end if
    call stop_timer ('mesh')

    !! Create the partition type
    call start_timer ('partition')
    if (params%is_sublist('partition-layout')) then
      plist => params%sublist('partition-layout')
      allocate(this%partition)
      call this%partition%init (this%mesh, plist)
    else
      call LS_fatal ('missing "partition-layout" sublist parameter')
    end if
    call stop_timer ('partition')

    !! Create the discretization object.
    call start_timer ('mfd-discretization')
    allocate(this%disc)
    call this%disc%init (this%mesh, this%partition)
    call stop_timer ('mfd-discretization')

    !! Create the heat conduction model.
    call start_timer ('hc-model')
    if (params%is_sublist('hc-model')) then
      plist => params%sublist('hc-model')
      allocate(this%model)
      call this%model%init (this%disc, plist)
    else
      call LS_fatal ('missing "hc-model" sublist parameter')
    end if
    call stop_timer ('hc-model')

    !! Create the heat conduction solver.
    call start_timer ('hc-solver')
    if (params%is_sublist('hc-solver')) then
      plist => params%sublist('hc-solver')
      allocate(this%solver)
      call this%solver%init (this%model, plist)
    else
      call LS_fatal ('missing "hc-solver" sublist parameter')
    end if
    call stop_timer ('hc-solver')

    !! Simulation control parameters
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

    !! Generate the initial temperature field
    call start_timer ('initial-state')
    if (params%is_sublist('initial-temperature')) then
      plist => params%sublist('initial-temperature')
      allocate(temp(this%mesh%ncell))
      call compute_mesh_func (this%mesh, plist, temp)
    else
      call LS_fatal ('missing "initial-temperature" sublist parameter')
    end if

    !! Define the initial heat conduction state
    call this%solver%set_initial_state (t_init, temp, this%dt_init)
    call stop_timer ('initial-state')

    call stop_timer ('initialization')

  end subroutine init

  subroutine run (this, stat, errmsg)

    class(HC_sim), intent(inout) :: this
    integer, intent(out) :: stat
    character(:), allocatable, intent(out) :: errmsg

    integer :: n, status
    real(r8) :: t, hnext
    real(r8), allocatable :: u(:)
    character(80) :: string(2)

    call start_timer ('integration')

    allocate(u(this%model%num_dof()))

    !! Write the initial solution.
    t = this%solver%time()
    call this%solver%get_solution_copy (u)
    call this%write_solution (t, u)

    call LS_info ('')
    write(string(1),'(a,es12.5)') 'Beginning integration at T = ', t
    call LS_info (string(1))

    hnext = this%dt_init
    do n = 1, size(this%tout)
      if (this%tout(n) <= this%solver%time()) then
        t = this%tout(n)
        call this%solver%get_interpolated_solution (t, u)
        call this%write_solution (t, u)
      else
        call this%solver%integrate (hnext, status, tout=this%tout(n), hmin=this%dt_min)
        t = this%solver%time()
        select case (status)
        case (SOLVED_TO_TOUT)
          call this%solver%get_interpolated_solution (this%tout(n), u)
          call this%write_solution (this%tout(n), u)
        case (STEP_FAILED)
          call this%solver%get_solution_copy (u)
          call this%write_solution (t, u)
          stat = -1
          errmsg = 'failed to take a step'
          return
        case (STEP_SIZE_TOO_SMALL)
          call this%solver%get_solution_copy (u)
          call this%write_solution (t, u)
          stat = -1
          errmsg = 'next time step is too small'
          return
        case (BAD_INPUT)
          stat = -1
          errmsg = 'bad integrator input parameters'
          return
        case default
          stat = -1
          errmsg = 'unknown integrator return status'
          return
        end select
        call this%solver%write_metrics (string)
        call LS_info ('')
        call LS_info (string(1))
        call LS_info (string(2))
      end if
    end do

    call LS_info ('')
    write(string(1),'(a,es12.5,a)') 'Completed integration to T = ', this%tout(size(this%tout))
    call LS_info (string(1))

    stat = 0
    call stop_timer ('integration')

  end subroutine run

  !! This auxiliary subroutine writes the solution to a GMV format viz file.
  !! There are lots of better things we could do here (write back to an Exodus
  !! file, VTK, or an HDF5 file), but procedures for writing GMV output were
  !! immediately available.

  subroutine write_solution (this, t, u)

    use unstr_mesh_gmv

    class(HC_sim), intent(inout) :: this
    real(r8), intent(in) :: t
    real(r8), intent(in), target :: u(:)

    character(16) :: filename
    real(r8), pointer :: hcell(:), tcell(:), tface(:)

    call start_timer ('output')

    write(filename,'(a,i4.4)') 'out.gmv.', this%nfile
    this%nfile = this%nfile + 1

    call this%model%get_cell_heat_view (u, hcell)
    call this%model%get_cell_temp_view (u, tcell)
    call this%model%get_face_temp_view (u, tface)

    call gmv_open (trim(filename))
    call gmv_write_unstr_mesh (this%mesh)
    call gmv_begin_variables (time=t)
    call gmv_write_cell_var (this%mesh, hcell, 'H')
    call gmv_write_cell_var (this%mesh, tcell, 'T')
    call gmv_end_variables
    call gmv_write_unstr_mesh_surf (this%mesh)
    call gmv_begin_surfvars
    call gmv_write_surf_var (this%mesh, tface, 'T')
    call gmv_end_surfvars
    call gmv_close

    call stop_timer ('output')

  end subroutine write_solution

end module HC_sim_type
