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
  use mesh_geom_type
  use logging_services
  implicit none
  private


  type, public :: vof_solver
    private
    type(unstr_mesh), pointer :: mesh        => null() ! reference only -- do not own
    type(mesh_geom),  pointer :: gmesh       => null() ! reference only -- do not own
    real(r8), pointer         :: fluidRho(:) => null() ! reference only -- do not own
    
    !! Pending/current state
    real(r8)                      :: t, dt
    integer, public               :: nmat          ! number of materials present globally
    !real(r8), public, allocatable :: velocity(:,:) ! face velocities
    real(r8), public, allocatable :: fluxing_velocity(:,:) ! face normal velocities
    integer,  allocatable, public :: matl_id(:)
    real(r8), allocatable, public :: vof(:,:)
  contains
    procedure :: init
    procedure :: set_initial_state
    ! procedure :: test_initial_state
    procedure :: advect_mass
    procedure, private :: advect_volume
    !procedure, private :: flux_renorm
    ! procedure :: time
    ! procedure :: get_interpolated_solution
    ! procedure :: get_solution_view
    ! procedure :: get_solution_copy
    ! procedure :: write_metrics
    ! procedure :: advance_state
    ! procedure :: commit_pending_state
  end type vof_solver

  real(r8), parameter :: cutvof = 1.0e-8_r8
contains
  
  subroutine init (this, mesh, gmesh, velocity_cc, fluidRho, params)

    use parameter_list_type
    use vof_init

    class(vof_solver), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(mesh_geom),  intent(in), target :: gmesh
    real(r8), intent(in), target :: velocity_cc(:,:), fluidRho(:)
    type(parameter_list) :: params

    integer :: stat
    type(parameter_list), pointer :: plist
    character(:), allocatable :: context,errmsg

    !this%velocity_cc => velocity_cc
    this%mesh     => mesh
    this%gmesh    => gmesh
    this%fluidRho => fluidRho
    allocate(this%fluxing_velocity(6,this%mesh%ncell))
    ! allocate(this%velocity(3,this%mesh%nface))

    context = 'processing ' // params%name() // ': '

    call params%get ('materials', this%matl_id, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)

    this%nmat = size(this%matl_id)

    allocate(this%vof(this%nmat,this%mesh%ncell))

  end subroutine init

  subroutine set_initial_state (this, plist)
    use parameter_list_type
    use vof_init
    !use vof_tools, only: distinct_matls, volume_of_matl

    class(vof_solver), intent(inout) :: this
    type(parameter_list), intent(in) :: plist

    !! Initialize the Vof
    call vof_initialize (this%mesh, plist, this%vof, this%matl_id, this%nmat)

    write(*,*) sum(this%vof(1,:)*this%mesh%volume(:))
    write(*,*) sum(this%vof(2,:)*this%mesh%volume(:))

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



  subroutine advect_mass(this,dt)
    use timer_tree_type
    class(vof_solver), intent(inout) :: this
    real(r8), intent(in) :: dt
    
    integer  :: stat
    real(r8) :: Vof_n(this%nmat, this%mesh%ncell)

    ! if (.not.fluid_to_move) return

    ! Start Volume Advection Timer
    call start_timer("Mass Advection")
    
    ! if (volume_track_interfaces .and. .not.associated(VT_Interface_Mask)) then
    !    allocate (VT_Interface_Mask(nfc,ncells), stat = stat)
    !    if (stat /= 0) call LS_fatal ('ADVECT_MASS: VT_Interface_Mask(nfc,ncells) allocation failed')
    ! end if
    
    ! need a copy of the time n Vof values.
    Vof_n = this%vof

    ! Advect material volumes.
    call this%advect_volume (Vof_n,dt)

    ! ! Update the mass & concentration distributions.
    ! if (volume_track_interfaces) then 
    !    call update_mass (this%fluxing_velocity, this%vof, Vof_n, VT_Interface_Mask)
    ! else
    !    call update_mass (this%fluxing_velocity, this%vof, Vof_n)
    ! end if

    ! ! Accumulate Momentum Advection Array
    ! call advect_momentum_accumulation (this%volume_flux)

    ! ! Find and store the accumulated inflow and outflow mass.
    ! call inflow_outflow_update ()

    ! Stop the Volume Advection Timer
    call stop_timer("Mass Advection")
    
  end subroutine advect_mass

  ! subroutine inflow_outflow_update (qin, qout)
  !   real(r8), intent(inout) :: qin,qout
    
  !   real(r8), dimension(this%mesh%ncell) :: Tmp,mask
  !   integer :: f,stat
  !   integer, parameter :: nfc = 6

  !   ! loop through all faces
  !   do f = 1, nfc
  !      mask = in_flow (f, this%fluxing_velocity) .or. out_flow(f, this%fluxing_velocity)

  !      if (any(mask)) then
  !         ! accumulate inflow volume - this%fluxing_velocity < zero
  !         tmp = min(dt*cell%face_area(f)*this%fluxing_velocity(f,:),0.0_r8)
  !         qin = qin + abs(pgslib_global_sum(tmp, mask = mask))

  !         ! accumulate outflow volume - this%fluxing_velocity > zero
  !         tmp = max(dt*cell%face_area(f)*this%fluxing_velocity(f,:),0.0_r8)
  !         qout = qout + pgslib_global_sum(tmp, mask = mask)
  !      end if
  !   end do
  ! end subroutine inflow_outflow_update
  
  subroutine advect_volume (this, vof_n, dt)
    use timer_tree_type
    use volume_track_module, only: volume_track

    class(vof_solver), intent(inout) :: this
    real(r8),          intent(in)    :: Vof_n(this%nmat,this%mesh%ncell), dt

    integer  :: status, p, vps, volume_track_subcycles
    real(r8) :: adv_dt
    real(r8) :: Volume_Flux_Sub(this%nmat, 6, this%mesh%ncell) ! keeps track of volume changes in every subcycle
    real(r8) :: volume_flux_tot(this%nmat, 6, this%mesh%ncell)
    
    ! Start the volume advection timer.
    call start_timer("Volume Tracking")
    
    ! Zero the total flux array.
    volume_flux_tot = 0.0_r8

    ! Set the advection timestep. this can be improved.
    volume_track_subcycles = 1
    adv_dt = dt/real(volume_track_subcycles,r8)
    
    do p = 1,volume_track_subcycles
      ! Get the donor fluxes.
      volume_flux_sub = volume_track (adv_dt, this%mesh, this%gmesh, this%vof, this%fluxing_velocity, this%nmat, this%fluidRho)
      
      ! Normalize the donor fluxes.
      !write(*,*) 'warning: flux renorm deactivated'
      call flux_renorm (this%fluxing_velocity, Vof_n, Volume_Flux_Tot, Volume_Flux_Sub, adv_dt, this%mesh)
      
      !write(*,*) 'current volumes: ',this%vof(:,501) * this%mesh%volume(501)
      ! Compute the acceptor fluxes.
      call flux_acceptor (volume_flux_sub, this%gmesh)

      ! Compute BC (inflow) fluxes.
      call flux_bc (this%fluxing_velocity, vof_n, volume_flux_sub, adv_dt, this%mesh, this%gmesh)
      
      ! Add the volume fluxes from this subcycle (Volume_Flux_Sub) to the
      ! total flux array (Volume_Flux_Tot), and update the volume fraction
      ! array (Vof).
      call volume_advance (Volume_Flux_Sub, volume_flux_tot, this%vof, this%mesh%volume)
      
      ! ! Make sure volume fractions of a particular material are within
      ! ! the allowed range (0 <= Vof <= 1) and that all materials sum to one.
      ! call this%vof_bounds ()

      !write(*,*) 'completed subcycle, dt =',adv_dt
    end do

    ! Stop the volume advection timer.
    call stop_timer("Volume Tracking")

  end subroutine advect_volume

  !=======================================================================
  ! Purpose(s):
  !
  !   Scan all faces with an outward flux and determine if any
  !   material is over-exhausted from this cell.  If so lower
  !   the fluxes until the material is just exhausted.  Then
  !   loop over the faces and balance the individual material
  !   fluxes with the total face flux. The sum of the material
  !   volume fluxes (Volume_Flux_Sub) for each face should sum
  !   to the total volume flux for that face.
  !
  !=======================================================================
  subroutine flux_renorm (Fluxing_Velocity, Vof_n, Volume_Flux_Tot, Volume_Flux_Sub, adv_dt, mesh)
    use unstr_mesh_type

    real(r8), intent(in)    :: fluxing_velocity(:,:), vof_n(:,:), adv_dt, volume_flux_tot(:,:,:)
    real(r8), intent(inout) :: volume_flux_sub(:,:,:)
    class(unstr_mesh), intent(in) :: mesh

    integer :: n,ierr

    do n = 1,mesh%ncell
      call renorm_cell (volume_flux_sub(:,:,n), volume_flux_tot(:,:,n), fluxing_velocity(:,n), vof_n(:,n), adv_dt, &
           adv_dt*Fluxing_Velocity(:,n)*mesh%area(mesh%cnode(:,n)), mesh%volume(n), ierr)
      if (ierr /= 0) then
        write(*,*) 'cell id:       ',n
        write(*,*) 'cell centroid: ',sum(mesh%x(:,mesh%cnode(:,n)), dim=2) / 8.0_r8
        call LS_fatal ('FLUX_RENORM: cannot reassign face flux to any other material')
      end if
    end do ! cell loop
    
  end subroutine flux_renorm

  subroutine renorm_cell (volume_flux_sub, volume_flux_tot, fluxing_velocity, vof_n, adv_dt, total_face_flux, cell_volume, ierr)
    !use fluid_data_module,    only: isImmobile

    real(r8), intent(in)    :: fluxing_velocity(:), vof_n(:), adv_dt, volume_flux_tot(:,:), total_face_flux(:), cell_volume
    real(r8), intent(inout) :: volume_flux_sub(:,:)
    integer,  intent(out)   :: ierr

    real(r8) :: Ratio, total_material_flux, cumulative_outward_material_flux, total_flux_through_face,total_flux_through_face_not_maxed
    integer  :: norm_iter, f, m, number_not_maxed
    logical  :: flux_reduced, maxed(size(vof_n))

    maxed = .false.
    ierr = 0
    
    ! Loop over the renorm_loop a maximum of nmat - 1 times, to resolve all instances
    ! of fluxing more material than was in the cell at the beginning of the timestep
    do norm_iter = 1, size(vof_n)+1
      flux_reduced = .false.

      ! We stay in this loop until the sum of material fluxes from each face of
      ! the cell equals the total face volume flux.  Where the cumulative sum of 
      ! individual material fluxes (from this and previous volume_track_subcycles) 
      ! exceeds the volume of a particular material originally within a 
      ! cell, we decrease those fluxes to equal the volume of material still available
      ! to be fluxed, and increase other fluxes appropriately.  If this increase
      ! leads to the over-exhaustion of other materials, we work our way through
      ! this loop again, and again, and again, and ... , until we're done.

      ! The first step is to determine if any material is being over-exhausted from 
      ! a cell.  If so mark it as MAXED and lower the Volume_Flux_Sub's so that the 
      ! material is just exhausted.
      do m = 1,size(vof_n)
        ! volume of material m attempting to leave the cell in this volume_track_subcycle
        total_material_flux = sum(Volume_Flux_Sub(m,:)) 
        if (total_material_flux == 0.0_r8) cycle
        
        ! If the cumulative_outward_material_flux
        ! exceeds the amount of material originally in the cell (from vof_n),
        ! calculate the 'Ratio' of fluid material volume still allowed to be
        ! fluxed to the flux volume, and note that we're not 'Done'
        !if (.not.isImmobile(m)) then
          ! cumulative volume of material m that left the cell in previous subcycles
          cumulative_outward_material_flux = sum(max(volume_flux_tot(m,:),0.0_r8))

          ! ratio between the volume of original (from beginning of flow cycle) material remaining in the cell
          ! (material that has entered is disregarded)
          ! and the volume of material m attempting to leave the cell in this volume track cycle
          ratio = (vof_n(m)*cell_volume - cumulative_outward_material_flux) / total_material_flux
        ! else ! if none of this material was originally in the cell or material isImmobile
        !   ratio = 0.0_r8
        ! end if

        if (Ratio < 1.0_r8) then
          ! lower the fluxes to match the material volume within the cell,
          ! and flag the material number as maxed
          flux_reduced = .true.
          maxed(m) = .true.
          volume_flux_sub(m,:) = ratio * volume_flux_sub(m,:)
        end if
      end do ! material loop
      
      if (flux_reduced) then
        ! This cell had one/more fluxes reduced.  For each of the faces, if the sum
        ! of material fluxes is less than Total_Face_Flux, multiply all non-maxed 
        ! fluxes by another 'Ratio' (this time > 1) that restores the flux balance.  
        ! This may in turn over-exhaust one or more of these materials, and so from
        ! the bottom of this loop, we head back to the top.
        do f = 1,6
          ! Calculate the total flux volume through the cell face (is this already
          ! available elsewhere?), and if the flux volume is greater than zero, 
          ! then concern ourselves with adjusting individual material fluxes.
          if (Total_Face_Flux(f) > cutvof*cell_volume) then
            ! Add up the sum of material fluxes at a face (Sum), and the sum of 
            ! un-maxed material fluxes (Sum_not_maxed).
            total_flux_through_face           = sum(Volume_Flux_Sub(:,f))
            total_flux_through_face_not_maxed = sum(volume_flux_sub(:,f), mask=.not.maxed)

            ! Ratio as defined below, when used to multiply the non-maxed fluxes at 
            ! a face, will restore the flux balance.
            if (total_flux_through_face_not_maxed > 0.0_r8) then
              Ratio = (Total_Face_Flux(f) + total_flux_through_face_not_maxed - total_flux_through_face) &
                   /                       total_flux_through_face_not_maxed
              where (.not.maxed) Volume_Flux_Sub(:,f) = Ratio * Volume_Flux_Sub(:,f)
            else
              number_not_maxed = count(.not.Maxed) ! .and. .not.isImmobile)
              if (number_not_maxed == 0) then
                ierr = 1
                return
              end if

              Ratio = (Total_Face_Flux(f) - total_flux_through_face) / real(number_not_maxed,r8)
              where (.not.Maxed) & ! .and. .not.isImmobile)
                   Volume_Flux_Sub(:,f) = Ratio
            end if

          end if
        end do ! face loop
      else
        exit
      end if
      
    end do ! renorm loop

  end subroutine renorm_cell
  
  !=======================================================================
  ! Purpose(s):
  !
  !   Compute inflow volume fluxes.
  !
  !=======================================================================
  subroutine flux_bc (Fluxing_Velocity, Vof_n, Volume_Flux_Sub, adv_dt, mesh, gmesh)
    use unstr_mesh_type

    type(unstr_mesh), intent(in)    :: mesh
    type(mesh_geom),  intent(in)    :: gmesh
    real(r8),         intent(in)    :: Fluxing_Velocity(:,:), Vof_n(:,:), adv_dt
    real(r8),         intent(inout) :: Volume_Flux_Sub(:,:,:)

    real(r8) :: flux_vol
    integer  :: f,n,fid

    do n = 1,mesh%ncell ! loop over all cells and faces
      do f = 1,6
        fid = mesh%cface(f,n)
        
        ! Calculate the inflow flux volume
        Flux_Vol = min(adv_dt*Fluxing_Velocity(f,n)*mesh%area(fid), 0.0_r8)
        
        ! if this is an inflow face, calculate the flux and delta volume
        if (gmesh%cneighbor(f,n)==-1 .and. flux_vol < -cutvof*mesh%volume(n)) then
          volume_flux_sub(:,f,n) = 0.0_r8 ! Zero out the subcycle volume fluxes

          ! If inflow material specified as a BC, assign it.
          ! otherwise, assume what's flowing in is more of what was in the cell at the beginning of the timestep
          ! if (BC_Mat(fid) /= NULL_I) then
          !   Volume_Flux_Sub(BC_Mat(fid),f,n) = Flux_Vol
          ! else
          Volume_Flux_Sub(:,f,n) = Flux_Vol*Vof_n(:,n)
          ! end if
        end if
      end do
    end do
    
  end subroutine flux_bc
  
  subroutine flux_acceptor (volume_flux_sub,gmesh)
    use unstr_mesh_type

    real(r8), intent(inout) :: volume_flux_sub(:,:,:)
    type(mesh_geom),  intent(in) :: gmesh

    real(r8) :: acceptor_flux
    integer :: m,n,f,nf,nc
    
    do m = 1,size(volume_flux_sub, dim=1)   ! loop through materials
      do n = 1,size(volume_flux_sub, dim=3) ! loop through cells
        do f = 1,6                          ! loop through faces
          nf = gmesh%fneighbor(f,n)
          nc = gmesh%cneighbor(f,n)

          if (nc>0) then ! neighbor exists (not boundary cell)
            acceptor_flux = volume_flux_sub(m,nf,nc)
            if (acceptor_flux > 0.0_r8) volume_flux_sub(m,f,n) = - acceptor_flux
          end if
          
        end do
      end do
    end do
    
    !write(*,*) 'final flux: ',volume_flux_sub(:,:,501)
    
  end subroutine flux_acceptor

  subroutine volume_advance (volume_flux_sub, volume_flux_tot, vof, volume)
    real(r8), intent(in)    :: volume_flux_sub(:,:,:), volume(:)
    real(r8), intent(inout) :: volume_flux_tot(:,:,:), vof(:,:)

    integer :: f,m
    
    volume_flux_tot = volume_flux_tot + volume_flux_sub

    do m = 1,size(vof, dim=1)
      !if (.not.isImmobile(m)) then
        do f = 1,6
          vof(m,:) = vof(m,:) - volume_flux_sub(m,f,:) / volume(:)
        end do
      !end if
    end do
    
  end subroutine volume_advance
  
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
