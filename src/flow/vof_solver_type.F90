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
  implicit none
  private


  type, public :: vof_solver
    private
    type(unstr_mesh), pointer :: mesh => null() ! reference only -- do not own

    !! Pending/current state
    real(r8)                      :: t, dt
    integer, public               :: nmat          ! number of materials present globally
    !real(r8), public, allocatable :: velocity(:,:) ! face velocities
    real(r8), public, allocatable :: fluxing_velocity(:) ! face normal velocities
    real(r8), allocatable :: Volume_flux(:,:,:)
    integer, dimension(:), allocatable, public :: matl_id
    real(r8), dimension(:,:), allocatable, public :: vof
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

contains

  subroutine init (this, mesh, velocity_cc, params)

    use parameter_list_type
    use vof_init

    class(vof_solver), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    real(r8), dimension(:,:), intent(in), target :: velocity_cc
    type(parameter_list) :: params

    integer :: stat
    type(parameter_list), pointer :: plist
    character(:), allocatable :: context,errmsg

    !this%velocity_cc => velocity_cc
    this%mesh        => mesh
    allocate(this%fluxing_velocity(this%mesh%nface))
    ! allocate(this%velocity(3,this%mesh%nface))

    context = 'processing ' // params%name() // ': '

    call params%get ('materials', this%matl_id, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)

    this%nmat = size(this%matl_id)

    allocate(this%vof(this%nmat,this%mesh%ncell))
    allocate(this%volume_flux(this%nmat, 6, this%mesh%ncell))

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



  subroutine advect_mass(this)
    class(vof_solver), intent(inout) :: this

    integer :: stat
    real(r8), dimension(  :,:), allocatable :: Vof_n

    ! Allocate working arrays
    allocate (Vof_n(this%nmat,this%mesh%ncell))

    ! if (volume_track_interfaces .and. .not.associated(VT_Interface_Mask)) then
    !    allocate (VT_Interface_Mask(nfc,ncells), stat = stat)
    !    if (stat /= 0) call LS_fatal ('ADVECT_MASS: VT_Interface_Mask(nfc,ncells) allocation failed')
    ! end if

    ! ! need a copy of the time n Vof values.
    ! Vof_n = this%vof

    ! Advect material volumes.
    call this%advect_volume (Vof_n)

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

    ! ! clean up
    ! deallocate (Vof_n)

    ! ! Stop the Volume Advection Timer
    ! call stop_timer("Mass Advection")

  end subroutine advect_mass

  ! subroutine inflow_outflow_update ()
  !   implicit none

  !   real(r8), dimension(this%mesh%ncell) :: Tmp,mask
  !   integer :: f,stat
  !   integer, parameter :: nfc = 6

  !   ! loop through all faces
  !   do f = 1, nfc
  !      Mask = IN_FLOW (f, this%fluxing_velocity) .OR. OUT_FLOW(f, this%fluxing_velocity)

  !      if (any(Mask)) then
  !         ! Accumulate inflow volume - this%fluxing_velocity < zero
  !         Tmp = MIN(dt*Cell%Face_Area(f)*this%fluxing_velocity(f,:),0.0_r8)
  !         qin = qin + ABS(PGSLib_Global_SUM(Tmp, MASK = Mask))

  !         ! Accumulate outflow volume - this%fluxing_velocity > zero
  !         Tmp = MAX(dt*Cell%Face_Area(f)*this%fluxing_velocity(f,:),0.0_r8)
  !         qout = qout + PGSLib_Global_SUM(Tmp, MASK = Mask)
  !      end if
  !   end do

  ! end subroutine inflow_outflow_update

  subroutine advect_volume (this, vof_n)
    use timer_tree_type
    
    class(vof_solver), intent(inout) :: this
    real(r8),          intent(in)    :: Vof_n(this%nmat,this%mesh%ncell)
    
    integer  :: status
    integer  :: p, vps
    real(r8) :: adv_dt
    real(r8) :: Volume_Flux_Sub(this%nmat, 6, this%mesh%ncell) ! keeps track of volume changes in every subcycle

    ! Start the volume advection timer.
    call start_timer("Volume Tracking")

    ! ! Zero the total flux array.
    ! this%Volume_Flux = 0.0_r8

    ! ! No subcycling if we're not volume tracking.
    ! vps = volume_track_subcycles

    ! ! Set the advection timestep.
    ! adv_dt = dt/real(vps,r8)

    ! do p = 1,vps
    !   ! Get the donor fluxes.
    !    volume_flux_sub = volume_track (this%mesh, this%vof, this%fluxing_velocity, Volume_Flux_Sub)

    !    ! Normalize the donor fluxes.
    !    call this%flux_renorm (Vof_n, Volume_Flux_Sub)

    !    ! Compute the acceptor fluxes.
    !    call flux_acceptor (Volume_Flux_Sub)

    !    ! Compute BC (inflow) fluxes.
    !    call this%flux_bc (Vof_n, Volume_Flux_Sub)

    !    ! Add the volume fluxes from this subcycle (Volume_Flux_Sub) to the
    !    ! total flux array (Volume_Flux_Tot), and update the volume fraction
    !    ! array (Vof).
    !    call this%volume_advance (Volume_Flux_Sub)

    !    ! Make sure volume fractions of a particular material are within
    !    ! the allowed range (0 <= Vof <= 1) and that all materials sum to one.
    !    call this%vof_bounds ()
    ! end do
    
    ! Stop the volume advection timer.
    call stop_timer("Volume Tracking")
    
  end subroutine advect_volume

  ! !=======================================================================
  ! ! Purpose(s):
  ! !
  ! !   Scan all faces with an outward flux and determine if any
  ! !   material is over-exhausted from this cell.  If so lower
  ! !   the fluxes until the material is just exhausted.  Then
  ! !   loop over the faces and balance the individual material
  ! !   fluxes with the total face flux. The sum of the material
  ! !   volume fluxes (Volume_Flux_Sub) for each face should sum
  ! !   to the total volume flux for that face.
  ! !
  ! !=======================================================================
  ! subroutine flux_renorm (Fluxing_Velocity, Vof_n, Volume_Flux_Tot, Volume_Flux_Sub, adv_dt)
  !   use cutoffs_module,       only: cutvof
  !   use fluid_data_module,    only: isImmobile
  !   use mesh_module,          only: Cell
  !   use parameter_module,     only: ncells, nfc, nmat
  !   use vof_data_module,      only: adv_dt

  !   real(r8), intent(in)    :: fluxing_velocity(:,:), vof_n(:,:), adv_dt
  !   real(r8), intent(inout) :: volume_flux_tot(:,:,:), volume_flux_sub(:,:,:)

  !   real(r8) :: Ratio, Sum, Cumul_Sum, Sum_not_maxed, Total_Face_Flux
  !   integer  :: norm_iter, f, m, n, number_not_maxed
  !   logical  :: Done_Renorm, maxed(nmat)

  !   CELLS: do n = 1,ncells

  !     Maxed = .False.

  !     ! Loop over the renorm_loop a maximum of nmat - 1 times, to resolve all instances
  !     ! of fluxing more material than was in the cell at the beginning of the timestep

  !     RENORM_LOOP: do norm_iter = 1, nmat+1

  !       Done_Renorm = .True.

  !       ! We stay in this loop until the sum of material fluxes from each face of
  !       ! the cell equals the total face volume flux.  Where the cumulative sum of 
  !       ! individual material fluxes (from this and previous volume_track_subcycles) 
  !       ! exceeds the volume of a particular material originally within a 
  !       ! cell, we decrease those fluxes to equal the volume of material still available
  !       ! to be fluxed, and increase other fluxes appropriately.  If this increase
  !       ! leads to the over-exhaustion of other materials, we work our way through
  !       ! this loop again, and again, and again, and ... , until we're done.

  !       ! The first step is to determine if any material is being over-exhausted from 
  !       ! a cell.  If so mark it as MAXED and lower the Volume_Flux_Sub's so that the 
  !       ! material is just exhausted.

  !       MAT_LOOP: do m = 1,nmat

  !         ! Remember that at this stage of the code Face_Flux'es are only positive
  !         ! or zero.  Sum is the volume of material m attempting to leave the cell
  !         ! in this volume_track_subcycle; Cumul_Sum is Sum plus the material that
  !         ! has already left in previous volume_track_subcycles.
  !         Sum = 0.0_r8
  !         Cumul_Sum = 0.0_r8
  !         do f = 1,nfc
  !           Sum = Sum + Volume_Flux_Sub(m,f,n)
  !           Cumul_Sum = Cumul_Sum + MAX(Volume_Flux_Tot(m,f,n),0.0_r8)
  !         end do
  !         if (Sum == 0.0_r8) CYCLE MAT_LOOP
  !         Cumul_Sum = Cumul_Sum + Sum

  !         ! If the CUMULATIVE sum of outward fluxes across faces (Cumul_Sum)
  !         ! exceeds the amount of material ORIGINALLY in the cell (from Vof_n),
  !         ! calculate the 'Ratio' of fluid material volume still allowed to be
  !         ! fluxed to the flux volume, and note that we're not 'Done'
  !         Ratio = 0.0_r8  ! if none of this material was originally in the cell
  !         ! Update the Ratio for fluid materials; if the material isImmobile, 
  !         ! leave Ratio = 0.0_r8
  !         if (.not. isImmobile(m)) &
  !              Ratio = (Vof_n(m,n)*Cell(n)%Volume - (Cumul_Sum-Sum)) / Sum
  !         if (Ratio < 1.0_r8) then
  !           Done_Renorm = .False.
  !           Maxed(m) = .True.
  !         end if

  !         ! If Ratio < 1, lower the fluxes to match the material volume within
  !         ! the cell, and flag the cell and material number with 'Maxed'.
  !         if (Ratio < 1.0_r8) then
  !           do f = 1,nfc
  !             Volume_Flux_Sub(m,f,n) = Ratio * Volume_Flux_Sub(m,f,n)
  !           end do
  !         end if

  !       end do MAT_LOOP

  !       if (Done_Renorm) exit RENORM_LOOP

  !       ! This cell had one/more fluxes reduced.  For each of the faces, if the sum
  !       ! of material fluxes is less than Total_Face_Flux, multiply all non-maxed 
  !       ! fluxes by another 'Ratio' (this time > 1) that restores the flux balance.  
  !       ! This may in turn over-exhaust one or more of these materials, and so from
  !       ! the bottom of this loop, we head back to the top.

  !       do f = 1,nfc

  !         ! Calculate the total flux volume through the cell face (is this already
  !         ! available elsewhere?), and if the flux volume is greater than zero, 
  !         ! then concern ourselves with adjusting individual material fluxes.

  !         Total_Face_Flux = adv_dt*Fluxing_Velocity(f,n)*Cell(n)%Face_Area(f)
  !         if (Total_Face_Flux > cutvof*Cell(n)%Volume) then

  !           ! Add up the sum of material fluxes at a face (Sum), and the sum of 
  !           ! un-maxed material fluxes (Sum_not_maxed).
  !           Sum = 0.0_r8
  !           Sum_not_maxed = 0.0_r8
  !           do m = 1,nmat
  !             Sum = Sum + Volume_Flux_Sub(m,f,n)
  !             if (.not. Maxed(m)) Sum_not_maxed = Sum_not_maxed + Volume_Flux_Sub(m,f,n)
  !           end do

  !           ! Ratio as defined below, when used to multiply the non-maxed fluxes at 
  !           ! a face, will restore the flux balance.
  !           if (Sum_not_maxed > 0.0_r8) then
  !             ! jms Note:  Ratio = (Total_Face_Flux - Maxed_Face_Flux) / Sum_not_maxed
  !             Ratio = 1.0_r8 + (Total_Face_Flux - Sum) / Sum_not_maxed
  !             do m = 1,nmat
  !               if (.not. Maxed(m)) Volume_Flux_Sub(m,f,n) = Ratio * Volume_Flux_Sub(m,f,n)
  !             end do
  !           else
  !             number_not_maxed = 0
  !             do m = 1,nmat
  !               if (.not. Maxed(m) .and. .not.isImmobile(m)) number_not_maxed = number_not_maxed + 1
  !             end do
  !             if (number_not_maxed == 0) &
  !                  call LS_fatal ('FLUX_RENORM: cannot reassign face flux to any other material')
              
  !             Ratio = (Total_Face_Flux - Sum) / number_not_maxed
  !             do m = 1,nmat
  !               if (.not. Maxed(m).and. .not.isImmobile(m)) Volume_Flux_Sub(m,f,n) = Ratio
  !             end do
  !           end if

  !         end if
  !       end do ! face loop
        
  !     end do RENORM_LOOP
  !   end do CELLS

  ! END SUBROUTINE FLUX_RENORM


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
