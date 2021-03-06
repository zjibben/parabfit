!!
!! FLOW_SOLVER_TYPE
!!
!! This module defines a class that encapsulates a (time) solver for the
!! discrete heat conduction model.  It is based on IDAESOL which is an
!! integrator for index-1 DAE systems that uses BDF2 time discretization.
!!
!! TODO: * make the fluxing velocity point to the array in NS_solver_type
!!       * make some sort of exact solution type and do away with the vof0
!!         array (currently only used for the deforming sphere test)
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

#include "f90_assert.fpp"

module vof_solver_type

  use kinds,  only: r8
  use consts, only: cutvof,nfc
  use unstr_mesh_type
  use mesh_geom_type
  use logging_services
  use surface_type
  implicit none
  private

  type, public :: vof_solver
    private

    type(unstr_mesh), pointer, public :: mesh          => null() ! reference only -- do not own
    type(mesh_geom),  pointer, public :: gmesh         => null() ! reference only -- do not own

    real(r8),         pointer :: fluidRho(:)           => null() ! reference only -- do not own
    real(r8),         pointer :: fluxing_velocity(:,:) => null() ! reference only -- do not own

    !! Pending/current state
    integer,                    public :: advect_method, volume_track_subcycles
    real(r8)                           :: t, dt
    integer,                    public :: nmat                 ! number of materials present globally
    !real(r8),      allocatable, public :: velocity(:,:)       ! face velocities
    integer,       allocatable, public :: matl_id(:)
    real(r8),      allocatable, public :: vof(:,:), vof0(:,:), volume_flux_tot(:,:,:)
    type(surface), allocatable, public :: intrec(:)             ! interface reconstruction
  contains
    procedure :: init
    procedure :: set_initial_state
    ! procedure :: test_initial_state
    procedure :: advect_mass
    procedure, private :: advect_volume
    procedure :: update_intrec_surf
    procedure :: print_vofs
    procedure :: print_error_norms
    !procedure, private :: flux_renorm
  end type vof_solver

  ! TODO: this needs to be removed and replaced with use of the velocity_bc and pressure_bc variables
  integer, allocatable :: BC_flag(:,:)

  logical, parameter :: using_mic = .false.

  ! local warning message flags
  integer,  save :: WLimit     = 10
  integer,  save :: WCountTot  = 0
  real(r8), save :: WMaxTot    = 0.0_r8
  integer,  save :: WCountMat  = 0
  real(r8), save :: WMaxMat    = 0.0_r8
  integer,  save :: WCountTotU = 0
  real(r8), save :: WMaxTotU   = 0.0_r8
  integer,  save :: WCountMatU = 0
  real(r8), save :: WMaxMatU   = 0.0_r8

  ! advection methods
  integer, parameter, public :: ONION_SKIN        = 1
  integer, parameter, public :: NESTED_DISSECTION = 2

contains

  subroutine init (this, mesh, gmesh, nmat, fluxing_velocity, fluidRho, params, velocity_bc)

    use consts, only: nfc
    use parameter_list_type
    use bndry_func_class

    class(vof_solver),        intent(out) :: this
    type(unstr_mesh), target, intent(in)  :: mesh
    type(mesh_geom),  target, intent(in)  :: gmesh
    ! WARNING: these aren't actually targets. need to pass in when they are used, not store them
    !          as type members
    real(r8),         target, intent(in)  :: fluxing_velocity(:,:), fluidRho(:)
    integer,                  intent(in)  :: nmat
    class(bndry_func), allocatable, intent(in)  :: velocity_bc
    type(parameter_list)                  :: params

    integer                       :: stat
    type(parameter_list), pointer :: plist
    character(:),     allocatable :: context,errmsg

    this%mesh  => mesh
    this%gmesh => gmesh

    this%fluidRho => fluidRho
    this%fluxing_velocity => fluxing_velocity

    context = 'processing ' // params%name() // ': '

    call params%get ('materials', this%matl_id, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)

    call params%get ('advection-method', this%advect_method, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)

    call params%get ('subcycles', this%volume_track_subcycles, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)

    ! when we prescribe the velocity field (and hence don't need to specify material properties),
    ! the nmat passed in will not have been initialized
    this%nmat = size(this%matl_id) !nmat

    allocate(this%vof(this%nmat,this%mesh%ncell), this%vof0(this%nmat,this%mesh%ncell), &
        this%volume_flux_tot(this%nmat,nfc,this%mesh%ncell), &
        this%intrec(this%nmat))

    this%vof = 0.0_r8
    this%volume_flux_tot = 0.0_r8

    ! this type of boundary handling ought to be refactored (it is already done nicely elsewhere)
    ! for right now, this just exists so that we can run simple tests,
    ! without removing the branches of code that deal with more complex scenarios
    allocate(BC_flag(nfc,this%mesh%ncell))
    BC_flag = 0

    ! when we prescribe the velocity field, velocity_bc will be unallocated
    if (allocated(velocity_bc)) call assign_vel_bc_flags (BC_flag, velocity_bc)

  contains

    subroutine assign_vel_bc_flags (BC_flag, velocity_bc)

      integer,           intent(inout) :: BC_flag(:,:)
      class(bndry_func), intent(in)    :: velocity_bc

      integer :: bndry_f, i, f, fid

      do bndry_f = 1,size(velocity_bc%index)
        fid = velocity_bc%index(bndry_f)

        i = gmesh%fcell(1,fid) ! id of cell attached to this face
        f = gmesh%flid(1,fid)  ! local id of the face

        BC_flag(f,i) = 2
      end do

    end subroutine assign_vel_bc_flags

  end subroutine init

  subroutine set_initial_state (this, plist)

    use parameter_list_type
    use vof_init

    class(vof_solver),    intent(inout) :: this
    type(parameter_list), intent(in)    :: plist

    !! Initialize the Vof
    call vof_initialize (this%mesh, this%gmesh, plist, this%vof, this%matl_id, this%nmat)
    this%vof0 = this%vof

    ! DEBUGGING ############
    call this%print_vofs ()
    ! ######################

  end subroutine set_initial_state

  subroutine advect_mass (this, dt, matl_is_void, dump_intrec)

    use timer_tree_type

    class(vof_solver), intent(inout) :: this
    real(r8),          intent(in)    :: dt
    logical,           intent(in)    :: matl_is_void(:), dump_intrec

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
    call this%advect_volume (Vof_n, dt, matl_is_void, dump_intrec)

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

  subroutine advect_volume (this, vof_n, dt, matl_is_void, dump_intrec)

    use consts,                 only: nfc
    use timer_tree_type
    use volume_track_module,    only: volume_track
    use volume_track_nd_module, only: volume_track_nd

    class(vof_solver), intent(inout) :: this
    real(r8),          intent(in)    :: Vof_n(:,:), dt
    logical,           intent(in)    :: matl_is_void(:), dump_intrec

    integer  :: p
    real(r8) :: adv_dt
    real(r8) :: volume_flux_sub(this%nmat, nfc, this%mesh%ncell) ! volume changes in every subcycle

    call start_timer ("volume tracking")

    this%volume_flux_tot = 0.0_r8                    ! Zero the total flux array
    adv_dt = dt/real(this%volume_track_subcycles,r8) ! Set the advection timestep

!!$omp parallel if(using_mic) default(private) shared(adv_dt, this, vof_n, volume_flux_sub, volume_track_subcycles)
    do p = 1,this%volume_track_subcycles
      ! Get the donor fluxes.
      select case (this%advect_method)
      case (ONION_SKIN)
        call volume_track (volume_flux_sub, adv_dt, this%mesh, this%gmesh, this%vof, &
            this%fluxing_velocity, this%nmat, this%fluidRho, this%intrec, &
            dump_intrec .and. p==this%volume_track_subcycles)
      case (NESTED_DISSECTION)
        call volume_track_nd (volume_flux_sub, adv_dt, this%mesh, this%gmesh, this%vof, &
            this%fluxing_velocity, this%fluidRho, this%intrec, &
            dump_intrec .and. p==this%volume_track_subcycles)
      case default
        call LS_fatal ("invalid advection method")
      end select

      !$omp parallel if(.not.using_mic) default(private) &
      !$omp shared(adv_dt, this, vof_n, volume_flux_sub, matl_is_void)
      ! Normalize the donor fluxes.
      call flux_renorm (this%fluxing_velocity, Vof_n, this%volume_flux_tot, Volume_Flux_Sub, &
          adv_dt, this%mesh)

      ! Compute the acceptor fluxes.
      call flux_acceptor (volume_flux_sub, this%gmesh)

      ! Compute BC (inflow) fluxes.
      call flux_bc (this%fluxing_velocity, vof_n, volume_flux_sub, adv_dt, this%mesh, this%gmesh)

      ! Add the volume fluxes from this subcycle to the total flux array
      ! and update the volume fraction array
      call volume_advance (Volume_Flux_Sub, this%volume_flux_tot, this%vof, this%mesh%volume)

      ! Ensure volume fractions of each material are within 0 and 1 and that all materials sum to one
      call vof_bounds (this%vof, this%volume_flux_tot, matl_is_void, this%mesh)
      !$omp end parallel
    end do
!!$omp end parallel
    call stop_timer ("volume tracking")

  end subroutine advect_volume

  ! Scan all faces with an outward flux and determine if any
  ! material is over-exhausted from this cell. If so lower
  ! the fluxes until the material is just exhausted. Then
  ! loop over the faces and balance the individual material
  ! fluxes with the total face flux. The sum of the material
  ! volume fluxes (Volume_Flux_Sub) for each face should sum
  ! to the total volume flux for that face.
  subroutine flux_renorm (Fluxing_Velocity, Vof_n, Volume_Flux_Tot, Volume_Flux_Sub, adv_dt, mesh)

    real(r8),          intent(in)    :: fluxing_velocity(:,:), vof_n(:,:), adv_dt, &
        volume_flux_tot(:,:,:)
    real(r8),          intent(inout) :: volume_flux_sub(:,:,:)
    class(unstr_mesh), intent(in)    :: mesh

    integer :: n,ierr, j
    real(r8), allocatable :: volume_flux_sub_cell(:,:)

    !$omp do
    do n = 1,mesh%ncell
      volume_flux_sub_cell = volume_flux_sub(:,:,n)
      call renorm_cell (volume_flux_sub(:,:,n), volume_flux_tot(:,:,n), vof_n(:,n), &
          adv_dt*Fluxing_Velocity(:,n)*mesh%area(mesh%cface(:,n)), mesh%volume(n), ierr)
      if (ierr /= 0) then
        write(*,'(a,i10)')     'cell id:       ',n
        write(*,'(a,3f14.10)') 'cell centroid: ',sum(mesh%x(:,mesh%cnode(:,n)), dim=2) / 8.0_r8
        write(*,'(a,3es14.4)') 'vof:           ',vof_n(:,n)*mesh%volume(n)
        write(*,'(a,3es14.4)') 'tot matl flx:  ',sum(volume_flux_sub_cell, dim=2)
        write(*,'(a,3es14.4)') 'prv matl flx:  ',sum(volume_flux_tot(:,:,n), dim=2)
        write(*,'(a,3es14.4)') 'tot flx:       ',sum(adv_dt*Fluxing_Velocity(:,n)*mesh%area(mesh%cface(:,n)))
        do j = 1,nfc
          write(*,'(a,3es14.4)')  'material flux: ',volume_flux_sub(:,j,n)
        end do
        write(*,*)
        do j = 1,nfc
          write(*,'(a,2es14.4)')  'volume flux:   ',sum(volume_flux_sub_cell(:,j)), &
              adv_dt*mesh%area(mesh%cface(j,n))*fluxing_velocity(j,n)
        end do
        call LS_fatal ('FLUX_RENORM: cannot reassign face flux to any other material')
      end if
    end do
    !$omp end do nowait

  end subroutine flux_renorm


  ! renorm_cell: ensure no materials are over-exhausted in a given cell
  !
  ! note 1:   We stay in this loop until the sum of material fluxes from each face of
  !           the cell equals the total face volume flux.  Where the cumulative sum of
  !           individual material fluxes (from this and previous volume_track_subcycles)
  !           exceeds the volume of a particular material originally within a
  !           cell, we decrease those fluxes to equal the volume of material still available
  !           to be fluxed, and increase other fluxes appropriately.  If this increase
  !           leads to the over-exhaustion of other materials, we work our way through
  !           this loop again, and again, and again, and ... , until we're done.
  subroutine renorm_cell (volume_flux_sub, volume_flux_tot, vof_n, total_face_flux, cell_volume, &
      ierr)

    use array_utils, only: isZero

    real(r8), intent(inout) :: volume_flux_sub(:,:)
    real(r8), intent(in)    :: volume_flux_tot(:,:), vof_n(:), total_face_flux(:), cell_volume
    integer,  intent(out)   :: ierr

    real(r8) :: Ratio, total_material_flux, cumulative_outward_material_flux, &
        total_flux_through_face, total_flux_through_face_not_maxed
    integer  :: norm_iter, f, m, number_not_maxed
    logical  :: flux_reduced, maxed(size(vof_n))

    maxed = .false.
    ierr = 0

    ! Loop over the renorm_loop a maximum of nmat - 1 times, to resolve all instances
    ! of fluxing more material than was in the cell at the beginning of the timestep
    do norm_iter = 1, size(vof_n)+1
      flux_reduced = .false.

      ! see note 1

      ! The first step is to determine if any material is being over-exhausted from
      ! a cell.  If so mark it as MAXED and lower the Volume_Flux_Sub's so that the
      ! material is just exhausted.
      do m = 1,size(vof_n)
        ! volume of material m attempting to leave the cell in this volume_track_subcycle
        total_material_flux = sum(Volume_Flux_Sub(m,:))
        if (isZero(total_material_flux)) cycle

        ! If the cumulative_outward_material_flux exceeds the amount of material
        ! originally in the cell (from vof_n), calculate the 'Ratio' of fluid material volume
        ! still allowed to be fluxed to the flux volume, and note that we're not 'Done'
        !if (.not.isImmobile(m)) then
        ! cumulative volume of material m that left the cell in previous subcycles
        cumulative_outward_material_flux = sum(volume_flux_tot(m,:), &
            mask=volume_flux_tot(m,:) > 0.0_r8)

        ! ratio between the volume of original (from beginning of flow cycle) material remaining in the cell
        ! (material that has entered is disregarded)
        ! and the volume of material m attempting to leave the cell in this volume track cycle
        ratio = (vof_n(m)*cell_volume - cumulative_outward_material_flux) / total_material_flux
        ! else ! if material isImmobile
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

      if (.not.flux_reduced) exit

      ! This cell had one/more fluxes reduced.  For each of the faces, if the sum
      ! of material fluxes is less than Total_Face_Flux, multiply all non-maxed
      ! fluxes by another 'Ratio' (this time > 1) that restores the flux balance.
      ! This may in turn over-exhaust one or more of these materials, and so from
      ! the bottom of this loop, we head back to the top.
      do f = 1,nfc
        ! Calculate the total flux volume through the cell face (is this already
        ! available elsewhere?), and if the flux volume is greater than zero,
        ! then concern ourselves with adjusting individual material fluxes.
        if (Total_Face_Flux(f) > cutvof*cell_volume) then
          ! calculate the sum of material fluxes at a face, and the sum of un-maxed material fluxes.
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

            where (.not.Maxed) & ! .and. .not.isImmobile)
                Volume_Flux_Sub(:,f) = (Total_Face_Flux(f) - total_flux_through_face) &
                / real(number_not_maxed,r8)
          end if

        end if
      end do ! face loop

    end do ! renorm loop

  end subroutine renorm_cell

  ! Compute inflow volume fluxes.
  subroutine flux_bc (Fluxing_Velocity, Vof_n, Volume_Flux_Sub, adv_dt, mesh, gmesh)

    type(unstr_mesh), intent(in)    :: mesh
    type(mesh_geom),  intent(in)    :: gmesh
    real(r8),         intent(in)    :: Fluxing_Velocity(:,:), Vof_n(:,:), adv_dt
    real(r8),         intent(inout) :: Volume_Flux_Sub(:,:,:)

    real(r8) :: flux_vol
    integer  :: f,n,fid

    !$omp do
    do n = 1,mesh%ncell ! loop over all cells and faces
      do f = 1,nfc
        fid = mesh%cface(f,n)

        ! Calculate the inflow flux volume
        Flux_Vol = min(adv_dt*Fluxing_Velocity(f,n)*mesh%area(fid), 0.0_r8)

        ! if this is an inflow face, calculate the flux and delta volume
        if (gmesh%cneighbor(f,n)==-1 .and. flux_vol < -cutvof*mesh%volume(n)) then

          ! If inflow material specified as a BC, assign it.
          ! otherwise, assume what's flowing in is more of what was in the cell at the beginning of the timestep
          ! if (BC_Mat(fid) /= NULL_I) then
          !   volume_flux_sub(:,f,n) = 0.0_r8 ! Zero out the subcycle volume fluxes
          !   Volume_Flux_Sub(BC_Mat(fid),f,n) = Flux_Vol
          ! else
          Volume_Flux_Sub(:,f,n) = Flux_Vol*Vof_n(:,n)
          ! end if
        end if
      end do
    end do
    !$omp end do nowait

  end subroutine flux_bc

  subroutine flux_acceptor (volume_flux_sub,gmesh)

    real(r8), intent(inout) :: volume_flux_sub(:,:,:)
    type(mesh_geom),  intent(in) :: gmesh

    real(r8) :: acceptor_flux
    integer :: m,n,f,nf,nc

    !$omp barrier

    !$omp do
    do n = 1,size(volume_flux_sub, dim=3)
      do f = 1,nfc
        nc = gmesh%cneighbor(f,n)   ! neighbor cell id

        if (nc>0) then              ! neighbor exists (not boundary cell)
          nf = gmesh%fneighbor(f,n) ! neighbor cell's local face id
          ! check how much of every material is being fluxed to me from my neighbor,
          ! if there is any, note it in my own volume_flux_sub
          do m = 1,size(volume_flux_sub, dim=1)
            acceptor_flux = volume_flux_sub(m,nf,nc)
            if (acceptor_flux > 0.0_r8) volume_flux_sub(m,f,n) = -acceptor_flux
          end do
        end if

      end do
    end do
    !$omp end do

  end subroutine flux_acceptor

  subroutine volume_advance (volume_flux_sub, volume_flux_tot, vof, volume)
    real(r8), intent(in)    :: volume_flux_sub(:,:,:), volume(:)
    real(r8), intent(inout) :: volume_flux_tot(:,:,:), vof(:,:)

    integer :: f,m,n

    !$omp do
    do n = 1,size(vof, dim=2)

      volume_flux_tot(:,:,n) = volume_flux_tot(:,:,n) + volume_flux_sub(:,:,n)

      do m = 1,size(vof, dim=1)
        !if (.not.isImmobile(m)) &
        vof(m,n) = vof(m,n) - sum(volume_flux_sub(m,:,n)) / volume(n)
      end do
    end do
    !$omp end do nowait

  end subroutine volume_advance

  ! Make sure volume fractions are within bounds:  0 <= Vof <= 1.
  ! If not, remove the overshoots (Vof > 1) and undershoots (Vof < 0).
  subroutine vof_bounds (vof, volume_flux_tot, matl_is_void, mesh)

    use array_utils, only: last_true_loc

    real(r8),         intent(inout) :: vof(:,:),Volume_Flux_Tot(:,:,:)
    logical,          intent(in)    :: matl_is_void(:)
    type(unstr_mesh), intent(in)    :: mesh

    integer  :: n, nmat

    nmat = size(vof,dim=1)

    !$omp do
    do n = 1,mesh%ncell
      call vof_bounds_cell (vof(:,n), volume_flux_tot(:,:,n), matl_is_void, mesh%volume(n), &
          nmat, BC_flag(:,n))
    end do
    !$omp end do nowait

  end subroutine vof_bounds

  subroutine vof_bounds_cell (vof, volume_flux_tot, matl_is_void, volume, nmat, BC_flag)

    use array_utils, only: last_true_loc

    real(r8),         intent(inout) :: vof(:), volume_flux_tot(:,:)
    real(r8),         intent(in)    :: volume
    logical,          intent(in)    :: matl_is_void(:)
    integer,          intent(in)    :: nmat, BC_flag(:)

    integer  :: m, void_m
    real(r8) :: Ftot, Ftot_m1, void_volume, delta_vol(size(vof,dim=1))

    ! make sure individual material vofs are within bounds
    do m = 1,nmat
      if (vof(m) > 1.0_r8-cutvof .and. vof(m) /= 1.0_r8) then
        ! if volume fraction is close to 1, round to 1
        if (matl_is_void(m)) then
          vof(m) = 1.0_r8
        else
          call adjust_flux_matl (volume_flux_tot(m,:), vof(m), 1.0_r8, volume, BC_flag)
          vof(m) = 1.0_r8
        end if
      else if (vof(m) < cutvof .and. vof(m) /= 0.0_r8) then
        ! if volume fraction is close to 0, round to 0
        if (matl_is_void(m)) then
          vof(m) = 0.0_r8
        else if (vof(m) /= 0.0_r8) then
          call adjust_flux_matl (volume_flux_tot(m,:), vof(m), 0.0_r8, volume, BC_flag)
          vof(m) = 0.0_r8
        end if
      end if
    end do

    ! make sure the sum of vofs is within bounds
    Ftot = sum(vof)
    if (Ftot == 1.0_r8) return

    ! Renormalize the liquid volume fractions

    ! Check to see if void is already in the cell.
    ! if so, grab the index of the last void material present
    void_m = last_true_loc (matl_is_void .and. vof > 0.0_r8)
    if (any(matl_is_void) .and. void_m > 0) then
      ! sum up the volume fraction taken up by void
      void_volume = sum(vof, matl_is_void .and. vof > 0.0_r8) ! can the vof be <0?
      if (Ftot > 1.0_r8) then ! if there is too much material in this cell
        ! Is there enough to balance the cell?
        if (void_volume > Ftot-1.0_r8) then ! There is enough void ...
          Ftot_m1 = Ftot - 1.0_r8
          do m = 1,nmat
            if (.not.matl_is_void(m)) cycle
            if (Vof(m) > Ftot_m1) then
              Vof(m) = Vof(m) - Ftot_m1
              exit
            else
              Ftot_m1 = Ftot_m1 - Vof(m)
              Vof(m) = 0.0_r8
            end if
          end do
        else ! There isn't enough void ...
          do m = 1,nmat
            if (matl_is_void(m)) then
              Ftot = Ftot - Vof(m)
              Vof(m) = 0.0_r8
            end if
          end do
          call adjust_flux_total (volume_flux_tot, delta_vol, Ftot, volume, BC_flag, matl_is_void)
          vof = vof + delta_vol / volume
        end if
      else ! Ftot < 1, and there's void already in the cell
        vof(void_m) = vof(void_m) + 1.0_r8 - Ftot
      end if

    else ! there's no void in this cell
      call adjust_flux_total (volume_flux_tot, delta_vol, Ftot, volume, BC_flag, matl_is_void)
      vof = vof + delta_vol / volume
    end if

  end subroutine vof_bounds_cell

  ! Adjust the material volume fluxes on the faces of a single cell
  ! to match the evaluated material volume to a target value.
  subroutine adjust_flux_matl (Volume_Flux_Tot, Current_Material_Vof, Target_Material_Vof, volume, &
      local_BC)

    real(r8), intent(inout) :: Volume_Flux_Tot(:)
    real(r8), intent(in)    :: Current_Material_Vof,Target_Material_Vof,volume
    integer,  intent(in)    :: local_BC(:)

    integer        :: f
    real(r8)       :: Total_Flow, Change_Fraction
    character(128) :: message

    ! Calculate the fractional change needed to adjust the current volume
    ! to the target value.  The same fractional increase/decrease is applied
    ! to incoming and outgoing flows.
    Total_Flow      = inflow  (volume_flux_tot, local_BC) + outflow (volume_flux_tot, local_BC)
    Change_Fraction = Target_Material_Vof - Current_Material_Vof

    if (Total_Flow /= 0.0_r8) then
      ! calculate the ratio between change in volume flux needed and current total volume being fluxed
      Change_Fraction = Change_Fraction*volume / Total_Flow

      ! This needs to be fixed for OpenMP
      ! ! Warning Message if the fractional change in volume is excessive
      ! if (abs(Change_Fraction) > WMaxMat .and. Total_Flow > 1.0e-4*volume) then
      !   WCountMat = 0
      !   WMaxMat   = abs(Change_Fraction)
      !   write(message,'(a,es11.3)') 'excessive volume adjustment: adjustment fraction=', Change_Fraction !' !: cell=', n, &
      !   !', matid=', MatId, '
      !   call LS_warn (message)
      ! else if (WCountMat < Wlimit .and. Total_Flow > 1.0e-4*volume) then
      !   WCountMat = WCountMat + 1
      !   write(message,'(a,es11.3)') 'excessive volume adjustment: adjustment fraction=', Change_Fraction !' !: cell=', n, &
      !   !', matid=', MatId, '
      !   call LS_warn (message)
      ! end if

      ! Adjust the volume changes.
      do f = 1,nfc
        ! Don't change dirichlet velocity BCs.
        if (local_BC(f)/=2) volume_flux_tot(f) = &
            (1.0_r8-sign(1.0_r8,volume_flux_tot(f))*Change_Fraction)*Volume_Flux_Tot(f)
      end do

    else ! there is some flow for this material
      if (Change_Fraction < 0.0_r8) then
        ! jms Note:   If the material is to be removed from the cell
        ! look for a face that doesn't have incoming material, and
        ! flux it out through that face
        do f = 1,nfc
          if (local_BC(f)/=2 .and. Volume_Flux_Tot(f) >= 0.0_r8) then
            volume_flux_tot(f) = volume_flux_tot(f) - change_fraction*volume
            exit
          end if
        end do
      else
        ! This needs to be fixed for OpenMP
        ! if (abs(change_fraction) > WMaxMatU) then
        !   WCountMatU = 0
        !   WMaxMatU = abs(change_fraction)
        !   write(message,'(2(a,i0),a,es11.3)') 'unable to adjust material volume' !: cell=', n, &
        !   !', matid=', MatId, ', desired change fraction=', Change_Fraction
        !   call LS_warn (message)
        ! else if (WCountMatU < Wlimit) then
        !   WCountMatU = WCountMatU + 1
        !   write(message,'(2(a,i0),a,es11.3)') 'unable to adjust material volume' !: cell=', n, &
        !   !', matid=', MatId, ', desired change fraction=', Change_Fraction
        !   call LS_warn (message)
        ! end if
      end if
    end if

  end subroutine adjust_flux_matl

  ! Adjust the material fluxes on the faces of a single cell to match
  ! the evaluated total material volume to a target value
  subroutine adjust_flux_total (volume_flux_tot, delta_vol, current_vof, volume, local_BC, &
      matl_is_void)

    real(r8), intent(inout) :: Volume_Flux_Tot(:,:)
    real(r8), intent(out)   :: Delta_Vol(:)
    real(r8), intent(in)    :: Current_Vof, volume
    integer,  intent(in)    :: local_BC(:)
    logical,  intent(in)    :: matl_is_void(:)

    integer        :: f, m, v, nmat
    real(r8)       :: Inflow_Volume, Outflow_Volume, Volume_Change, Total_Flow, Change_Fraction
    character(128) :: message

    delta_vol = 0.0_r8
    nmat = size(volume_flux_tot,dim=1)

    ! Determine total incoming and outgoing volumes changes
    Inflow_Volume = 0.0_r8; Outflow_Volume = 0.0_r8
    do m = 1, nmat
      if (.not.matl_is_void(m)) then
        inflow_volume  = inflow_volume  + inflow  (volume_flux_tot(m,:), local_BC)
        outflow_volume = outflow_volume + outflow (volume_flux_tot(m,:), local_BC)
      end if
    end do

    ! Calculate the fractional change needed to adjust the current volume
    ! to the target value.  The same fractional increase/decrease is applied
    ! to incoming and outgoing flows.
    Total_Flow = Inflow_Volume + Outflow_Volume
    Volume_Change = 1.0_r8 - Current_Vof
    if (Total_Flow /= 0.0_r8) then
      Change_Fraction = Volume_Change*volume/Total_Flow

      ! This needs to be fixed for OpenMP
      ! ! Warning Message if the fractional change in volume is excessive
      ! if (abs(Change_Fraction) > WMaxTot .and. Total_Flow > 1.0e-4*volume) then
      !   WCountTot = 0
      !   WMaxTot = ABS(Change_Fraction)
      !   write(message,'(a,es11.3)') 'excessive volume adjustment: adjustment fraction=', Change_Fraction !: cell=', n, &
      !   call LS_warn (message)
      ! else if (WCountTot < Wlimit .and. Total_Flow > 1.0e-4*volume) then
      !   WCountTot = WCountTot + 1
      !   write(message,'(a,es11.3)') 'excessive volume adjustment: adjustment fraction=', Change_Fraction !' !: cell=', n, &
      !   call LS_warn (message)
      ! endif

      ! Adjust the Volume Changes by face and accumulate them by cell.
      Delta_Vol = 0.0_r8
      do m = 1, nmat
        if (matl_is_void(m)) cycle
        do f = 1,nfc
          if (local_BC(f)==2) cycle ! Don't change dirichlet velocity BCs.
          ! adjust material transfers
          Delta_Vol(m) = Delta_Vol(m) + Change_Fraction*abs(Volume_Flux_Tot(m,f))
          Volume_Flux_Tot(m,f) = Volume_Flux_Tot(m,f) &
              * (1.0_r8 - sign(1.0_r8,volume_flux_tot(m,f))*Change_Fraction)
        end do
      end do ! material loop
    else ! there is no total flow through this cell
      ! This needs to be fixed for OpenMP
      ! if (abs(Volume_Change) > WMaxTotU) then
      !   WCountTotU = 0
      !   WMaxTotU   = abs(Volume_Change)
      ! else if (WCountTotU < Wlimit) then
      !   WCountTotU = WCountTotU + 1
      ! end if
      ! write(message,'(a,i0,a,es11.3)') 'unable to adjust total volume' !: cell=', n, &
      ! !', desired change fraction=', Volume_Change
      ! call LS_warn (message)
    end if

  end subroutine adjust_flux_total

  real(r8) function outflow (volume_flux_tot, local_BC)
    real(r8), intent(in) :: volume_flux_tot(:)
    integer,  intent(in) :: local_BC(:)
    outflow = sum(volume_flux_tot, mask = local_bc/=2 .and. volume_flux_tot>0.0_r8)
  end function outflow

  real(r8) function inflow (volume_flux_tot, local_BC)
    real(r8), intent(in) :: volume_flux_tot(:)
    integer,  intent(in) :: local_BC(:)
    inflow = -sum(volume_flux_tot, mask = local_bc/=2 .and. volume_flux_tot<0.0_r8)
  end function inflow

  subroutine update_intrec_surf (this)

    use locate_plane_module
    use polyhedron_type
    use multimat_cell_type
    use int_norm_module, only: interface_normal
    use hex_types,       only: hex_f, hex_e

    class(vof_solver), intent(inout) :: this

    real(r8)               :: int_norm(3,this%nmat,this%mesh%ncell), vofint
    integer                :: i,ni,locate_plane_niters, ierr
    type(locate_plane_hex) :: plane_cell
    type(multimat_cell)    :: ndcell
    type(polyhedron)       :: poly

    do ni = 1,size(this%intrec)
      call this%intrec(ni)%purge ()
    end do

    ! compute interface normal vectors for all the materials.
    int_norm = interface_normal (this%vof, this%mesh, this%gmesh, this%advect_method==ONION_SKIN)

    do i = 1,this%mesh%ncell
      select case (this%advect_method)
      case (ONION_SKIN)
        do ni = 1,size(this%intrec)
          vofint = min(max(sum(this%vof(1:ni,i)), 0.0_r8), 1.0_r8)
          if (cutvof < vofint .and. vofint < 1.0_r8 - cutvof) then
            call plane_cell%init (int_norm(:,ni,i), vofint, this%mesh%volume(i), &
                this%mesh%x(:,this%mesh%cnode(:,i)))
            call plane_cell%locate_plane (locate_plane_niters)
            call poly%init (ierr, i, this%mesh, this%gmesh)
            !call poly%init (ierr, plane_cell%node, hex_f, hex_e)
            if (ierr /= 0) call LS_fatal ('plane reconstruction dump failed')
            !call this%intrec(ni)%append (poly%intersection_verts (plane_cell%P), i)
          end if
        end do
      case (NESTED_DISSECTION)
        ! send cell data to the multimat_cell type
        call ndcell%init (ierr, i, this%mesh, this%gmesh)
        if (ierr /= 0) call LS_fatal ('plane reconstruction dump failed')

        ! partition the cell based on vofs and norms
        call ndcell%partition (this%vof(:,i), int_norm(:,:,i))

        ! dump the interface reconstruction, if requested
        do ni = 1,size(this%intrec)
          call this%intrec(ni)%append (ndcell%interface_polygons(ni), i)
        end do
      case default
        call LS_fatal('unknown partition type')
      end select
    end do

  end subroutine update_intrec_surf

  subroutine print_vofs (this)

    class(vof_solver), intent(in) :: this

    integer :: n

    do n = 1,size(this%vof, dim=1)
      write(*,'(a,i2,es20.10)') 'vof ', n, sum(this%vof(n,:)*this%mesh%volume)
    end do

  end subroutine print_vofs

  subroutine print_error_norms (this)

    class(vof_solver), intent(in) :: this

    real(r8) :: error(this%mesh%ncell)

    ! calculate the error
    ! TODO: compare against a given exact solution.
    !       right now just initial condition for time-reversing tests
    error = sum(abs(this%vof - this%vof0), dim=1) / size(this%vof, dim=1)

    write(*,'(a,3es14.4)') 'L1, L2, Linf norm: ', l1_norm(error), l2_norm(error), linf_norm(error)

  end subroutine print_error_norms

  ! TODO: these norms are rather general and could be placed elsewhere
  real(r8) function l1_norm (error)
    use array_utils, only: meanArithmetic
    real(r8), intent(in) :: error(:)
    l1_norm = meanArithmetic(error)
  end function l1_norm

  real(r8) function l2_norm (error)
    use array_utils, only: meanArithmetic
    real(r8), intent(in) :: error(:)
    l2_norm = sqrt(meanArithmetic(error**2))
  end function l2_norm

  real(r8) function linf_norm (error)
    real(r8), intent(in) :: error(:)
    linf_norm = maxval(error)
  end function linf_norm

end module vof_solver_type
