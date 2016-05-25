!=======================================================================
! Purpose:
!
!   Define the main volume-tracking procedures.
!
!   Public Interface:
!
!     * call VOLUME_TRACK (Vof, Fluxing_Velocity, Volume_Flux_Sub)
!       Control routine for a piecewise-linear volume tracking of
!       material interfaces, in which material interfaces are
!       reconstructed as planes from local volume fraction data.
!
! Author(s): Stewart J. Mosso (sjm@lanl.gov)
!            Douglas B. Kothe (dbk@lanl.gov)
!            Matthew Williams (mww@lanl.gov)
!=======================================================================

module volume_track_module

  use kinds,  only: r8
  use consts, only: cutvof, nfc
  use unstr_mesh_type
  use mesh_geom_type
  implicit none
  private

  ! material_volume_flux is only public for unit testing
  public :: volume_track, material_volume_flux

  logical,  parameter :: using_mic = .false. ! flag to select mic-optimized segments of code (currently unneeded)

contains

  !=======================================================================
  ! Purpose(s):
  !
  !   Control routine for a volume-tracking interface treatment,
  !   in which material interfaces are reconstructed as planes
  !   from local volume fraction data.  The plane interface
  !   locations are then used to find material advection volumes
  !   at all relevant cell faces.
  !
  !=======================================================================
  subroutine volume_track (volume_flux_sub, adv_dt, mesh, gmesh, Vof, fluxing_velocity, nmat, fluidRho, &
       intrec, dump_intrec) !result(volume_flux_sub)
    use timer_tree_type
    use hex_types,       only: cell_data
    use int_norm_module, only: interface_normal, interface_normal_cell
    use surface_type
    !use devices,   only: using_mic
    
    real(r8),         intent(in)  :: adv_dt, vof(:,:), fluxing_velocity(:,:), fluidrho(:)
    type(unstr_mesh), intent(in)  :: mesh
    type(mesh_geom),  intent(in)  :: gmesh
    integer,          intent(in)  :: nmat
    real(r8),         intent(out) :: volume_flux_sub(:,:,:)
    type(surface),    intent(out) :: intrec(:) ! surface reconstruction for dumping purposes
    logical,          intent(in)  :: dump_intrec
    !real(r8)                      :: volume_flux_sub(nmat, 6, mesh%ncell)

    integer               :: ninterfaces, i, ni, nmat_in_cell !(mesh%ncell)
    real(r8), allocatable :: int_norm_global(:,:,:)
    real(r8)              :: int_norm(3,nmat) !, int_norm_global(3,nmat,mesh%ncell)
    type(cell_data)       :: cell

    !!$omp barrier

    !nmat_in_cell = count(vof > 0.0_r8, dim=1)      ! count the number of materials in each cell
    !ninterfaces  = maxval(nmat_in_cell) - 1         ! number of interfaces to process
    ninterfaces = nmat-1 !size(vof, dim=1) - 1

    if (dump_intrec) then
      do ni = 1,ninterfaces
        call intrec(ni)%purge ()
      end do
    end if
    
    ! get the volume flux in every cell
    if (using_mic) then
      !$omp do
      do i = 1,mesh%ncell
        nmat_in_cell = count(vof(:,i) > 0.0_r8)
        int_norm = interface_normal_cell (vof, i, mesh, gmesh, .true.)
        call cell%init (mesh%x(:,mesh%cnode(:,i)), mesh%volume(i), mesh%area(mesh%cface(:,i)), &
             mesh%normal(:,mesh%cface(:,i)), mesh%cfpar(i))
        volume_flux_sub(:,:,i) = cell_volume_flux (adv_dt, cell, fluidRho(i), vof(:,i), &
             int_norm(:,:), nmat_in_cell, nmat, &
             fluxing_velocity(:,i), ninterfaces, dump_intrec, intrec)
      end do
      !$omp end do nowait
    else
      allocate(int_norm_global(3,nmat,mesh%ncell))
      int_norm_global = interface_normal (vof, mesh, gmesh, .true.) ! compute interface normal vectors for all materials
      call start_timer ("reconstruct/advect")

      !$omp parallel do default(private) &
      !$omp shared(volume_flux_sub, mesh, adv_dt, fluidRho, vof, int_norm_global, nmat, &
      !$omp        fluxing_velocity, ninterfaces,dump_intrec,intrec) !,nmat_in_cell,gmesh)
      do i = 1,mesh%ncell
        nmat_in_cell = count(vof(:,i) > 0.0_r8)
        call cell%init (mesh%x(:,mesh%cnode(:,i)), mesh%volume(i), mesh%area(mesh%cface(:,i)), &
             mesh%normal(:,mesh%cface(:,i)), mesh%cfpar(i))
        volume_flux_sub(:,:,i) = cell_volume_flux (adv_dt, cell, fluidRho(i), vof(:,i), &
             int_norm_global(:,:,i), nmat_in_cell, nmat, &
             fluxing_velocity(:,i), ninterfaces, dump_intrec, intrec)
      end do
      !$omp end parallel do

      deallocate(int_norm_global)
      call stop_timer ("reconstruct/advect")
    end if

  end subroutine volume_track

  ! get the volume flux for every material in the given cell
  function cell_volume_flux (adv_dt, cell, fluidRho, vof, int_norm, nmat_in_cell, nmat, &
       face_fluxing_velocity, ninterfaces, dump_intrec, intrec)
    use consts,              only: alittle
    use hex_types,           only: cell_data, hex_f, hex_e
    use locate_plane_module, only: locate_plane_hex
    use flux_volume_module,  only: flux_vol_quantity, flux_vol_vertices
    use array_utils,         only: last_true_loc
    use surface_type
    use polyhedron_type
    use logging_services
    use timer_tree_type

    real(r8),        intent(in)    :: adv_dt, int_norm(:,:), vof(:), fluidRho, face_fluxing_velocity(:)
    integer,         intent(in)    :: ninterfaces, nmat_in_cell, nmat
    type(cell_data), intent(in)    :: cell
    type(surface),   intent(inout) :: intrec(:)
    logical,         intent(in)    :: dump_intrec
    real(r8)                       :: cell_volume_flux(nmat,nfc)

    real(r8)                       :: Vofint, vp, dvol
    real(r8)                       :: flux_vol_sum(nfc)
    integer                        :: ni,f,locate_plane_niters,nlast, nint, ierr
    logical                        :: is_mixed_donor_cell
    type(locate_plane_hex)         :: plane_cell
    type(polyhedron)               :: poly
    type(flux_vol_quantity)        :: Flux_Vol

    !call start_timer ("ra_cell")

    cell_volume_flux = 0.0_r8
    flux_vol_sum = 0.0_r8
    !nint = count(vof > 0.0_r8)
    ! Here, I am not certain the conversion from pri_ptr to direct material indices worked properly.
    ! This will be clear when trying 3 or more materials. -zjibben
    
    ! Loop over the interfaces in priority order
    do ni = 1,ninterfaces
      ! check if this is a mixed material cell
      ! First accumulate the volume fraction of this material and materials with lower priorities.
      ! Force 0.0 <= Vofint <= 1.0
      Vofint = min(max(sum(vof(1:ni)), 0.0_r8), 1.0_r8)
      is_mixed_donor_cell = cutvof < Vofint .and. Vofint < (1.0_r8 - cutvof) .and. .not.fluidRho<alittle

      ! locate each interface plane by computing the plane constant
      if (is_mixed_donor_cell) then
        call plane_cell%init (int_norm(:,ni), vofint, cell%volume, cell%node)
        call plane_cell%locate_plane (locate_plane_niters)
        if (dump_intrec) then
          call poly%init (ierr, plane_cell%node, hex_f, hex_e)
          if (ierr /= 0) call LS_fatal ('failed plane reconstruction dump')
          call intrec(ni)%append (poly%intersection_verts (plane_cell%P))
        end if
      end if
      
      ! calculate delta advection volumes for this material at each donor face and accumulate the sum
      cell_volume_flux(ni,:) =  material_volume_flux (flux_vol_sum, plane_cell, cell, &
           is_mixed_donor_cell, face_fluxing_velocity, adv_dt, vof(ni))
    end do ! interface loop

    ! get the id of the last material
    nlast = last_true_loc (vof > 0.0_r8)

    !call start_timer ("last_loop")
    ! Compute the advection volume for the last material.
    do f = 1,nfc
      ! Recalculate the total flux volume for this face.
      Flux_Vol%Vol = adv_dt*face_fluxing_velocity(f)*cell%face_area(f)
      if (abs(flux_vol%vol) > 0.5_r8 * cell%volume) then
        write(*,*) adv_dt,flux_vol%vol,cell%volume,flux_vol%vol/cell%volume
        call LS_fatal('advection timestep too large')
      end if
      if (Flux_Vol%Vol <= cutvof*cell%volume) cycle
      
      ! For donor cells containing only one material, assign the total flux.
      if (nmat_in_cell==1) then
        cell_volume_flux(nlast,f) = Flux_Vol%Vol
      else
        ! The volume flux of the last material shouldn't be less than
        ! zero nor greater than the volume of this material in the donor cell.
        dvol = min(max(abs(Flux_Vol%Vol - Flux_Vol_Sum(f)), 0.0_r8), Vof(nlast)*cell%volume)

        ! Store the last material's volume flux.
        if (dvol > cutvof*cell%volume) cell_volume_flux(nlast,f) = dvol
      end if
    end do ! face loop
    !call stop_timer ("last_loop")

    !call stop_timer ("ra_cell")

  end function cell_volume_flux

  ! calculate the flux of one material in a cell
  function material_volume_flux (flux_vol_sum, plane_cell, cell, is_mixed_donor_cell, face_fluxing_velocity, adv_dt, vof)
    use hex_types,              only: cell_data
    use locate_plane_module,    only: locate_plane_hex
    use truncate_volume_module, only: truncate_volume, face_param, truncvol_data
    use flux_volume_module,     only: flux_vol_quantity, flux_vol_vertices
    use timer_tree_type

    real(r8),               intent(inout) :: flux_vol_sum(:)
    real(r8),               intent(in)    :: face_fluxing_velocity(:), adv_dt, vof
    type(locate_plane_hex), intent(in)    :: plane_cell
    type(cell_data),        intent(in)    :: cell
    logical,                intent(in)    :: is_mixed_donor_cell
    real(r8)                              :: material_volume_flux(nfc)
    
    integer                 :: f,ff, idbg
    real(r8)                :: vp, dvol
    type(truncvol_data)     :: trunc_vol(nfc)
    type(flux_vol_quantity) :: Flux_Vol

    !call start_timer ("material_flux_vol")

    material_volume_flux = 0.0_r8

    do f = 1,nfc
      ! Flux volumes
      Flux_Vol%Fc  = f
      Flux_Vol%Vol = adv_dt*face_fluxing_velocity(f)*cell%face_area(f)
      if (Flux_Vol%Vol <= cutvof*cell%volume) cycle
      
      ! calculate the vertices describing the volume being truncated through the face
      call flux_vol_vertices (f, cell, is_mixed_donor_cell, face_fluxing_velocity(f)*adv_dt, Flux_Vol)
      
      if (is_mixed_donor_cell) then
        ! Now compute the volume truncated by interface planes in each flux volumes.
        do ff = 1,nfc
          trunc_vol(ff) = face_param (plane_cell, 'flux_cell', ff, flux_vol%xv)
        end do
        
        ! For mixed donor cells, the face flux is in Int_Flux%Advection_Volume.
        Vp = truncate_volume(plane_cell, trunc_vol)
      else
        ! For clean donor cells, the entire Flux volume goes to the single donor material.
        if (vof >= (1.0_r8-cutvof)) then
          Vp = abs(flux_vol%vol)
        else
          Vp = 0.0_r8
        end if
      end if
      
      ! If Vp is close to 0 set it to 0.  If it is close
      ! to 1 set it to 1. This will avoid numerical round-off.
      if (Vp > (1.0_r8-cutvof)*abs(Flux_Vol%Vol)) Vp = abs(Flux_Vol%Vol)

      ! Make sure that the current material-integrated advection
      ! volume hasn't decreased from its previous value.  (This
      ! can happen if the interface significantly changed its
      ! orientation and now crosses previous interfaces.)  Also
      ! limit Volume_Flux_Sub to take no more than the material
      ! occupied volume in the donor cell.
      dvol = min(max(Vp - Flux_Vol_Sum(f), 0.0_r8), Vof*cell%volume)
      
      ! Now gather advected volume information into Volume_Flux_Sub
      material_volume_flux(f) = dvol
      Flux_Vol_Sum(f) = Flux_Vol_Sum(f) + dvol
    end do ! face loop

    !call stop_timer ("material_flux_vol")

  end function material_volume_flux
  
end module volume_track_module
