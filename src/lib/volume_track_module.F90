!=======================================================================
! Purpose:
!
!   Define the main volume-tracking procedures.
!
!   Public Interface:
!
!     * call VOLUME_TRACK (Vof, Fluxing_Velocity, Volume_Flux_Sub)
!
!       Control routine for a piecewise-linear volume tracking of
!       material interfaces, in which material interfaces are
!       reconstructed as planes from local volume fraction data.
!
! Contains: VOLUME_TRACK
!           INT_NORMAL
!
! Author(s): Stewart J. Mosso (sjm@lanl.gov)
!            Douglas B. Kothe (dbk@lanl.gov)
!            Matthew Williams (mww@lanl.gov)
!=======================================================================

module volume_track_module
  use kinds, only: r8
  use unstr_mesh_type
  use mesh_geom_type
  implicit none
  private

  public :: volume_track, material_flux_unit_test

  integer, parameter :: ndim = 3
  integer, parameter :: nfc = 6 ! number of faces per cell
  integer, parameter :: nvc = 8 ! number of vertices per cell
  integer, parameter :: alittle = epsilon(1.0_r8)
  integer, parameter :: cutvof = 1.0e-8_r8

contains
  
  logical function material_flux_unit_test (plane_cell, face_fluxing_velocity, fluxex)
    use kinds, only: r8
    use hex_types, only: cell_data
    use locate_plane_module
    
    type(locate_plane_hex), intent(in) :: plane_cell
    real(r8),               intent(in) :: fluxex(:), face_fluxing_velocity(:)
    
    type(cell_data) :: cell
    real(r8)        :: mat_vol_flux(6), flux_vol_sum(6)

    call cell%init (plane_cell%node)
    
    mat_vol_flux = 0.0_r8
    flux_vol_sum = 0.0_r8
    
    call material_volume_flux (mat_vol_flux, flux_vol_sum, plane_cell, cell, .true., face_fluxing_velocity, 1.0_r8, 0.5_r8)
    
    write(*,*) 'mat_vol_flux', mat_vol_flux,'exact: ',fluxex

    material_flux_unit_test = all(mat_vol_flux==fluxex)
    
  end function material_flux_unit_test
  
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
  function volume_track (adv_dt, mesh, gmesh, Vof, fluxing_velocity, nmat, fluidRho) result(volume_flux_sub)
    use timer_tree_type
    use hex_types, only: cell_data
    
    real(r8),         intent(in) :: adv_dt, vof(:,:), fluxing_velocity(:,:), fluidrho(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    integer,          intent(in) :: nmat
    real(r8)                     :: volume_flux_sub(nmat, 6, mesh%ncell)
    
    integer         :: ninterfaces, i, nmat_in_cell(mesh%ncell)
    real(r8)        :: int_norm(3,nmat,mesh%ncell)
    type(cell_data) :: cell
    logical         :: verbose
    
    call start_timer ("Reconstruct/Advect")   ! Start the volume track timer

    nmat_in_cell = count(vof > 0.0_r8, dim=1)      ! count the number of materials in each cell
    ninterfaces = maxval(Nmat_In_Cell) - 1         ! number of interfaces to process
    int_norm = interface_normal (Vof, mesh, gmesh) ! compute interface normal vectors for all the materials.
    
    ! write(*,'(a,3es13.4)') 'norm1:     ',int_norm(:,1,8495)
    ! write(*,'(a,3es13.4)') 'norm2:     ',int_norm(:,2,8495)
    ! write(*,'(a,6es13.4)') 'flux_vel: ',fluxing_velocity(:,8495)
    
    ! get the volume flux in every cell
    !$omp parallel do default(private) shared(volume_flux_sub,mesh,adv_dt,fluidRho,vof,int_norm,nmat_in_cell,nmat,fluxing_velocity,ninterfaces)
    do i = 1,mesh%ncell
      verbose = .false.
      call cell%init (mesh%x(:,mesh%cnode(:,i)), mesh%volume(i), mesh%area(mesh%cface(:,i)), mesh%normal(:,mesh%cface(:,i)), mesh%cfpar(i))
      volume_flux_sub(:,:,i) = cell_volume_flux (adv_dt, cell, fluidRho(i), vof(:,i), int_norm(:,:,i), nmat_in_cell(i), nmat, &
            fluxing_velocity(:,i), ninterfaces )
    end do ! cell loop
    !$omp end parallel do
        
    call stop_timer ("Reconstruct/Advect") ! Stop the volume track timer

  end function volume_track

  ! get the volume flux for every material in the given cell
  function cell_volume_flux (adv_dt, cell, fluidRho, vof, int_norm, nmat_in_cell, nmat, face_fluxing_velocity, ninterfaces)
    use hex_types,               only: cell_data
    use locate_plane_module,     only: locate_plane_hex
    use flux_volume_module,      only: flux_vol_quantity, flux_vol_vertices
    use logging_services
    
    real(r8), intent(in)        :: adv_dt, int_norm(:,:), vof(:), fluidRho, face_fluxing_velocity(6)
    integer,  intent(in)        :: ninterfaces, nmat_in_cell, nmat
    type(cell_data), intent(in) :: cell
    real(r8)                    :: cell_volume_flux(nmat,6)

    real(r8)                :: Vofint, vp, dvol
    real(r8)                :: flux_vol_sum(nfc)
    integer                 :: ni,f,locate_plane_niters,nlast
    logical                 :: is_mixed_donor_cell
    type(locate_plane_hex)  :: plane_cell
    type(flux_vol_quantity) :: Flux_Vol
    
    cell_volume_flux = 0.0_r8
    flux_vol_sum = 0.0_r8

    ! Here, I am not certain the conversion from pri_ptr to direct material indices worked properly. -zjibben
    
    ! Loop over the interfaces in priority order
    do ni = 1,ninterfaces
      ! check if this is a mixed material cell
      ! First accumulate the volume fraction of this material and materials with lower priorities.
      ! Force 0.0 <= Vofint <= 1.0
      Vofint = min(max(sum(vof(1:ni)), 0.0_r8), 1.0_r8)
      is_mixed_donor_cell = (0.0_r8 + cutvof) < Vofint .and. Vofint < (1.0_r8 - cutvof) .and. .not.fluidRho<alittle

      ! locate each interface plane by computing the plane constant
      if (is_mixed_donor_cell) then
        call plane_cell%init (int_norm(:,ni), vofint, cell%volume, cell%node)
        call plane_cell%locate_plane (locate_plane_niters)
        !plane_cell%normal = -plane_cell%normal ! flux routines want the negative of the normal for some reason -zjibben
      end if
      
      ! calculate delta advection volumes for this material at each donor face and accumulate the sum
      call material_volume_flux (cell_volume_flux(ni,:), flux_vol_sum, plane_cell, cell, &
           is_mixed_donor_cell, face_fluxing_velocity, adv_dt, vof(ni))
    end do ! interface loop

    ! get the id of the last material
    nlast = last_true_loc(vof > 0.0_r8)
    if (nlast < 1) write(*,*) 'ERROR'
    
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

  end function cell_volume_flux

  ! calculate the flux of one material in a cell
  subroutine material_volume_flux (mat_vol_flux, flux_vol_sum, plane_cell, cell, is_mixed_donor_cell, face_fluxing_velocity, adv_dt, vof)
    use hex_types,              only: cell_data
    use locate_plane_module,    only: locate_plane_hex
    use truncate_volume_module, only: truncate_volume, face_param, truncvol_data
    use flux_volume_module,     only: flux_vol_quantity, flux_vol_vertices
    
    real(r8),               intent(out)   :: mat_vol_flux(:)
    real(r8),               intent(inout) :: flux_vol_sum(:)
    real(r8),               intent(in)    :: face_fluxing_velocity(:), adv_dt, vof
    type(locate_plane_hex), intent(in)    :: plane_cell
    type(cell_data),        intent(in)    :: cell
    logical,                intent(in)    :: is_mixed_donor_cell
    
    integer                 :: f,ff, idbg
    real(r8)                :: vp, dvol
    type(truncvol_data)     :: trunc_vol(nfc)
    type(flux_vol_quantity) :: Flux_Vol

    mat_vol_flux = 0.0_r8

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
      mat_vol_flux(f) = dvol
      Flux_Vol_Sum(f) = Flux_Vol_Sum(f) + dvol
    end do ! face loop
    
  end subroutine material_volume_flux
  
  !=======================================================================
  ! Purpose(s):
  !
  !   This routine calculates the interface normal for each
  !   material by computing the gradient of material volume fractions
  !
  !=======================================================================
  function interface_normal (Vof, mesh, gmesh)
    real(r8), intent(IN)  :: Vof(:,:)
    real(r8)              :: interface_normal(3,size(vof,1),size(vof,2))
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh

    integer :: m, n, mi, mid(2), nmat, nmat_in_cell

    nmat = size(vof, dim=1)
    
    ! Calculate the volume fraction gradient for each material
    ! recall interface normal is in opposite direction of gradient
    do m = 1, nmat
      interface_normal(:,m,:) = -gradient (Vof(m,:), mesh, gmesh)
    end do
    
    ! mmfran 07/22/11 --- begin changes
    ! bug fix for material order independent results
    ! consider special case the cell contains only two materials
    ! if cell has two materials only, their normal should be consistent

    !$omp parallel do default(private) shared(interface_normal,vof)
    do n = 1,size(vof, dim=2) ! loop through cells
      nmat_in_cell = count(vof(:,n) > 0.0_r8)
      
      ! if there are more than 2 materials in the cell
      ! Sum the gradients in priority order.  This is equivalent to calculating the
      !  interface normal for a material composed of the first few materials
      if (nmat_in_cell > 2) then
        do m = 2, nmat
          if (Vof(m,n) > alittle) then
            Interface_Normal(:,m,n) = Interface_Normal(:,m,n) + Interface_Normal(:,1,n)
          else
            Interface_Normal(:,m,n) =                           Interface_Normal(:,1,n)
          end if
        end do
      else if (nmat_in_cell == 2) then ! if there are only two materials in the cell, ensure the normals are consistent
        mid = 0
        ! identify the material IDs in the cells
        mi=1; m=1
        do while (mi<=2)
          if (Vof(m,n) > 0.0_r8) then
            mid(mi) = m
            mi = mi+1
          endif
          m = m+1
        end do
        
        Interface_Normal(:,mid(2),n) = - Interface_Normal(:,mid(1),n)
      end if
    end do ! cell loop
    !$omp end parallel do
        
    ! Eliminate noise from the gradient
    where (abs(Interface_Normal) < alittle) Interface_Normal = 0.0_r8
        
    ! normalize the gradient
    do m = 1,nmat
      call normalize (interface_normal(:,m,:))
    end do
    
  end function interface_normal

  subroutine normalize (n)
    real(r8), intent(inout) :: n(:,:)

    real(r8) :: tmp(size(n, dim=2))

    ! calculate the denominator
    tmp = sqrt(n(1,:)**2 + n(2,:)**2 + n(3,:)**2)

    ! invert when it won't blow up
    where (abs(tmp) > alittle)
      tmp = 1.0_r8/Tmp
    elsewhere
      tmp = 1.0_r8
    end where

    ! normalize
    n(1,:) = n(1,:) * tmp
    n(2,:) = n(2,:) * tmp
    n(3,:) = n(3,:) * tmp

    ! set tiny components to zero
    where (abs(n) < 1.0d-6) n = 0.0_r8

    ! set all components to 1.0 in pure cells
    where (n(1,:) == 0.0_r8 .and. n(2,:) == 0.0_r8 .and. n(3,:) == 0.0_r8)
      n(1,:) = 1.0_r8
      n(2,:) = 1.0_r8
      n(3,:) = 1.0_r8
    end where
    
  end subroutine normalize

  !=======================================================================
  ! Purpose(s):
  ! 
  !   Compute the cell-centered gradient (dPhi_dx,dPhi_dy,dPhi_dz) of a
  !   cell-centered scalar quantity Phi with the green-gauss method:
  !       With a discrete approximation to Gauss's theorem, whereby
  !       the integral of (dPhi_dx,dPhi_dy,dPhi_dz) over the cell
  !       volume is converted to an integral of Phi over the cell
  !       surface area vector.  The area integral is approximated
  !       as a discrete sum over cell faces.  This control volume
  !       formulation is discretely conservative, i.e., adjacent
  !       face contributions will telescope, leaving only boundary
  !       contributions.
  !
  !  note 1: Truchas does this with linear prop. It could also possibly
  !          be done locally just through the 2 cells sharing the face,
  !          using some flux function.
  !
  !=======================================================================
  function gradient (Phi, mesh, gmesh)
    real(r8), intent(in)         :: Phi(:)
    real(r8)                     :: gradient(3,size(phi))
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh

    integer  :: f,i
    real(r8) :: Phi_f(mesh%nface)

    ! compute face average of phi
    ! this is calculated through the average of node values at a face
    ! the node values are calculated through the average of neighboring cell values
    ! this could be modified and moved inside the following loop, although some effort would be duplicated that way
    phi_f = face_avg (vertex_avg (Phi, mesh, gmesh), mesh)

    ! green-gauss method
    ! Loop over faces, accumulating the product
    ! Phi_f*Face_Normal for each area vector component
    gradient = 0.0_r8
    
    !$omp parallel do default(private) shared(gradient,phi_f,mesh,gmesh)
    do i = 1,mesh%ncell
      ! Accumulate the dot product
      do f = 1,nfc
        Gradient(:,i) = Gradient(:,i) + gmesh%outnorm(:,f,i)*mesh%area(mesh%cface(f,i))*Phi_f(mesh%cface(f,i))
      end do

      ! Normalize by cell volume
      Gradient(:,i) = Gradient(:,i) / mesh%volume(i)
      
      ! Eliminate noise
      where (abs(gradient(:,i)) <= alittle) gradient(:,i) = 0.0_r8
    end do
    !$omp end parallel do
        
  end function gradient

  function face_avg (x_vtx, mesh) result(x_face)
    real(r8),         intent(in) :: x_vtx(:)
    type(unstr_mesh), intent(in) :: mesh
    real(r8)                     :: x_face(mesh%nface)

    integer :: f

    !$omp parallel do default(private) shared(x_face,x_vtx,mesh)
    do f = 1,mesh%nface
      ! Interpolate vertex values to this face (see note 1)
      x_face(f) = sum(x_vtx(mesh%fnode(:,f))) / 4.0_r8
    end do
    !$omp end parallel do

  end function face_avg
  
  !=======================================================================
  ! Purpose(s):
  !
  !   Compute a vertex-averaged quantity X_vertex from a cell-centered
  !   quantity X_cell.  X_vertex is an inverse-volume-weighted average
  !   of X_cell: X_vertex = SUM(X_cell/Vol)/SUM(1/Vol).
  !
  !=======================================================================
  function vertex_avg (x_cell, mesh, gmesh) result(x_vertex)
    real(r8),         intent(in) :: X_cell(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    real(r8)                     :: X_vertex(mesh%nnode)

    integer  :: i,n
    real(r8) :: sum_rvol(mesh%nnode)

    ! calculate the inverse sum of inverse volumes
    X_vertex = 0.0_r8; sum_rvol = 0.0_r8

    !$omp parallel do default(private) shared(X_vertex, sum_rvol, X_cell, gmesh, mesh)
    do n = 1,mesh%nnode
      do i = 1,8
        if (gmesh%vcell(i,n) /= -1) then
          X_vertex(n) = X_vertex(n) + X_cell(gmesh%vcell(i,n)) / mesh%volume(gmesh%vcell(i,n))
          sum_rvol(n) = sum_rvol(n) + 1.0_r8 / mesh%volume(gmesh%vcell(i,n))
        end if
      end do
    end do
    !$omp end parallel do

    ! X_vertex: Vertex-averaged X_cell
    X_vertex = X_vertex / sum_rvol
  end function vertex_avg

  ! return the index of the last true element of a logical array
  ! this function is used in a few different modules -- it may be a good idea to bring it into one module
  function last_true_loc (mask)
    logical, intent(in) :: mask(:)
    integer :: last_true_loc
    
    integer :: i

    do i = size(mask),1,-1
      if (mask(i)) then
        last_true_loc = i
        return
      end if
    end do

    last_true_loc = 0
    
  end function last_true_loc
  
end module volume_track_module
