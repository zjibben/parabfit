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
  subroutine volume_track (volume_flux_sub, adv_dt, mesh, gmesh, Vof, fluxing_velocity, nmat, fluidRho) !result(volume_flux_sub)
    use timer_tree_type
    use hex_types, only: cell_data
    
    real(r8),         intent(in) :: adv_dt, vof(:,:), fluxing_velocity(:,:), fluidrho(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    integer,          intent(in) :: nmat
    real(r8),         intent(out) :: volume_flux_sub(nmat, 6, mesh%ncell)
    
    integer         :: ninterfaces, i, nmat_in_cell
    !real(r8)        :: int_norm(3,nmat,mesh%ncell)
    real(r8)        :: int_norm(3,nmat)
    type(cell_data) :: cell

    !call start_timer ("reconstruct/advect")   ! Start the volume track timer

    !nmat_in_cell = count(vof > 0.0_r8, dim=1)      ! count the number of materials in each cell
    ninterfaces = size(vof, dim=1) !maxval(nmat_in_cell) - 1         ! number of interfaces to process
    !int_norm = interface_normal (vof, mesh, gmesh) ! compute interface normal vectors for all the materials.
    
    ! get the volume flux in every cell
    !!$omp parallel do default(private) shared(volume_flux_sub,mesh,adv_dt,fluidRho,vof,int_norm) &
    !!$omp parallel do default(private) shared(volume_flux_sub,mesh,adv_dt,fluidRho,vof,gmesh) &
    !!$omp shared(nmat_in_cell,nmat,fluxing_velocity,ninterfaces)
    !$omp do
    do i = 1,mesh%ncell
      nmat_in_cell = count(vof(:,i) > 0.0_r8)
      int_norm = interface_normal_cell (vof, i, mesh, gmesh)
      call cell%init (mesh%x(:,mesh%cnode(:,i)), mesh%volume(i), mesh%area(mesh%cface(:,i)), &
           mesh%normal(:,mesh%cface(:,i)), mesh%cfpar(i))
      volume_flux_sub(:,:,i) = cell_volume_flux (adv_dt, cell, fluidRho(i), vof(:,i), &
           int_norm(:,:), nmat_in_cell, nmat, &
           fluxing_velocity(:,i), ninterfaces )
      ! if (i==32713) then
      !   write(*,*) 'fl',volume_flux_sub(:,:,i)
      !   write(*,*) 'dt',adv_dt
      !   write(*,*) 'n',int_norm
      !   write(*,*) 'v',fluxing_velocity(:,i)*adv_dt
      !   write(*,*) 'vof',vof(:,i)
      ! end if
    end do
    !$omp end do nowait
    !!$omp end parallel do

    !call stop_timer ("reconstruct/advect") ! Stop the volume track timer
    
  end subroutine volume_track

  ! get the volume flux for every material in the given cell
  function cell_volume_flux (adv_dt, cell, fluidRho, vof, int_norm, nmat_in_cell, nmat, face_fluxing_velocity, ninterfaces)
    use hex_types,           only: cell_data
    use locate_plane_module, only: locate_plane_hex
    use flux_volume_module,  only: flux_vol_quantity, flux_vol_vertices
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
      end if
      
      ! calculate delta advection volumes for this material at each donor face and accumulate the sum
      call material_volume_flux (cell_volume_flux(ni,:), flux_vol_sum, plane_cell, cell, &
           is_mixed_donor_cell, face_fluxing_velocity, adv_dt, vof(ni))
    end do ! interface loop

    ! get the id of the last material
    nlast = last_true_loc(vof > 0.0_r8)
    
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
    use timer_tree_type

    real(r8), intent(IN)  :: Vof(:,:)
    real(r8)              :: interface_normal(3,size(vof,1),size(vof,2))
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh

    integer :: m, n, mi, mid(2), nmat, nmat_in_cell

    call start_timer ("normal calculation")

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
    
    !$omp parallel do default(private) shared(interface_normal,vof,mesh,gmesh,nmat)
    do n = 1,size(vof, dim=2) ! loop through cells
      ! ! Calculate the volume fraction gradient for each material
      ! ! recall interface normal is in opposite direction of gradient
      ! do m = 1,nmat
      !   interface_normal(:,m,n) = -gradient_cell (Vof(m,:), n, mesh, gmesh)
      ! end do
      
      ! if there are more than 2 materials in the cell
      ! Sum the gradients in priority order.  This is equivalent to calculating the
      !  interface normal for a material composed of the first few materials
      nmat_in_cell = count(vof(:,n) > 0.0_r8)
      if (nmat_in_cell > 2) then
        do m = 2,nmat
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
        
        interface_normal(:,mid(2),n) = -interface_normal(:,mid(1),n)
      end if
      
      do m = 1,nmat
        call eliminate_noise (interface_normal(:,m,n))
        call normalize (interface_normal(:,m,n))
      end do
    end do ! cell loop
    !$omp end parallel do

    call stop_timer ("normal calculation")
    
  end function interface_normal

  function interface_normal_cell (Vof, n, mesh, gmesh) result(interface_normal)
    use timer_tree_type

    real(r8),         intent(in) :: Vof(:,:)
    integer,          intent(in) :: n
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    real(r8)                     :: interface_normal(3,size(vof,1))
    
    integer :: m, mi, mid(2), nmat, nmat_in_cell

    nmat = size(vof, dim=1)
        
    ! mmfran 07/22/11 --- begin changes
    ! bug fix for material order independent results
    ! consider special case the cell contains only two materials
    ! if cell has two materials only, their normal should be consistent
    
    ! Calculate the volume fraction gradient for each material
    ! recall interface normal is in opposite direction of gradient
    do m = 1,nmat
      interface_normal(:,m) = -gradient_cell (Vof(m,:), n, mesh, gmesh)
    end do

    ! if there are more than 2 materials in the cell
    ! Sum the gradients in priority order.  This is equivalent to calculating the
    !  interface normal for a material composed of the first few materials
    nmat_in_cell = count(vof(:,n) > 0.0_r8)
    if (nmat_in_cell > 2) then
      do m = 2,nmat
        if (Vof(m,n) > alittle) then
          Interface_Normal(:,m) = Interface_Normal(:,m) + Interface_Normal(:,1)
        else
          Interface_Normal(:,m) =                         Interface_Normal(:,1)
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

      interface_normal(:,mid(2)) = -interface_normal(:,mid(1))
    end if

    do m = 1,nmat
      call eliminate_noise (interface_normal(:,m))
      call normalize (interface_normal(:,m))
    end do
    
  end function interface_normal_cell

  subroutine normalize (n)
    real(r8), intent(inout) :: n(:)

    real(r8) :: tmp
    
    tmp = sqrt(sum(n**2))              ! calculate the denominator
    if (abs(tmp) > alittle)  n = n/tmp ! invert when it won't blow up
    call eliminate_noise (n)           ! set tiny components to zero
    if (all(n == 0.0_r8)) n = 1.0_r8   ! set all components to 1.0 in pure cells
    
  end subroutine normalize

  pure subroutine eliminate_noise (a)
    real(r8), intent(inout) :: a(3)

    if (abs(a(1)) < alittle) a(1) = 0.0_r8
    if (abs(a(2)) < alittle) a(2) = 0.0_r8
    if (abs(a(3)) < alittle) a(3) = 0.0_r8
  end subroutine eliminate_noise

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
    use timer_tree_type

    real(r8), intent(in)         :: Phi(:)
    real(r8)                     :: gradient(3,size(phi))
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh

    integer  :: f,i
    real(r8) :: Phi_f(mesh%nface)
    
    call start_timer("gradient")

    ! compute face average of phi
    ! this is calculated through the average of node values at a face
    ! the node values are calculated through the average of neighboring cell values
    ! this could be modified and moved inside the following loop, although some effort would be duplicated that way
    phi_f = face_avg_global (vertex_avg_global (Phi, mesh, gmesh), mesh)
    
    ! green-gauss method
    ! Loop over faces, accumulating the product
    ! Phi_f*Face_Normal for each area vector component

    !$omp parallel do default(private) shared(gradient,phi_f,mesh,gmesh)
    do i = 1,mesh%ncell
      ! Accumulate the dot product
      gradient(:,i) = 0.0_r8
      do f = 1,nfc
        Gradient(:,i) = Gradient(:,i) + gmesh%outnorm(:,f,i)*mesh%area(mesh%cface(f,i))*Phi_f(mesh%cface(f,i))
      end do

      Gradient(:,i) = Gradient(:,i) / mesh%volume(i) ! Normalize by cell volume
      call eliminate_noise (gradient(:,i))
    end do
    !$omp end parallel do
    
    call stop_timer("gradient")

  end function gradient

  function gradient_cell (phi, i, mesh, gmesh) !result(gradient)
    use timer_tree_type

    real(r8),         intent(in) :: phi(:)
    integer,          intent(in) :: i
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    real(r8)                     :: gradient_cell(ndim)

    integer  :: f
    
    ! green-gauss method
    ! Loop over faces, accumulating the product
    ! Phi_f*Face_Normal for each area vector component
    
    ! Accumulate the dot product
    gradient_cell = 0.0_r8
    do f = 1,nfc
      gradient_cell = gradient_cell + gmesh%outnorm(:,f,i)*mesh%area(mesh%cface(f,i))*face_avg (phi,mesh%cface(f,i), mesh, gmesh)
    end do

    gradient_cell = gradient_cell / mesh%volume(i) ! Normalize by cell volume
    call eliminate_noise (gradient_cell)

  end function gradient_cell
  
  real(r8) function face_avg (x_cell, f, mesh, gmesh)
    real(r8),         intent(in) :: x_cell(:)
    integer,          intent(in) :: f
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh

    integer  :: v

    face_avg = 0.0_r8
    do v = 1,4 ! loop through each vertex on the face
      face_avg = face_avg + vertex_avg (x_cell, mesh%fnode(v,f), mesh, gmesh)
    end do
    face_avg = face_avg / 4.0_r8
    
  end function face_avg
  
  function face_avg_global (x_vtx, mesh) result(x_face)
    real(r8),         intent(in) :: x_vtx(:)
    type(unstr_mesh), intent(in) :: mesh
    real(r8)                     :: x_face(mesh%nface)

    integer :: f,i

    !$omp parallel do default(private) shared(x_face,x_vtx,mesh)
    do f = 1,mesh%nface
      ! Interpolate vertex values to this face (see note 1)
      x_face(f) = sum(x_vtx(mesh%fnode(:,f))) / 4.0_r8
    end do
    !$omp end parallel do

  end function face_avg_global

  !=======================================================================
  ! Purpose(s):
  !
  !   Compute a vertex-averaged quantity X_vertex from a cell-centered
  !   quantity X_cell.  X_vertex is an inverse-volume-weighted average
  !   of X_cell: X_vertex = SUM(X_cell/Vol)/SUM(1/Vol).
  !
  !=======================================================================
  function vertex_avg_global (x_cell, mesh, gmesh) result(x_vertex)
    real(r8),         intent(in) :: x_cell(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    real(r8)                     :: x_vertex(mesh%nnode)

    integer  :: n

    ! calculate the inverse sum of inverse volumes
    !$omp parallel do default(private) shared(x_vertex, x_cell, mesh, gmesh)
    do n = 1,mesh%nnode
      x_vertex(n) = vertex_avg (x_cell, n, mesh, gmesh)
    end do
    !$omp end parallel do

  end function vertex_avg_global

  real(r8) function vertex_avg (x_cell, n, mesh, gmesh)
    real(r8),         intent(in) :: X_cell(:)
    integer,          intent(in) :: n
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    real(r8)                     :: X_vertex(mesh%nnode)

    integer  :: i,cellid
    real(r8) :: tmp1,tmp2

    ! calculate the inverse sum of inverse volumes
    tmp1 = 0.0_r8; tmp2 = 0.0_r8
    do i = 1,size(gmesh%vcell(:,n))
      cellid = gmesh%vcell(i,n)
      if (cellid>0) then
        tmp1 = tmp1 + X_cell(cellid) / mesh%volume(cellid)
        tmp2 = tmp2 +         1.0_r8 / mesh%volume(cellid)
      end if
    end do
    vertex_avg = tmp1/tmp2
    ! X_vertex(n) = sum( X_cell(gmesh%vcell(:,n)) / mesh%volume(gmesh%vcell(:,n)), mask=(gmesh%vcell(:,n)>0)) &
    !  /        sum(                   1.0_r8 / mesh%volume(gmesh%vcell(:,n)), mask=(gmesh%vcell(:,n)>0))

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
