module volume_track_module
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
  use kinds, only: r8
  use unstr_mesh_type, only: unstr_mesh
  implicit none
  private

  public :: volume_track

  integer, parameter :: ndim = 3
  integer, parameter :: nfc = 6 ! number of faces per cell
  integer, parameter :: nvc = 8 ! number of vertices per cell
  integer, parameter :: alittle = epsilon(1.0_r8)
  integer, parameter :: cutvof = 1.0e-8_r8
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
  function volume_track (adv_dt, mesh, Vof, Fluxing_Velocity, nmat, fluidRho) result(volume_flux_sub)
    use flux_volume_module,      only: flux_vol_quantity, flux_vol_vertices
    use timer_tree_type
    use truncate_volume_module,  only: truncate_volume, face_param, truncvol_data
        use locate_plane_module, only: interface_data, locate_plane_hex
        
    real(r8), intent(in) :: adv_dt
    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in)    :: Vof(:,:)
    real(r8), intent(in)    :: fluxing_velocity(:)
    integer,  intent(in)    :: nmat
    real(r8), intent(in)    :: fluidRho(:)
    real(r8)                :: volume_flux_sub(nmat, 6, mesh%ncell)
    
    ! Local Variables
    integer  :: interfaces, m, n, f, ff, ni, v, p, nicells, i, icell, locate_plane_niters
    integer  :: Mat(mesh%ncell)
    real(r8) :: vp, g2, vofint(mesh%ncell)
    real(r8) :: Flux_Vol_Sum(nfc,mesh%ncell)
    !real(r8) :: Xv(3,nvc,mesh%ncell)
    real(r8) :: Grad(3,nmat,mesh%ncell)
    !real(r8), allocatable :: Pack_vert(:,:,:)

    type(locate_plane_hex) :: plane_cell
    logical :: Mask(mesh%ncell)
    type(FLUX_VOL_QUANTITY) :: Flux_Vol
    type(truncvol_data) :: trunc_vol(nfc)
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! Start the volume track timer
    call start_timer ("Reconstruct/Advect")

    volume_flux_sub = 0.0_r8

    ! Count the number of materials in each cell
    Mat = 0
    do m = 1,nmat
      where (Vof(m,:) > 0.0_r8) Mat = Mat + 1
    end do

    ! Number of interfaces to process, and the maximum number of materials in a cell
    interfaces = maxval(Mat) - 1

    ! Compute interface normal vectors for all the materials.
    Grad = interface_normal (Vof, mesh)

    ! Initialize the sum of volume fluxes on each donor face.
    Flux_Vol_Sum = 0.0_r8

    ! Loop over the interfaces in priority order; if interfaces <= 0,
    ! we don't ever go into this loop.
    do ni = 1, interfaces

      ! Accumulate the volume fraction of this material and
      ! the volume fraction of materials with lower priorities.
      Vofint = sum(vof(1:ni,:), dim=1)

      ! Force 0.0 <= Vofint <= 1.0
      Vofint = min(max(Vofint, 0.0_r8), 1.0_r8)

      ! Count the number of mixed donor cells
      Mask = Vofint > (0.0_r8 + cutvof) .and. Vofint < (1.0_r8 - cutvof)
      where(fluidRho < alittle) Mask = .false.
      nicells = COUNT(Mask)

      ! Pack interface data into the Interface derived type
      icell = 1
      do i = 1,mesh%ncell
        
        ! Locate each interface plane by computing rho, the plane constant.
        if (mask(i)) then
          call plane_cell%init (grad(:,ni,i), vofint(i), mesh%volume(i), mesh%x(:,mesh%cnode(:,i)))
          call plane_cell%locate_plane (locate_plane_niters)
        end if

        ! Calculate delta advection volumes (Volume_Flux_Sub) for this material
        ! at each donor face and accumulate the sum.
        do f = 1,nfc
          ! Flux volumes
          if (Flux_Vol%Vol <= cutvof*mesh%volume(i)) then
            Flux_Vol%Fc = 0
            Flux_Vol%Vol = 0.0_r8
          else
            Flux_Vol%Fc  = f
            Flux_Vol%Vol = adv_dt*Fluxing_Velocity(mesh%cface(f,i))*mesh%area(mesh%cface(f,i))
          end if

          call flux_vol_vertices (f, mesh, Mask(i), Fluxing_Velocity, Flux_Vol, i, adv_dt)

          if (mask(i)) then
            ! Now compute the volume truncated by interface planes in each flux volumes.
            do ff = 1,nfc
              call FACE_PARAM (plane_cell, 'flux_cell', ff, trunc_vol, transpose(flux_vol%xv))
            end do
            
            ! For clean donor cells, the entire Flux volume goes to the single donor material.
            ! For mixed donor cells, the face flux is in Int_Flux%Advection_Volume.
            Vp = truncate_volume(plane_cell, trunc_vol)
          else
            if (Vof(ni,i) >= (1.0_r8-cutvof)) then
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
          if (abs(Flux_Vol%Vol) > 0.0_r8) then
            G2 = min(max(Vp - Flux_Vol_Sum(f,i), 0.0_r8), Vof(ni,i)*mesh%volume(i))
          else
            G2 = 0.0_r8
          end if

          ! Now gather advected volume information into Volume_Flux_Sub
          if (Flux_Vol%Vol > 0.0_r8) then
            Volume_Flux_Sub(ni,f,i) = G2
            Flux_Vol_Sum(f,i) = Flux_Vol_Sum(f,i) + G2
          end if

        end do ! face loop

        if (mask(i)) icell = icell+1
      end do ! cell loop

      ! ! Deallocate interface cell vertices if interface areas have been computed
      ! if (allocated(Pack_vert)) deallocate (Pack_vert)

    end do ! interface loop

    ! If I'm right, then the INTERFACE_LOOP, which looped over the number of interfaces
    ! in each cell, did so one less times than the number of materials in each cell.  And
    ! so this last loop looks at the lowest-priority material, and so-to-speak, reconciles
    ! the books.  Yes?   - Markus

    ! Make sure the volume flux sum has been initialized if there are no interfaces.
    if (interfaces == 0) Flux_Vol_Sum = 0.0_r8

    ! Compute the advection volume for the last material.
    do i = 1,mesh%ncell
      do f = 1,nfc
        ! Recalculate the total flux volume for this face.
        Flux_Vol%Vol = adv_dt*Fluxing_Velocity(mesh%cface(f,i))*mesh%area(mesh%cface(f,i))
        if (Flux_Vol%Vol <= cutvof*mesh%volume(i)) Flux_Vol%Vol = 0.0_r8
        
        ! The volume flux of the last material shouldn't be less than
        ! zero nor greater than the volume of this material in the donor cell.
        if (abs(Flux_Vol%Vol) > 0.0_r8) then
          G2 = min(max(abs(Flux_Vol%Vol - Flux_Vol_Sum(f,i)), 0.0_r8), Vof(nmat,i)*mesh%volume(i))
        else
          G2 = 0.0_r8
        end if
        
        ! Store the last material's volume flux.
        if (G2 > cutvof*mesh%volume(i)) then
          Volume_Flux_Sub(nmat,f,i) = G2
          Flux_Vol_Sum(f,i) = Flux_Vol_Sum(f,i) + G2
        end if
        
        ! For donor cells containing only one material, assign the total flux.
        if (Mat(i) == 1 .and. Flux_Vol%Vol > 0.0_r8) then
          Volume_Flux_Sub(nmat,f,i) = Flux_Vol%Vol
          Flux_Vol_Sum(f,i) = Flux_Vol%Vol
        end if
      end do ! face loop
    end do ! cell loop

    ! Stop the volume track timer
    call stop_timer ("Reconstruct/Advect")

  end function volume_track

  function interface_normal (Vof, mesh)
    !=======================================================================
    ! Purpose(s):
    !
    !   This routine calculates the interface normal for each
    !   material by computing the gradient of material volume fractions
    !
    !=======================================================================

    ! Arguments
    real(r8), intent(IN)  :: Vof(:,:)
    real(r8)              :: interface_normal(3,size(vof,1),size(vof,2))
    type(unstr_mesh), intent(in) :: mesh

    ! Local Variables
    integer :: m, n, mi, mid(2), mat(size(vof,2)), nmat
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    nmat = size(vof, dim=1)
    
    ! Calculate the volume fraction gradient for each material
    do m = 1, nmat
      interface_normal(:,m,:) = gradient (Vof(m,:), mesh)
    end do

    ! The interface normal is in the opposite sense of the gradient
    interface_normal = -interface_normal

    ! mmfran 07/22/11 --- begin changes
    ! bug fix for material order independent results
    ! consider special case the cell contains only two materials
    ! if cell has two materials only, their normal should be consistent
    Mat = 0
    do m=1,nmat
      where(Vof(m,:) > 0.0_r8) Mat = Mat + 1
    enddo

    do n = 1,size(vof, dim=2)
      ! if there are more than 2 materials in the cell
      ! Sum the gradients in priority order.  This is equivalent to calculating the
      !  interface normal for a material composed of the first few materials
      if (Mat(n) >= 3) then
        do m = 2, nmat
          if (Vof(m,n) > alittle) Interface_Normal(:,m,n) = Interface_Normal(:,m,n) + Interface_Normal(:,1,n)
          if (Vof(m,n) < alittle) Interface_Normal(:,m,n) =                           Interface_Normal(:,1,n)
        end do
      else if (mat(n) == 2) then ! if there are only two materials in the cell, ensure the normals are consistent
        mid = 0
        ! identify the material ID in the cells
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
    Tmp = SQRT(n(1,:)**2 + n(2,:)**2 + n(3,:)**2)

    ! invert when it won't blow up
    where (ABS(Tmp) > alittle)
      Tmp = 1.0_r8/Tmp
    elsewhere
      Tmp = 1.0_r8
    end where

    ! normalize
    n(1,:) = n(1,:) * Tmp
    n(2,:) = n(2,:) * Tmp
    n(3,:) = n(3,:) * Tmp

    ! set tiny components to zero
    where (ABS(n) < 1.0d-6) n = 0.0_r8

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
  function gradient (Phi, mesh)
    real(r8), intent(in)         :: Phi(:)
    real(r8)                     :: gradient(3,size(phi))
    type(unstr_mesh), intent(in) :: mesh

    integer  :: f,i
    real(r8) :: Phi_f(mesh%ncell)
    real(r8) :: Phi_vtx(mesh%nnode)

    ! Compute Phi_vtx, a vertex value of Phi
    phi_vtx = vertex_avg (Phi, mesh)

    ! green-gauss method
    ! Loop over faces, accumulating the product
    ! Phi_f*Face_Normal for each area vector component
    Gradient = 0.0_r8
    do f = 1,nfc
      ! Interpolate vertex values to this face (see note 1)
      do i = 1,mesh%ncell
        phi_f(i) = sum(phi_vtx(mesh%cnode(:,i))) / 4.0_r8
      end do

      ! Accumulate the dot product
      Gradient(1,:) = Gradient(1,:) + mesh%normal(1,mesh%cface(f,:))*Phi_f
      Gradient(2,:) = Gradient(2,:) + mesh%normal(2,mesh%cface(f,:))*Phi_f
      Gradient(3,:) = Gradient(3,:) + mesh%normal(3,mesh%cface(f,:))*Phi_f
    end do

    ! Normalize by cell volume
    Gradient(1,:) = Gradient(1,:) / mesh%volume(:)

    ! Eliminate noise
    Gradient(1,:) = MERGE(0.0_r8, Gradient(1,:), ABS(Gradient(1,:)) <= alittle)
    Gradient(2,:) = MERGE(0.0_r8, Gradient(2,:), ABS(Gradient(2,:)) <= alittle)
    Gradient(3,:) = MERGE(0.0_r8, Gradient(3,:), ABS(Gradient(3,:)) <= alittle)
  end function gradient

  !=======================================================================
  ! Purpose(s):
  !
  !   Compute a vertex-averaged quantity X_vertex from a cell-centered
  !   quantity X_cell.  X_vertex is an inverse-volume-weighted average
  !   of X_cell: X_vertex = SUM(X_cell/Vol)/SUM(1/Vol).
  !
  !=======================================================================
  function vertex_avg (x_cell, mesh) result(x_vertex)
    real(r8),         intent(in) :: X_cell(:)
    type(unstr_mesh), intent(in) :: mesh
    real(r8)                     :: X_vertex(mesh%nnode)

    integer  :: i,n, vcell(8,mesh%nnode) ! the node cell array vcell(k,j) is the cell number of local cell k of node j
    real(r8) :: rsum_rvol(mesh%nnode)

    ! get cells neighboring a vertex
    vcell = cells_neighboring_vertices (mesh)

    ! calculate the inverse sum of inverse volumes
    X_vertex = 0.0_r8; rsum_rvol = 0.0_r8
    do n = 1,mesh%nnode
      do i = 1,8
        if (vcell(i,n) /= -1) then
          X_vertex(n) = X_vertex(n) + X_cell(vcell(i,n)) / mesh%volume(vcell(i,n))
          rsum_rvol(n) = rsum_rvol(n) + 1.0_r8 / mesh%volume(vcell(i,n))
        end if
      end do
    end do
    rsum_rvol = 1.0_r8 / rsum_rvol

    ! X_vertex: Vertex-averaged X_cell
    X_vertex = X_vertex * rsum_rvol
  end function vertex_avg

  ! this function goes to each node and generates a list of cells containing that node
  ! there is probably a smarter way of doing this
  function cells_neighboring_vertices (mesh)
    type(unstr_mesh), intent(in) :: mesh
    integer :: cells_neighboring_vertices(8,mesh%nnode)

    integer :: i,n,j

    cells_neighboring_vertices = -1 

    do n = 1,mesh%nnode  ! loop through all nodes
      j = 1
      do i = 1,mesh%nnode ! at every node, loop through all cells
        if (any(mesh%cnode(:,i) == n)) then ! any of this cell's nodes are this node, add this cell to the list
          cells_neighboring_vertices(j,n) = i
          j = j + 1
        end if
      end do
    end do

  end function cells_neighboring_vertices

end module volume_track_module
