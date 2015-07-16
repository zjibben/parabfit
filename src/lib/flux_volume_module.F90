module flux_volume_module
  !=======================================================================
  ! Purpose(s):
  !
  !   Define the variables and procedures associated with 
  !   computing the advection flux volume.
  !
  ! Public Interface(s):
  !
  !   * call FLUX_VOL_VERTICES (Fluxing_Velocity, Flux_Vol)
  !
  !     Compute the vertices that describe the flux volume.
  !
  ! Contains: FLUX_VOL_VERTICES
  !           SCREWED_VOLUME
  !
  ! Author(s): Douglas B. Kothe (LANL Group T-3, dbk@lanl.gov)
  !            S. Jay Mosso (LANL Group X-HM, sjm@lanl.gov)
  !
  !=======================================================================
  use kinds, only: r8
  use logging_services
  implicit none
  private

  public :: flux_vol_quantity, flux_vol_vertices
  
  ! Maximum allowable iterations for the flux volume vertex location
  integer,  parameter :: flux_vol_iter_max = 10
  integer,  parameter :: ndim = 3
  integer,  parameter :: nvc = 8
  integer,  parameter :: nvf = 4
  integer,  parameter :: nfc = 6
  real(r8), parameter :: alittle = epsilon(1.0_r8)
  real(r8), parameter :: cutvof = 1.0e-8_r8

  ! Define flux_vol_quantity structure
  type flux_vol_quantity
    integer  :: Fc           ! face number through which we this cell is fluxing
    real(r8) :: Vol          ! Volume of the flux
    real(r8) :: Xv(ndim,nvc) ! Vertices of the flux volume
  end type flux_vol_quantity

contains
  
  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>
  !=======================================================================
  ! Purpose(s):
  !
  !   Given the value of Flux_Vol (the volume of material that moves
  !   through the current advection cell face), what are the vertices
  !   that describe this volume.  Four of these vertices will be the ones
  !   that describe the advection cell face.  The other four vertices will
  !   lie approximately "DIST" away from the advection cell face.  These
  !   are only approximately "DIST" away because the cell cross-sectional
  !   area may increase or decrease as one moves away from the advection
  !   cell face.  The value used is varied from "DIST" such that the 
  !   vertices describe a hexagonal volume that matches the value of
  !   Flux_Vol.
  !
  !=======================================================================
  subroutine flux_vol_vertices (face, cell, is_mixed_donor_cell, dist, Flux_Vol)
    use hex_types, only: cell_data
    use cell_geometry, only: eval_hex_volumes

    integer,                 intent(in)    :: face
    type(cell_data),         intent(in)    :: cell
    logical,                 intent(in)    :: is_mixed_donor_cell
    real(r8),                intent(in)    :: dist
    type(flux_vol_quantity), intent(inout) :: flux_vol

    integer :: n, v, ia, ib, e, iter, &
         ! Edge_ends(2,nvf,nfc) = reshape( [ &
         ! 3,2,  4,1,  7,6,  8,5,   & ! face 1 edges
         ! 1,4,  2,3,  5,8,  6,7,   & ! face 2 edges
         ! 1,2,  4,3,  5,6,  8,7,   & ! face 3 edges
         ! 2,1,  3,4,  6,5,  7,8,   & ! face 4 edges
         ! 1,5,  2,6,  3,7,  4,8,   & ! face 5 edges
         ! 5,1,  6,2,  7,3,  8,4 ], & ! face 6 edges
         ! Edge_ends(2,nvf,nfc) = reshape( [ &
         ! 1,4,  2,3,  5,8,  6,7,   & ! face 1 edges
         ! 3,2,  4,1,  7,6,  8,5,   & ! face 2 edges
         ! 3,4,  2,1,  7,8,  6,5,   & ! face 3 edges
         ! 4,3,  1,2,  8,7,  5,6,   & ! face 4 edges
         ! 3,7,  4,8,  1,5,  2,6,   & ! face 5 edges
         ! 7,3,  8,4,  5,1,  6,2 ], & ! face 6 edges
         ! [2,nvf,nfc] )
         Edge_ends(2,nvf,nfc) = reshape( [ &
         3,2,  4,1,  7,6,  8,5,   & ! face 2 edges
         1,4,  2,3,  5,8,  6,7,   & ! face 1 edges
         4,3,  1,2,  8,7,  5,6,   & ! face 4 edges
         3,4,  2,1,  7,8,  6,5,   & ! face 3 edges
         3,7,  4,8,  1,5,  2,6,   & ! face 5 edges
         7,3,  8,4,  5,1,  6,2 ], & ! face 6 edges
         [2,nvf,nfc] )
    real(r8)       :: Percnt(nvf), tmp(8)
    real(r8)       :: Uedge(ndim,nvf)
    real(r8)       :: Volume, Mult, ndotuedge
    logical        :: converged
    character(256) :: errmsg

    ! first gather the vertex coordinates and store them in flux_vol%Xv
    Flux_Vol%Xv = cell%node
    
    ! need only go through this exercise for mixed donor cells and for outgoing fluxes
    if (.not.is_mixed_donor_cell .or. dist*cell%face_area(face) <= cutvof*cell%volume) return

    ! compute the edge unit vectors
    do e = 1, nvf
      ia  = Edge_ends(1,e,face) ! front node
      ib  = Edge_ends(2,e,face) ! back node
      Uedge(:,e) = Flux_Vol%Xv(:,ib) - Flux_Vol%Xv(:,ia)
      
      ndotuedge = sum(cell%face_normal(:,face) * Uedge(:,e))
      Percnt(e) = -Dist/(ndotuedge+alittle)
    end do

    if (any(Percnt < 0.0_r8) .or. any(Percnt > 1.0_r8)) then
      write(errmsg,'(a,3es13.6)') 'FLUX_VOL_VERTICES: invalid flux volume or inverted element; cell centroid =', &
           sum( cell%node, dim=2 ) / real(nvc,r8)
      call LS_fatal (errmsg)
    end if

    ! iterate to find the 4 other vertices that define the back end of the flux volume
    mult = 1.0_r8
    do iter = 1,flux_vol_iter_max
      ! loop over edges to adjust the vertices
      do e = 1,nvf
        ia = Edge_ends(1,e,face)
        ib = Edge_ends(2,e,face)
        Flux_Vol%Xv(:,ib) = Flux_Vol%Xv(:,ia) + Mult*Percnt(e)*Uedge(:,e)
      end do
      
      ! compute the flux volume bounded by the computed vertices
      call eval_hex_volumes (flux_vol%xv, volume, tmp)
      ! volume = screwed_volume (Flux_Vol%Xv)

      ! compare this volume with the actual Flux_Vol
      if (abs(flux_vol%vol - volume) < cutvof*cell%volume) then
        return
      else ! increment multiplier for next iteration
        Mult = Mult * Flux_Vol%Vol/volume
      end if
    end do

    ! print out a warning message if we iterated up to the maximum
    write(errmsg,'(a,i0,a,es12.5,a,i0)') &
         'Flux volume vertex iteration did not converge in ', flux_vol_iter_max,&
         '. Maximum flux volume difference is ', abs(Flux_Vol%Vol-Volume) !,&
    !' in cell '!, cell_index
    call LS_warn (errmsg)
  end subroutine flux_vol_vertices

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>
  
  ! !=======================================================================
  ! ! Purpose(s):
  ! !
  ! !   Compute a hexahedral cell volume, according to the prescription
  ! !   given by J. Dukowicz, JCP 74: 493-496 (1988).
  ! !
  ! !=======================================================================
  ! real(r8) function screwed_volume (node)
  !   use hex_types, only: face_node
    
  !   real(r8), intent(in) :: node(:,:)
    
  !   real(r8) :: X1(ndim), X2(ndim), X3(ndim)
  !   integer  :: f, v1, v2, v3, v4, v5, v6, n, v(6)
  !   character(256) :: errmsg
    
  !   Screwed_Volume = 0.0_r8

  !   ! Loop over the six faces of the flux volume
  !   do f = 1,6
  !     ! select case (f)
  !     ! case (1) ! Side 1 (vertices 4-8-7-3)
  !     !   v1 = 8; v2 = 4; v3 = 7; v4 = 8; v5 = 3; v6 = 4
  !     ! case (2) ! Side 2 (vertices 1-2-6-5)
  !     !   v1 = 6; v2 = 2; v3 = 5; v4 = 6; v5 = 1; v6 = 2
  !     ! case (3) ! Side 3 (vertices 4-1-5-8)
  !     !   v1 = 5; v2 = 1; v3 = 8; v4 = 5; v5 = 4; v6 = 1
  !     ! case (4) ! Side 4 (vertices 3-7-6-2)
  !     !   v1 = 7; v2 = 3; v3 = 6; v4 = 7; v5 = 2; v6 = 3
  !     ! case (5) ! Side 5 (vertices 4-3-2-1)
  !     !   v1 = 3; v2 = 4; v3 = 2; v4 = 3; v5 = 1; v6 = 4
  !     ! case (6) ! Side 6 (vertices 8-5-6-7)
  !     !   v1 = 6; v2 = 5; v3 = 7; v4 = 6; v5 = 8; v6 = 5
  !     ! end select
      
  !     select case (f)
  !     case (1)
  !       v = face_node([2,1,3,2,4,1],f)
  !     case (2,3,6)
  !       v = face_node([4,3,1,4,2,3],f)
  !     case (4,5)
  !       v = face_node([1,4,2,1,3,4],f)
  !     end select
      

  !     X1(:) = node(:,v(1)) + node(:,v(2))
  !     X2(:) = node(:,v(3)) + node(:,v(4))
  !     X3(:) = node(:,v(5)) + node(:,v(6))

  !     Screwed_Volume = Screwed_Volume + X1(1)*(X2(2)*X3(3) - X3(2)*X2(3)) &
  !          +                            X1(2)*(X3(1)*X2(3) - X2(1)*X3(3)) &
  !          +                            X1(3)*(X2(1)*X3(2) - X3(1)*X2(2))
  !   end do

  !   Screwed_Volume = Screwed_Volume / 12.0_r8

  !   ! Punt if Volume < 0
  !   if (screwed_volume < 0.0_r8) then
  !     write(errmsg,'(a,i0,a)') 'SCREWED_VOLUME: cell '!, cell_index, ' contains a negative flux volume'
  !     call LS_fatal (errmsg)
  !   end if
  ! end function screwed_volume

end module flux_volume_module
