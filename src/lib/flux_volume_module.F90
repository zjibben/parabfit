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

  public :: FLUX_VOL_QUANTITY, FLUX_VOL_VERTICES

  ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

  ! Maximum allowable iterations for the flux volume vertex location
  integer, public, save :: flux_vol_iter_max

  integer,  parameter :: ndim = 3
  integer,  parameter :: nvc = 8
  integer,  parameter :: nvf = 4
  integer,  parameter :: nfc = 6
  real(r8), parameter :: alittle = epsilon(1.0_r8)
  real(r8), parameter :: cutvof = 1.0e-8_r8
  
  ! Define FLUX structure
  type FLUX_VOL_QUANTITY
  
     ! Face number through which we this cell is fluxing
     integer :: Fc

     ! Volume of the flux
     real(r8) :: Vol

     ! Vertices of the flux volume
     real(r8), dimension(nvc,ndim) :: Xv
     
  end type FLUX_VOL_QUANTITY

contains

  ! <><><><><><><><><><><><> PUBLIC ROUTINES <><><><><><><><><><><><><><><>
  
  subroutine flux_vol_vertices (face, mesh, Mask, Fluxing_Velocity, Flux_Vol, cell_index, adv_dt)
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
    use unstr_mesh_type, only: unstr_mesh

    ! Arguments
    integer,  intent(IN) :: face
    type(unstr_mesh), intent(in) :: mesh
    logical,  intent(IN) :: Mask
    real(r8), intent(IN) :: Fluxing_Velocity(:)
    type(FLUX_VOL_QUANTITY), intent(INOUT) :: Flux_Vol
    integer, intent(in) :: cell_index
    real(r8), intent(in) :: adv_dt

    ! Local Variables
    integer :: n, v, ia, ib, e, iter
    real(r8) :: Percnt(nvf)
    real(r8) :: Vtx(nvc)
    real(r8) :: Uedge(nvf,ndim)
    real(r8) :: Volume, Mult, Tmp, Dist
    integer :: &
         Edge_ends(2,nvf,nfc) = reshape( [ &
         3,2,  4,1,  7,6,  8,5,   & ! face one edges
         1,4,  2,3,  5,8,  6,7,   & ! face two edges
         1,2,  4,3,  5,6,  8,7,   & ! face three edges
         2,1,  3,4,  6,5,  7,8,   & ! face four edges
         1,5,  2,6,  3,7,  4,8,   & ! face five edges
         5,1,  6,2,  7,3,  8,4 ], & ! face six edges
         [2,nvf,nfc] )
    logical        :: converged
    character(256) :: errmsg

    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    ! First gather the vertex coordinates and store them in Vtx
    do v = 1,nvc
       Flux_Vol%Xv(v,:) = mesh%x(:,mesh%cnode(v,cell_index))
    end do
        
    ! Need only go through this exercise for Mask = .true., and for outgoing fluxes
    if (.not.Mask .or. &
         Fluxing_Velocity(mesh%cface(face,cell_index))*mesh%area(mesh%cface(face,cell_index))*adv_dt &
         <= cutvof*mesh%Volume(cell_index)) return

    ! An initial guess for the depth of the flux volume is 'Dist'
    Dist = Fluxing_Velocity(mesh%cface(face,cell_index)) * adv_dt

    ! Compute the edge unit vectors
    do e = 1, nvf
       ia  = Edge_ends(1,e,face)
       ib  = Edge_ends(2,e,face)

       Tmp = 0.0_r8
       do n = 1,ndim
          Uedge(e,n) = Flux_Vol%Xv(ib,n) - Flux_Vol%Xv(ia,n)
          Tmp = Tmp + mesh%normal(n,mesh%cface(face,cell_index))/mesh%area(mesh%cface(face,cell_index)) * Uedge(e,n)
       end do
       Percnt(e) = -Dist/(Tmp+alittle)
    end do

    do e = 1, nvf
       if (Percnt(e) < 0.0_r8 .or. Percnt(e) > 1.0_r8) then
          write(errmsg,'(a,3es13.6)') 'FLUX_VOL_VERTICES: invalid flux volume or inverted element; cell centroid =', &
               sum( mesh%x(:,mesh%cnode(:,cell_index)), dim=2 ) / real(nvc,r8)
          call LS_fatal (errmsg)
       end if
    end do

    ! Initialize the scale factor to unity
    Mult = 1.0_r8

    converged = .false.

    ! Iterate to find the 4 other vertices that define the back end of the flux volume
    do iter = 1,flux_vol_iter_max

       ! Loop over edges to adjust the vertices
       do e = 1,nvf

          ia = Edge_ends(1,e,face)
          ib = Edge_ends(2,e,face)

          do n = 1,ndim
             Flux_Vol%Xv(ib,n) = Flux_Vol%Xv(ia,n) + Mult*Percnt(e)*Uedge(e,n)
          end do

       end do

       ! Compute the flux volume bounded by the computed vertices
       volume = screwed_volume (Flux_Vol)
       
       ! Punt if Volume < 0
       if (Volume < 0.0_r8) then
          write(errmsg,'(a,i0,a)') 'SCREWED_VOLUME: cell ', cell_index, ' contains a negative flux volume'
          call LS_fatal (errmsg)
       end if
       
       ! Now compare this volume with the actual Flux_Vol
       if (ABS(Flux_Vol%Vol - Volume) < cutvof*mesh%Volume(cell_index)) then
          converged = .true.
          exit
       else
          Mult = Mult * Flux_Vol%Vol/Volume
       endif

    end do

    ! Print out a warning message if we iterated up to the maximum
    if (.not.converged) then
       write(errmsg,'(a,i0,a,es12.5,a,i0)') &
            'Flux volume vertex iteration did not converge in ', &
            flux_vol_iter_max, '. Maximum flux volume difference is ', &
            ABS(Flux_Vol%Vol-Volume), ' in cell ', cell_index
       call TLS_warn (errmsg)
    end if
    

  end subroutine flux_vol_vertices

  ! <><><><><><><><><><><><> PRIVATE ROUTINES <><><><><><><><><><><><><><><>

  real(r8) function screwed_volume (Flux_Vol)
    !=======================================================================
    ! Purpose(s):
    !
    !   Compute a hexahedral cell volume, according to the prescription
    !   given by J. Dukowicz, JCP 74: 493-496 (1988).
    !
    !=======================================================================

    ! Arguments
    type(FLUX_VOL_QUANTITY), intent(IN) :: Flux_Vol

    ! Local Variables
    real(r8), dimension(ndim) :: X1, X2, X3
    integer :: f, v1, v2, v3, v4, v5, v6, n
    ! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

    Screwed_Volume = 0.0_r8

    ! Loop over the six faces of the flux volume
    do f = 1,6

       select case (f)
          case (1)
             ! Side 1 (vertices 4-8-7-3)
             v1 = 8; v2 = 4; v3 = 7; v4 = 8; v5 = 3; v6 = 4
          case (2)
             ! Side 2 (vertices 1-2-6-5)
             v1 = 6; v2 = 2; v3 = 5; v4 = 6; v5 = 1; v6 = 2
          case (3)
             ! Side 3 (vertices 4-1-5-8)
             v1 = 5; v2 = 1; v3 = 8; v4 = 5; v5 = 4; v6 = 1
          case (4)
             ! Side 4 (vertices 3-7-6-2)
             v1 = 7; v2 = 3; v3 = 6; v4 = 7; v5 = 2; v6 = 3
          case (5)
             ! Side 5 (vertices 4-3-2-1)
             v1 = 3; v2 = 4; v3 = 2; v4 = 3; v5 = 1; v6 = 4
          case (6)
             ! Side 6 (vertices 8-5-6-7)
             v1 = 6; v2 = 5; v3 = 7; v4 = 6; v5 = 8; v6 = 5
       end select

       do n = 1,ndim
          X1(n) = Flux_Vol%Xv(v1,n) + Flux_Vol%Xv(v2,n)
          X2(n) = Flux_Vol%Xv(v3,n) + Flux_Vol%Xv(v4,n)
          X3(n) = Flux_Vol%Xv(v5,n) + Flux_Vol%Xv(v6,n)
       end do

       Screwed_Volume = Screwed_Volume + X1(1)*(X2(2)*X3(3) - X3(2)*X2(3)) &
            +                            X1(2)*(X3(1)*X2(3) - X2(1)*X3(3)) &
            +                            X1(3)*(X2(1)*X3(2) - X3(1)*X2(2))

    end do

    Screwed_Volume = Screwed_Volume / 12.0_r8

  end function screwed_volume

end module flux_volume_module
