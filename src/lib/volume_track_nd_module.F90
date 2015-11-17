!!
!! VOLUME_TRACK_ND
!!
!! This module provides volume tracking capabilities
!! using a nested dissection method.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! Last revised 4 Nov 2012.
!!

module volume_track_nd_module
  use kinds, only: r8
  implicit none
  private

  public :: volume_track_nd, unit_test_nd

contains

  subroutine unit_test_nd () !(plane_cell, face_fluxing_velocity, fluxex)
    ! do something here!!
  end subroutine unit_test_nd

  subroutine volume_track_nd (volume_flux_sub, adv_dt, mesh, gmesh, Vof, fluxing_velocity, nmat, fluidRho, &
       intrec, dump_intrec)
    use plane_type
    use multimat_cell_type
    use unstr_mesh_type
    use mesh_geom_type
    use surface_type
    use hex_types,       only: hex_f,hex_e
    use int_norm_module, only: interface_normal
    use consts,          only: cutvof,nfc
    use logging_services ! DEBUGGING
    
    real(r8),         intent(out) :: volume_flux_sub(:,:,:)
    real(r8),         intent(in)  :: adv_dt, vof(:,:), fluxing_velocity(:,:), fluidrho(:)
    type(unstr_mesh), intent(in)  :: mesh
    type(mesh_geom),  intent(in)  :: gmesh
    integer,          intent(in)  :: nmat
    type(surface),    intent(out) :: intrec(:)
    logical,          intent(in)  :: dump_intrec

    type(multimat_cell)  :: cell
    type(plane)          :: flux_plane
    real(r8)             :: xf(3), int_norm(3,nmat,mesh%ncell)
    integer              :: i,f,ni,ninterfaces

    ninterfaces = nmat-1

    ! compute interface normal vectors for all the materials.
    ! TODO: the norm gets calculated differently for nested dissection??
    int_norm = interface_normal (vof, mesh, gmesh)

    if (dump_intrec) then
      do ni = 1,ninterfaces
        call intrec(ni)%purge ()
      end do
    end if
    
    ! calculate the flux volumes for each face
    do i = 1,mesh%ncell
    !i = 6443
      ! send cell data to the multimat_cell type
      call cell%init (mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, mesh%volume(i), gmesh%outnorm(:,:,i))
      
      ! partition the cell based on vofs and norms
      call cell%partition (vof(:,i), int_norm(:,:,i))

      if (dump_intrec) then
        ! TODO: dump interface reconstruction
      end if

      
      volume_flux_sub(:,:,i) = cell%outward_volflux (adv_dt, fluxing_velocity(:,i))

      ! ! calculate how much material is being fluxed out of each face
      ! do f = 1,nfc
      !   if (fluxing_velocity(f,i) < cutvof*mesh%volume(i)) then
      !     volume_flux_sub(:,f,i) = 0.0_r8
      !   else
      !     ! find the plane equation for the back end of the flux volume
      !     ! WARNING: in general, this could be non-planar, just like cell faces
      !     flux_plane%normal = -gmesh%outnorm(:,f,i)
      !     xf = sum(mesh%x(:,mesh%fnode(:,mesh%cface(f,i)))) / 6.0_r8 ! face center
      !     flux_plane%rho  = sum(xf*flux_plane%normal) - adv_dt * fluxing_velocity(f,i)

      !     ! find the volume of the materials behind the flux plane
      !     volume_flux_sub(:,f,i) = cell%volumes_behind_plane (flux_plane)
      !   end if
      ! end do
      
    end do

    !   !i = 6443
    !   write(*,*)
    !   write(*,*) 'vol ',mesh%volume(i)
    !   write(*,*) 'vof ',vof(:,i)
    !   write(*,*) 'vols',cell%mat_poly(1)%volume (), cell%mat_poly(2)%volume ()
    ! do f = 1,nfc
    !   write(*,('(3es12.4)')) fluxing_velocity(f,i), adv_dt * fluxing_velocity(f,i) * mesh%area(mesh%cface(f,i))
    !   write(*,('(2es12.4)')) volume_flux_sub(:,f,i)
    !   write(*,*)
    ! end do

    ! call LS_fatal ("stop here")

  end subroutine volume_track_nd

  ! function cell_volume_flux (cell, adv_dt, fluxing_velocity, mesh, cell_norm)
  !   use consts,             only: nfc,cutvof
  !   use multimat_cell_type
  !   use unstr_mesh_type
    
  !   type(multimat_cell), intent(in) :: cell
  !   real(r8),            intent(in) :: adv_dt, fluxing_velocity(:), cell_norm(:,:)
  !   type(unstr_mesh),    intent(in) :: mesh
  !   real(r8)                        :: cell_volume_flux(:,:)

  !   real(r8)             :: xf(3)
  !   type(plane)          :: flux_plane
  !   integer              :: f
    
  !   ! calculate how much material is being fluxed out of each face
  !   do f = 1,nfc
  !     if (fluxing_velocity(f) < cutvof*mesh%volume(i)) then
  !       ! no flux out this face
  !       cell_volume_flux(:,f) = 0.0_r8
  !     else
  !       ! find the plane equation for the back end of the flux volume
  !       ! WARNING: in general, this could be non-planar, just like cell faces
  !       flux_plane%normal = -cell_norm(:,f)
  !       xf = sum(mesh%x(:,mesh%fnode(:,mesh%cface(f,i)))) / 6.0_r8 ! face center
  !       flux_plane%rho  = sum(xf*flux_plane%normal) + adv_dt * fluxing_velocity(f)

  !       ! find the volume of the materials behind the flux plane
  !       cell_volume_flux(:,f) = cell%volumes_behind_plane (flux_plane)
  !     end if
  !   end do

  ! end function cell_volume_flux

end module volume_track_nd_module
