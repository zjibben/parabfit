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
      ! send cell data to the multimat_cell type
      call cell%init (mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, mesh%volume(i), gmesh%outnorm(:,:,i))

      ! partition the cell based on vofs and norms
      call cell%partition (vof(:,i), int_norm(:,:,i))

      if (dump_intrec) then
        ! TODO: dump interface reconstruction
      end if

      volume_flux_sub(:,:,i) = cell%outward_volflux (adv_dt, fluxing_velocity(:,i))
    end do

  end subroutine volume_track_nd
  
end module volume_track_nd_module
