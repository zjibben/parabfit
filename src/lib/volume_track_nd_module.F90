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

  subroutine unit_test_nd ()
    ! do something here!!
  end subroutine unit_test_nd

  subroutine volume_track_nd (volume_flux_sub, adv_dt, mesh, gmesh, Vof, fluxing_velocity, nmat, fluidRho, &
       intrec, dump_intrec)
    use plane_type
    use multimat_cell_type
    use unstr_mesh_type
    use mesh_geom_type
    use surface_type
    use int_norm_module, only: interface_normal
    
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
    integer              :: i,f

    ! compute interface normal vectors for all the materials.
    ! does the normal get calculated differently for nested dissection?
    int_norm = interface_normal (vof, mesh, gmesh)

    ! calculate the flux volumes for each face
    do i = 1,mesh%ncell
      ! send cell data to the multimat_cell type
      call cell%init (mesh%x(:,mesh%cnode(:,i)), mesh%volume(i), mesh%area(mesh%cface(:,i)), &
           mesh%normal(:,mesh%cface(:,i)), mesh%cfpar(i))

      ! partition the cell based on vofs and norms
      call cell%partition (vof(:,i), int_norm(:,:,i))

      if (dump_intrec) then
        ! do something here!!
      end if

      ! calculate how much material is being fluxed out of each face
      do f = 1,6
        ! find the plane equation for the back end of the flux volume
        ! WARNING: in general, this plane could be non-planar, just like cell faces
        flux_plane%normal = -gmesh%outnorm(:,f,i)
        xf = sum(mesh%x(:,mesh%fnode(:,mesh%cface(f,i)))) / 6.0_r8 ! face center
        flux_plane%rho  = sum(xf*flux_plane%normal) + adv_dt * fluxing_velocity(f,i)

        ! find the volume of the materials behind the flux plane
        volume_flux_sub(:,f,i) = cell%volumes_behind_plane (flux_plane)
      end do

    end do
    
  end subroutine volume_track_nd

end module volume_track_nd_module
