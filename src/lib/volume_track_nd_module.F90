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

  public :: volume_track_nd, volume_track_nd_unit_test

contains

  subroutine volume_track_nd_unit_test () !(plane_cell, face_fluxing_velocity, fluxex)
    use consts, only: nfc
    use unstr_mesh_type
    use unstr_mesh_factory
    use mesh_geom_type
    use surface_type

    integer, parameter :: nmat = 2, ncell = 9
    integer            :: i,j
    type(unstr_mesh)   :: mesh
    type(mesh_geom)    :: gmesh
    type(surface)      :: intrec(nmat)
    real(r8)           :: volflux(nmat,nfc,ncell), volflux_ex(nmat,nfc,ncell), &
         fluxing_velocity(nfc,ncell), vof(nmat,ncell), fluidRho(ncell)

    ! initialize a simple 2D mesh
    mesh = new_unstr_mesh ([0.0_r8, 0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8, 1.0_r8], [3,3,1])
    call gmesh%init (mesh)
    
    ! initialize a vof on that mesh
    fluidRho = 0.0_r8 ! this is currently ignored by volume_track_nd
    
    do j = 1,3
      vof(1,cell_index (1,j,1)) = 1.0_r8; vof(2,cell_index (1,j,1)) = 0.0_r8
      vof(1,cell_index (2,j,1)) = 0.5_r8; vof(2,cell_index (2,j,1)) = 0.5_r8
      vof(1,cell_index (3,j,1)) = 0.0_r8; vof(2,cell_index (3,j,1)) = 1.0_r8
    end do

    fluxing_velocity = 0.0_r8
    fluxing_velocity(3,:) = -1.0_r8
    fluxing_velocity(4,:) =  1.0_r8

    ! do some tests
    write(*,*)
    write(*,*) 'POLYHEDRON'
    write(*,*) '===================================================='

    call volume_track_nd (volflux, 0.01_r8, mesh, gmesh, vof, fluxing_velocity, nmat, fluidRho, intrec, .false.)

    volflux_ex = 0.0_r8
    
    do i = 1,3
      do j = 1,3
        write(*,'(2es14.4)') volflux(:,4,cell_index(i,j,1))
      end do
      write(*,*)
    end do

    ! do i = 1,3
    !   do j = 1,3
    !     write(*,'(2es14.4)') volflux(:,3,cell_index(i,j,1))
    !   end do
    !   write(*,*)
    ! end do

    ! do i = 1,3
    !   do j = 1,3
    !     write(*,'(2es14.4)') volflux(:,2,cell_index(i,j,1))
    !   end do
    !   write(*,*)
    ! end do

    ! do i = 1,3
    !   do j = 1,3
    !     write(*,'(2es14.4)') volflux(:,1,cell_index(i,j,1))
    !   end do
    !   write(*,*)
    ! end do
    
    write(*,*) '===================================================='
    write(*,*)

  contains
    
    ! note this is a duplicate of a private function in unstr_mesh_factory
    integer function cell_index (i, j, k)
      integer, intent(in) :: i, j, k
      integer, parameter  :: nx(3) = [3,3,1]
      cell_index = i + ((j-1) + (k-1)*nx(2))*nx(1)
    end function cell_index

  end subroutine volume_track_nd_unit_test

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
    !i=4
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

    ! i = 4
    ! write(*,*)
    ! write(*,'(a,3f14.4)') 'vof  ',vof(:,i)
    ! do f = 1,nmat
    !   write(*,'(a,3f14.4)') 'norm ',int_norm(:,f,i)
    ! end do
    ! do f = 1,nfc
    !   write(*,'(a,i3,3es14.4)') 'volflux ',f,volume_flux_sub(:,f,i)
    ! end do
    !call LS_fatal ("stop here")
    !write(*,*) 'completed iteration'
    
  end subroutine volume_track_nd
  
end module volume_track_nd_module
