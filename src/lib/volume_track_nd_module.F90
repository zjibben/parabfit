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
  use logging_services
  implicit none
  private

  public :: volume_track_nd

contains

  subroutine volume_track_nd (volume_flux_sub, adv_dt, mesh, gmesh, vof, fluxing_velocity, fluidRho,&
       intrec, dump_intrec)

    use unstr_mesh_type
    use mesh_geom_type
    use surface_type
    use timer_tree_type
    use int_norm_module, only: interface_normal
    ! ! DEBUGGING ##############
    ! use consts, only: cutvof
    ! use curvature
    ! ! ########

    real(r8),         intent(out) :: volume_flux_sub(:,:,:)
    real(r8),         intent(in)  :: adv_dt, vof(:,:), fluxing_velocity(:,:), fluidrho(:)
    type(unstr_mesh), intent(in)  :: mesh
    type(mesh_geom),  intent(in)  :: gmesh
    type(surface),    intent(out) :: intrec(:)
    logical,          intent(in)  :: dump_intrec

    real(r8) :: int_norm(3,size(vof,dim=1),mesh%ncell),k,lnorm(3),kex
    integer  :: i,m,j

    if (dump_intrec) then
      do m = 1,size(intrec)
        call intrec(m)%purge ()
      end do
    end if

    ! compute interface normal vectors for all the materials
    int_norm = interface_normal (vof, mesh, gmesh, .false.)

    ! calculate the flux volumes for each face
    call start_timer ("reconstruct/advect")

    j=0;lnorm=0.0_r8; kex=1.0_r8/0.25_r8
    !!$omp parallel do schedule(dynamic,100)
    do i = 1,mesh%ncell
      volume_flux_sub(:,:,i) = cell_outward_volflux (mesh%x(:,mesh%cnode(:,i)), mesh%volume(i), &
          mesh%area(mesh%cface(:,i)), gmesh%outnorm(:,:,i), vof(:,i), int_norm(:,:,i), dump_intrec, &
          intrec, adv_dt, fluxing_velocity(:,i), fluidRho(i), i)
      ! if (vof(1,i) > cutvof .and. vof(1,i) < 1.0_r8-cutvof .and. all(gmesh%cneighbor(:,i)>0)) then
      !   k = meanCurvature(int_norm(:,1,:), vof(1,:), mesh%area(mesh%cface(:,i)), mesh%volume(i), i, mesh, gmesh)
      !   lnorm(1) = lnorm(1) + abs(k-kex)    !l1
      !   lnorm(2) = lnorm(2) + (k-kex)**2    !l2
      !   lnorm(3) = max(lnorm(3),abs(k-kex)) !linf
      !   j = j+1
      !   ! print '(a,i6,a,3es15.4,a,2f10.4)', "cell ",i, "   (",gmesh%xc(:,i),")  curvature: ",k, abs(k-kex)
      !   ! print *
      ! end if
      ! !if (j>=1) exit
    end do
    !!$omp end parallel do
    ! print '(4(a,es10.2))', &
    !     'done, l1:', lnorm(1)/real(j,r8), &
    !     ' l2:', sqrt(lnorm(2)/real(j,r8)), &
    !     ' linf:', lnorm(3), &
    !     ' ex:',kex
    !stop

    call stop_timer ("reconstruct/advect")

  end subroutine volume_track_nd

  function cell_outward_volflux (x, vol, farea, outnorm, vof, int_norm, dump_intrec, intrec, adv_dt,&
       fluxing_velocity, fluidRho, cell_id)

    use consts,    only: nfc
    use hex_types, only: hex_f,hex_e
    use multimat_cell_type
    use surface_type

    real(r8),      intent(in)    :: x(:,:), vol, farea(:), outnorm(:,:), vof(:), int_norm(:,:), &
        adv_dt, fluxing_velocity(:), fluidRho
    type(surface), intent(inout) :: intrec(:)
    logical,       intent(in)    :: dump_intrec
    integer,       intent(in)    :: cell_id
    real(r8)                     :: cell_outward_volflux(size(vof),nfc)

    type(multimat_cell) :: cell
    integer :: m, ierr

    ! send cell data to the multimat_cell type
    call cell%init (ierr, x, hex_f, hex_e, outnorm, vol)
    if (ierr /= 0) call LS_fatal ('cell_outward_volflux failed: could not initialize cell')

    ! partition the cell based on vofs and norms
    call cell%partition (vof, int_norm)

    ! dump the interface reconstruction, if requested
    if (dump_intrec) then
      do m = 1,size(intrec)
        call intrec(m)%append (cell%interface_polygons(m), cell_id)
      end do
    end if

    ! calculate the outward volume flux
    cell_outward_volflux = cell%outward_volflux (adv_dt, fluxing_velocity, farea, ierr)
    if (ierr /= 0) then
      write(*,*) 'cell_outward_volflux failed'
      write(*,'(a,10es20.10)') 'vof: ',vof
      write(*,'(a, 6es20.10)') 'flx: ',fluxing_velocity
      write(*,*) 'dt: ',adv_dt
      do m = 1,size(int_norm, dim=2)
        write(*,'(a,i3,3es20.10)') 'int_norm: ', m, int_norm(:,m)
      end do

      write(*,*) 'cell dimensions: '
      call cell%geom%print_data ()

      call LS_fatal ("cell_outward_volflux failed")
    end if

  end function cell_outward_volflux

end module volume_track_nd_module
