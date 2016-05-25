module volume_track_nd_test

  use kinds, only: r8
  use volume_track_nd_module
  use logging_services
  implicit none
  private

  public :: volume_track_nd_unit_test
  
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
    write(*,*) 'VOLUME TRACK ND'
    write(*,*) '===================================================='

    call volume_track_nd (volflux, 0.01_r8, mesh, gmesh, vof, fluxing_velocity, fluidRho, intrec, .false.)

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

end module volume_track_nd_test
