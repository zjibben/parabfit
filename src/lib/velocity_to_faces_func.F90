!!
!! VELOCITY_TO_FACES
!!
!! This module contains the subroutine, and supporting subroutines,
!! to project cell centered velocities to face velocites
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! June 2015
!!

module velocity_to_faces_func
  use kinds, only: r8
  use logging_services
  implicit none
  private

  public :: velocity_to_faces

contains

  subroutine velocity_to_faces(fluxing_velocity, velocity_cc, mesh, gmesh, t, prescribed, prescribed_case)
    use unstr_mesh_type
    use prescribed_velocity_fields
    use mesh_geom_type
    implicit none

    real(r8),          intent(inout)        :: fluxing_velocity(:,:)
    real(r8),          intent(in)           :: velocity_cc(:,:)
    class(unstr_mesh), intent(in)           :: mesh
    class(mesh_geom),  intent(in)           :: gmesh
    real(r8),          intent(in)           :: t
    logical,           intent(in), optional :: prescribed
    integer,           intent(in), optional :: prescribed_case

    ! local variables
    logical  :: prescribedh
    integer  :: f,n
    real(r8) :: face_center(3),cell_center(3),outward_normal(3)

    if (present(prescribed)) then
      prescribedh = prescribed
    else
      prescribedh = .false.
    end if

    if (prescribedh) then
      if (.not.present(prescribed_case)) call LS_fatal ('velocity_to_faces expected prescribed field case')

      !$omp parallel do default(private) shared(fluxing_velocity,mesh,gmesh,t,prescribed_case)
      do n = 1,mesh%ncell
        do f = 1,6
          ! face center is the average of the face node positions
          face_center = sum(mesh%x(:,mesh%fnode(:,mesh%cface(f,n))), dim=2) / 4.0_r8 
          
          fluxing_velocity(f,n) = sum(prescribed_velocity (face_center, t, prescribed_case) * gmesh%outnorm(:,f,n) )
        end do
      end do
      !$omp end parallel do

    else
      call LS_fatal ('need to code center to face projection routine')
    end if

  end subroutine velocity_to_faces

  ! subroutine calculate_fluxing_velocities (fluxing_velocity, velocity_cc, mesh, t, prescribed, prescribed_case)
  !   ! calculates velocity magnitudes normal to faces
  !   use unstr_mesh_type
  !   use prescribed_velocity_fields
  !   implicit none
    
  !   real(r8), dimension(  :), intent(inout) :: fluxing_velocity
  !   real(r8), dimension(:,:), intent(in)    :: velocity_cc
  !   class(unstr_mesh)       , intent(in)    :: mesh
  !   real(r8)                , intent(in)    :: t
  !   logical, intent(in), optional :: prescribed
  !   integer, intent(in), optional :: prescribed_case

  !   ! local variables
  !   logical :: prescribedh
  !   integer :: i
    
  !   if (present(prescribed)) then
  !      prescribedh = prescribed
  !   else
  !      prescribedh = .false.
  !   end if

  !   if (prescribedh) then
  !      if (.not.present(prescribed_case)) call LS_fatal ('velocity_to_faces expected prescribed field case')
  !      do i = 1,mesh%ncell
  !         fluxing_velocity(i) = sum( prescribed_velocity (mesh%x(:,i), t, prescribed_case) * mesh%normal(:,i)/mesh%area(i) )
  !      end do
  !      return
  !   end if
    
  ! end subroutine calculate_fluxing_velocities
  
end module velocity_to_faces_func
