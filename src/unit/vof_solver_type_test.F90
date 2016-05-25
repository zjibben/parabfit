module vof_solver_type_test

  use kinds,  only: r8
  use consts, only: cutvof,nfc
  use vof_solver_type
  use unstr_mesh_type
  use mesh_geom_type
  use logging_services
  use surface_type
  implicit none
  private

  public :: parallel_interfaces_test, intersecting_interfaces_test

contains

  subroutine parallel_interfaces_test ()

    use unstr_mesh_factory
    use array_utils, only: int2str

    type(unstr_mesh), target :: mesh
    type(mesh_geom),  target :: gmesh
    type(vof_solver)         :: vof_slv
    integer                  :: m

    ! generate a 3x3x1 regular mesh
    mesh = new_unstr_mesh ([0.0_r8, 0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8, 1.0_r8], [3,3,1])
    call gmesh%init (mesh)

    ! initialize the vof solver
    vof_slv%nmat = 3
    allocate(vof_slv%matl_id(vof_slv%nmat), vof_slv%vof(vof_slv%nmat,mesh%ncell), &
        vof_slv%intrec(vof_slv%nmat-1))
    vof_slv%matl_id = [1,2,3]
    vof_slv%advect_method = ONION_SKIN
    vof_slv%mesh  => mesh
    vof_slv%gmesh => gmesh

    ! initialize the vof
    vof_slv%vof(1,cell_index(1,1,1)) = 1.0_r8
    vof_slv%vof(1,cell_index(1,2,1)) = 1.0_r8
    vof_slv%vof(1,cell_index(1,3,1)) = 1.0_r8
    vof_slv%vof(1,cell_index(2,1,1)) = 1.0_r8/3.0_r8
    vof_slv%vof(1,cell_index(2,2,1)) = 1.0_r8/3.0_r8
    vof_slv%vof(1,cell_index(2,3,1)) = 1.0_r8/3.0_r8
    vof_slv%vof(1,cell_index(3,1,1)) = 0.0_r8
    vof_slv%vof(1,cell_index(3,2,1)) = 0.0_r8
    vof_slv%vof(1,cell_index(3,3,1)) = 0.0_r8

    vof_slv%vof(2,cell_index(1,1,1)) = 0.0_r8
    vof_slv%vof(2,cell_index(1,2,1)) = 0.0_r8
    vof_slv%vof(2,cell_index(1,3,1)) = 0.0_r8
    vof_slv%vof(2,cell_index(2,1,1)) = 1.0_r8/3.0_r8
    vof_slv%vof(2,cell_index(2,2,1)) = 1.0_r8/3.0_r8
    vof_slv%vof(2,cell_index(2,3,1)) = 1.0_r8/3.0_r8
    vof_slv%vof(2,cell_index(3,1,1)) = 0.0_r8
    vof_slv%vof(2,cell_index(3,2,1)) = 0.0_r8
    vof_slv%vof(2,cell_index(3,3,1)) = 0.0_r8

    vof_slv%vof(3,cell_index(1,1,1)) = 0.0_r8
    vof_slv%vof(3,cell_index(1,2,1)) = 0.0_r8
    vof_slv%vof(3,cell_index(1,3,1)) = 0.0_r8
    vof_slv%vof(3,cell_index(2,1,1)) = 1.0_r8/3.0_r8
    vof_slv%vof(3,cell_index(2,2,1)) = 1.0_r8/3.0_r8
    vof_slv%vof(3,cell_index(2,3,1)) = 1.0_r8/3.0_r8
    vof_slv%vof(3,cell_index(3,1,1)) = 1.0_r8
    vof_slv%vof(3,cell_index(3,2,1)) = 1.0_r8
    vof_slv%vof(3,cell_index(3,3,1)) = 1.0_r8

    ! reconstruct planes
    call vof_slv%update_intrec_surf ()

    ! plot plane reconstruction
    do m = 1,vof_slv%nmat-1
      call vof_slv%intrec(m)%write_ply ('par_'//trim(int2str(m))//'.ply')
    end do

    write(*,*) 'parallel interfaces reconstruction dumped'

  end subroutine parallel_interfaces_test

  subroutine intersecting_interfaces_test ()

    use unstr_mesh_factory
    use array_utils, only: int2str

    type(unstr_mesh), target :: mesh
    type(mesh_geom),  target :: gmesh
    type(vof_solver)         :: vof_slv
    integer                  :: m

    ! generate a 3x3x1 regular mesh
    mesh = new_unstr_mesh ([0.0_r8, 0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8, 1.0_r8], [3,3,1])
    call gmesh%init (mesh)

    ! initialize the vof solver
    vof_slv%nmat = 3
    allocate(vof_slv%matl_id(vof_slv%nmat), vof_slv%vof(vof_slv%nmat,mesh%ncell), &
        vof_slv%intrec(vof_slv%nmat))
    vof_slv%matl_id = [1,2,3]
    vof_slv%advect_method = ONION_SKIN
    vof_slv%mesh  => mesh
    vof_slv%gmesh => gmesh

    ! initialize the vof
    vof_slv%vof(1,cell_index(1,1,1)) = 1.0_r8
    vof_slv%vof(1,cell_index(1,2,1)) = 1.0_r8
    vof_slv%vof(1,cell_index(1,3,1)) = 1.0_r8
    vof_slv%vof(1,cell_index(2,1,1)) = 1.0_r8/2.0_r8
    vof_slv%vof(1,cell_index(2,2,1)) = 1.0_r8/2.0_r8
    vof_slv%vof(1,cell_index(2,3,1)) = 1.0_r8/2.0_r8
    vof_slv%vof(1,cell_index(3,1,1)) = 0.0_r8
    vof_slv%vof(1,cell_index(3,2,1)) = 0.0_r8
    vof_slv%vof(1,cell_index(3,3,1)) = 0.0_r8

    vof_slv%vof(2,cell_index(1,1,1)) = 0.0_r8
    vof_slv%vof(2,cell_index(1,2,1)) = 0.0_r8
    vof_slv%vof(2,cell_index(1,3,1)) = 0.0_r8
    vof_slv%vof(2,cell_index(2,1,1)) = 1.0_r8/2.0_r8
    vof_slv%vof(2,cell_index(2,2,1)) = 1.0_r8/4.0_r8
    vof_slv%vof(2,cell_index(2,3,1)) = 0.0_r8
    vof_slv%vof(2,cell_index(3,1,1)) = 1.0_r8
    vof_slv%vof(2,cell_index(3,2,1)) = 1.0_r8/2.0_r8
    vof_slv%vof(2,cell_index(3,3,1)) = 0.0_r8

    vof_slv%vof(3,cell_index(1,1,1)) = 0.0_r8
    vof_slv%vof(3,cell_index(1,2,1)) = 0.0_r8
    vof_slv%vof(3,cell_index(1,3,1)) = 0.0_r8
    vof_slv%vof(3,cell_index(2,1,1)) = 0.0_r8
    vof_slv%vof(3,cell_index(2,2,1)) = 1.0_r8/4.0_r8
    vof_slv%vof(3,cell_index(2,3,1)) = 1.0_r8/2.0_r8
    vof_slv%vof(3,cell_index(3,1,1)) = 0.0_r8
    vof_slv%vof(3,cell_index(3,2,1)) = 1.0_r8/2.0_r8
    vof_slv%vof(3,cell_index(3,3,1)) = 1.0_r8

    ! reconstruct planes
    call vof_slv%update_intrec_surf ()

    ! plot plane reconstruction
    do m = 1,vof_slv%nmat-1
      call vof_slv%intrec(m)%write_ply ('int_'//trim(int2str(m))//'.ply')
    end do

    write(*,*) 'intersecting interfaces reconstruction dumped'

  end subroutine intersecting_interfaces_test

  ! note this is a duplicate of a private function in unstr_mesh_factory
  integer function cell_index (i, j, k)
    integer, intent(in) :: i, j, k
    integer, parameter  :: nx(3) = [3,3,1]
    cell_index = i + ((j-1) + (k-1)*nx(2))*nx(1)
  end function cell_index


end module vof_solver_type_test
