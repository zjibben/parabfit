module surface_type_test

  use kinds, only: r8
  use surface_type
  implicit none
  private

  public :: surface_unit_test

contains

  subroutine surface_unit_test ()

    use polygon_type
    use polyhedron_type
    use plane_type
    use array_utils, only: normalize
    use hex_types,   only: hex_f,hex_e,cube_v

    type(surface)    :: surf
    type(polyhedron) :: cube
    type(plane)      :: P
    type(polygon)    :: element
    integer          :: f,nV, ierr

    ! generate a cube and use its faces as the definition of a surface, then print to file
    call cube%init (ierr, cube_v, hex_f, hex_e)

    do f = 1,cube%parent%nfaces
      nV = count(cube%parent%face_vid(:,f) /= 0) ! number of vertices on this face
      call element%init (cube%parent%x(:,cube%parent%face_vid(1:nV,f)))

      !call surf%append (element, f)
    end do

    call surf%write_ply ('surf.ply')

    call surf%purge ()
    P%normal = [4.0_r8, 1.0_r8, 1.0_r8]
    P%normal = normalize(P%normal)
    P%rho    = 0.5_r8 / sqrt(3.0_r8)
    !call surf%append (cube%intersection_verts (P), 1)
    call surf%write_ply ('surf.ply')

  end subroutine surface_unit_test

end module surface_type_test
