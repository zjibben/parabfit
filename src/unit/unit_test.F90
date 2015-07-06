program unit_test
  use locate_plane_module
  use volume_track_module, only: material_flux_unit_test
  implicit none

  type(locate_plane_hex) :: locate_plane

  ! vof plane reconstruction unit test
  call locate_plane%unit_test ()

  call material_flux_unit_test (locate_plane)

end program unit_test
