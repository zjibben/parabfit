program unit_test
  use locate_plane_module
  implicit none

  type(locate_plane_hex) :: locate_plane

  ! vof plane reconstruction unit test
  call locate_plane%unit_test ()
  
end program unit_test
