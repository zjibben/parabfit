module polygon_type_test

  use kinds, only: r8
  use polygon_type
  implicit none
  private

  public :: polygon_unit_test

contains

  subroutine polygon_unit_test ()

    use array_utils, only: xrange

    type(polygon) :: pg
    integer, allocatable :: verts(:)

    write(*,*)
    write(*,*) 'POLYGON'
    write(*,*) '===================================================='

    ! call pg%init (reshape([&
    !     5.9375E-01_r8, 6.1179385469195579627665893E-01_r8, 0.0_r8,&
    !     5.9375E-01_r8, 6.1179385471097824655828390E-01_r8, 0.0_r8,&
    !     6.2500E-01_r8, 6.1455516841502266789376563E-01_r8, 0.0_r8,&
    !     6.2500E-01_r8, 6.1455516843936131010650570E-01_r8, 0.0_r8],&
    !     [3,4]))

    ! verts = xrange(1,4)
    ! call pg%order (verts)

    ! call pg%print_data ()
    ! write(*,*) verts

    call pg%init (reshape([&
        3.5830184068639514E-01_r8, 3.750023012842075842E-01_r8, 2.65625000000000000E-01_r8,&
        3.5830184068639514E-01_r8, 3.750000000000000000E-01_r8, 2.65537669284851618E-01_r8,&
        3.5830184068639514E-01_r8, 3.750023012853020421E-01_r8, 2.65625000000000000E-01_r8], [3,3]))
    call pg%order ()

    !call LS_fatal ("temporarily killed testing")

    write(*,*) '===================================================='
    write(*,*)

  end subroutine polygon_unit_test

end module polygon_type_test
