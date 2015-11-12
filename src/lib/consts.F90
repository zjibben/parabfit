module consts
  use kinds, only: r8
  public

  real(r8), parameter :: alittle = epsilon(1.0_r8)
  real(r8), parameter :: alpha   = 1e-9_r8 ! this one could probably go in a less general module
  real(r8), parameter :: cutvof  = 1e-8_r8

  integer, parameter :: ndim = 3
  integer, parameter :: nfc  = 6

end module consts
