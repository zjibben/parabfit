!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Copyright (c) Los Alamos National Security, LLC.  This file is part of the
!! Pececillo code (C17093) and is subject to a modified BSD license; see the
!! file LICENSE located in the top-level directory of this distribution.
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! routine for calculating gravityheads for orthogonal meshes
!

module body_force_module

  use kinds, only: r8
  use logging_services
  implicit none
  private

  public :: calc_gravityhead, calc_gravityhead_bndry

  real(r8), parameter :: NULL_R = huge(1.0_r8)
  real(r8), parameter :: mech_energy_bound = NULL_R ! TODO: this should be user-set

contains

  function calc_gravityhead (fluidRho, vel_cc, xc, xf, body_force) result(gh)

    real(r8), intent(in)  :: fluidRho(:), vel_cc(:,:), xc(:,:), xf(:), body_force(:)
    real(r8)              :: gh(2)

    logical :: gl

    gl = gravity_limited(body_force, vel_cc, xc)

    gh(1) = merge(-fluidRho(1)*dot_product(xc(:,1)-xf, body_force), 0.0_r8, .not.gl)
    gh(2) = merge(-fluidRho(2)*dot_product(xc(:,2)-xf, body_force), 0.0_r8, .not.gl)

  end function calc_gravityhead

  logical function gravity_limited (body_force, vel_cc, xc)

    real(r8), intent(in) :: body_force(:), vel_cc(:,:), xc(:,:)

    if (mech_energy_bound == NULL_R) then
      gravity_limited = .false.
    else
      gravity_limited = &
          cell_mech_energy(vel_cc(:,1), xc(:,1), body_force) > mech_energy_bound .or. &
          cell_mech_energy(vel_cc(:,2), xc(:,2), body_force) > mech_energy_bound
    end if

  end function gravity_limited

  ! TODO: find a more elegant way of handling boundaries. maybe use an interface?
  function calc_gravityhead_bndry (fluidRho, vel_cc, xc, xf, body_force) result(gh)

    real(r8), intent(in)  :: fluidRho, vel_cc(:), xc(:), xf(:), body_force(:)
    real(r8)              :: gh(2)

    logical :: gl

    gl = gravity_limited_bndry (body_force, vel_cc, xc)

    gh(1) = merge(-fluidRho*dot_product(xc-xf, body_force), 0.0_r8, .not.gl)
    gh(2) = 0.0_r8

  end function calc_gravityhead_bndry

  logical function gravity_limited_bndry (body_force, vel_cc, xc)

    real(r8), intent(in) :: body_force(:), vel_cc(:), xc(:)

    if (mech_energy_bound == NULL_R) then
      gravity_limited_bndry = .false.
    else
      gravity_limited_bndry = cell_mech_energy(vel_cc, xc, body_force) > mech_energy_bound
    end if

  end function gravity_limited_bndry

  ! calculate the cell mechanical energy
  ! note the potential energy increases in the direction opposite to body_force
  real(r8) pure function cell_mech_energy (vel_cc, xc, body_force)
    use array_utils, only: normalize
    real(r8), intent(in) :: vel_cc(:), xc(:), body_force(:)
    cell_mech_energy = 0.5_r8 * sum(max(0.0_r8,vel_cc*normalize(body_force)))**2 & ! kinetic
        + (-dot_product(xc,body_force))                                            ! potential
  end function cell_mech_energy

end module body_force_module
