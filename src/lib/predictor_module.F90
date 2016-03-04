!
! Define the Navier-Stokes predictor step
!

module predictor_module  

  use kinds, only: r8
  implicit none
  private

  public :: predictor

contains

  ! note 1: This is a crude way to account for solid material within the cell.
  !         In SOLVE_FOR_VELOCITY we also divide by FluidVof to account for
  !         the mass of fluid in the cell, more or less canceling out this
  !         term. Momentum advection is specifically excluded because the VOF
  !         is already accounting for the solid material.
  subroutine predictor (velocity_cc, gradP_dynamic_over_rho_cc, dt, material_density, volume_flux, fluidRho, &
      fluidRho_n, fluidVof, fluidVof_n, velocity_cc_n, fluxing_velocity, &
      inviscid, mesh, gmesh)

    use consts, only: ndim,nfc
    use unstr_mesh_type
    use mesh_geom_type

    real(r8),         intent(inout) :: velocity_cc(:,:)
    real(r8),         intent(in)    :: gradP_dynamic_over_rho_cc(:,:), dt, material_density(:), &
        volume_flux(:,:,:), fluidRho(:), fluidRho_n(:), fluidVof(:), fluidVof_n(:), &
        velocity_cc_n(:,:), fluxing_velocity(:,:)
    logical,          intent(in)    :: inviscid
    type(unstr_mesh), intent(in)    :: mesh
    type(mesh_geom),  intent(in)    :: gmesh

    real(r8) :: dMomentum(ndim,mesh%ncell), mass(ndim), mass_flux
    integer :: i, i_ngbr, f

    do i = 1,mesh%ncell
      dMomentum(:,i) = 0.0_r8

      ! pressure
      dMomentum(:,i) = dMomentum(:,i) - dt*fluidRho(i)*gradP_dynamic_over_rho_cc(:,i)

      ! TODO: explicit viscous stress
      
      ! TODO: turbulence
      
      ! multiply pressure gradient and viscous stress by fluidvof (see note 1)
      dMomentum(:,i) = dMomentum(:,i) * fluidVof(i)

      ! advect momentum
      do f = 1,nfc
        mass_flux = sum(material_density*volume_flux(:,f,i))
        if (mass_flux > 0.0_r8) then ! outflow
          dMomentum(:,i) = dMomentum(:,i) - mass_flux*velocity_cc(:,i) / mesh%volume(i)
        else ! inflow
          i_ngbr = gmesh%cneighbor(f,i)
          if (i_ngbr > 0) then ! not a boundary
            dMomentum(:,i) = dMomentum(:,i) - mass_flux*velocity_cc_n(:,i_ngbr) / mesh%volume(i)
          else
            dMomentum(:,i) = dMomentum(:,i) &
                - mass_flux*fluxing_velocity(f,i)*gmesh%outnorm(:,f,i) / mesh%volume(i)
          end if
        end if
      end do
      
      ! TODO: viscosity
      
      ! TODO: solidifying flow

      ! TODO: surface tension

      ! TODO: porous drag
      
      ! TODO: body forces

      ! final rhs part (here it is no longer dmomentum, but the new momentum value)
      ! TODO: * either rename dMomentum or change the behavior so it does not get the previous
      !         momentum value added in
      dMomentum(:,i) = dMomentum(:,i) + fluidRho_n(i)*fluidVof_n(i)*velocity_cc_n(:,i)
      
      ! update the velocity explicitly if there is no need for the implicit solve
      ! note that updating velocity_cc inside this loop necessitates the use of
      ! the copy velocity_cc_n above when accessing neighbors
      ! TODO: * expand the "mass matrix" with the other physics that gets added,
      !         possibly warranting breaking this entire section into its own subroutine.
      !       * do the implicit stuff after this loop
      !       * instead of a inviscid flag, just check viscosity values?
      if (inviscid) then
        mass = fluidRho(i)*fluidVof(i)
        velocity_cc(:,i) = merge(dMomentum(:,i) / mass, 0.0_r8, mask=mass>0.0_r8)
      end if
    end do

    ! TODO: solve the implicit system for viscous flows

  end subroutine predictor

end module predictor_module
