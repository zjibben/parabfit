!
! Define the Navier-Stokes predictor step
!

module predictor_module  

  use kinds,  only: r8
  use consts, only: ndim,nfc
  use logging_services
  implicit none
  private

  public :: predictor

contains

  ! note 1: This is a crude way to account for solid material within the cell.
  !         In SOLVE_FOR_VELOCITY we also divide by FluidVof to account for
  !         the mass of fluid in the cell, more or less canceling out this
  !         term. Momentum advection is specifically excluded because the VOF
  !         is already accounting for the solid material.
  subroutine predictor (velocity_cc, gradP_dynamic_over_rho_cc, dt, volume_flux, &
      fluidRho, fluidRho_n, vof, fluidVof, fluidVof_n, velocity_cc_n, fluxing_velocity, &
      viscous_implicitness, mprop, inviscid, solid_face, is_pure_immobile, velocity_bc, pressure_bc, mesh, gmesh, &
      params)

    use unstr_mesh_type
    use mesh_geom_type
    use csr_matrix_type
    use hypre_hybrid_type
    use bndry_func_class
    use matl_props_type

    real(r8),          intent(inout) :: velocity_cc(:,:)
    real(r8),          intent(in)    :: gradP_dynamic_over_rho_cc(:,:), dt, &
        volume_flux(:,:,:), fluidRho(:), fluidRho_n(:), vof(:,:), fluidVof(:), fluidVof_n(:), &
        velocity_cc_n(:,:), fluxing_velocity(:,:), viscous_implicitness
    type(matl_props),  intent(in)    :: mprop
    logical,           intent(in)    :: inviscid, solid_face(:), is_pure_immobile(:)
    class(bndry_func), intent(in)    :: velocity_bc, pressure_bc
    type(unstr_mesh),  intent(in)    :: mesh
    type(mesh_geom),   intent(in)    :: gmesh
    type(parameter_list), pointer    :: params

    real(r8) :: rhs(ndim*mesh%ncell), mass(ndim), mass_flux, cell_rhs(ndim), &
        vel_cc_new(ndim*mesh%ncell), &
        grad_vel_face_bndry(ndim,ndim,mesh%nface), viscosity_face_bndry(mesh%nface)
    integer :: i, i_ngbr, f, stat
    type(csr_matrix) :: lhs
    type(hypre_hybrid) :: solver

    ! DEBUGGING #################
    real(r8) :: tmp(ndim)
    ! ###########################

    grad_vel_face_bndry = 0.0_r8
    viscosity_face_bndry = 0.0_r8

    ! initialize the sparse left hand side matrix
    if (.not.inviscid .and. viscous_implicitness > 0.0_r8) then
      call init_disc (lhs, mesh%ncell, gmesh%cneighbor)
      call lhs%set_all (0.0_r8)
    end if
    rhs = 0.0_r8
    
    ! update boundary face values prior to the rest of the domain
    call apply_velocity_bcs (grad_vel_face_bndry, viscosity_face_bndry, rhs, velocity_bc, &
        velocity_cc, dt, viscous_implicitness, inviscid, mesh, gmesh)
    call apply_pressure_bcs (viscosity_face_bndry, rhs, pressure_bc, &
        velocity_cc, dt, viscous_implicitness, inviscid, mesh, gmesh)
    
    do i = 1,mesh%ncell
      cell_rhs = 0.0_r8

      ! pressure
      cell_rhs = cell_rhs - dt*fluidRho(i)*gradP_dynamic_over_rho_cc(:,i)

      ! TODO: viscosity
      if (.not.inviscid) &
          call apply_viscosity (cell_rhs, lhs, mprop, dt, viscous_implicitness, velocity_cc_n, vof, &
          fluidVof, fluidRho(i), grad_vel_face_bndry(:,:,mesh%cface(:,i)), &
          viscosity_face_bndry(mesh%cface(:,i)), mesh%area(mesh%cface(:,i)), mesh%volume(i), &
          solid_face(mesh%cface(:,i)), is_pure_immobile(i), i, mesh%cface(:,i), gmesh)

      ! TODO: turbulence
      
      ! account for solid material (see note 1)
      cell_rhs = cell_rhs * fluidVof(i)

      ! advect momentum
      do f = 1,nfc
        mass_flux = sum(mprop%density*volume_flux(:,f,i))
        ! upwind to determine which velocity to grab
        if (mass_flux > 0.0_r8) then ! outflow
          cell_rhs = cell_rhs - mass_flux*velocity_cc_n(:,i) / mesh%volume(i)
        else ! inflow
          i_ngbr = gmesh%cneighbor(f,i)
          if (i_ngbr > 0) then ! not a boundary
            cell_rhs = cell_rhs - mass_flux*velocity_cc_n(:,i_ngbr) / mesh%volume(i)
          else
            cell_rhs = cell_rhs &
                - mass_flux*fluxing_velocity(f,i)*gmesh%outnorm(:,f,i) / mesh%volume(i)
          end if
        end if
      end do
      
      ! TODO: solidifying flow

      ! TODO: surface tension

      ! TODO: porous drag
      
      ! TODO: body forces

      ! final rhs part (here it is no longer dmomentum, but the new momentum value)
      cell_rhs = cell_rhs + fluidRho_n(i)*fluidVof_n(i)*velocity_cc_n(:,i)

      ! if (cell_rhs(2) < 0.0_r8) then
      !   tmp = 0.0_r8
      !   do f = 1,nfc
      !     mass_flux = sum(mprop%density*volume_flux(:,f,i))
      !     ! upwind to determine which velocity to grab
      !     if (mass_flux > 0.0_r8) then ! outflow
      !       tmp = tmp - mass_flux*velocity_cc_n(:,i) / mesh%volume(i)
      !     else ! inflow
      !       i_ngbr = gmesh%cneighbor(f,i)
      !       if (i_ngbr > 0) then ! not a boundary
      !         tmp = tmp - mass_flux*velocity_cc_n(:,i_ngbr) / mesh%volume(i)
      !       else
      !         tmp = tmp &
      !             - mass_flux*fluxing_velocity(f,i)*gmesh%outnorm(:,f,i) / mesh%volume(i)
      !       end if
      !     end if
      !   end do


      !   write(*,*) 'negative predictor'
      !   write(*,*) i, cell_rhs(2)
      !   write(*,*) 'convc', tmp
      !   write(*,*) 'mom_n', fluidRho_n(i)*fluidVof_n(i)*velocity_cc_n(:,i)
      !   write(*,*) 'press', -dt*fluidRho(i)*gradP_dynamic_over_rho_cc(:,i)
      !   tmp = viscousExplicit (dt, viscous_implicitness, velocity_cc_n, grad_vel_face_bndry(:,:,mesh%cface(:,i)), &
      !       viscosityFaces (mprop, vof, fluidVof, viscosity_face_bndry(mesh%cface(:,i)), solid_face(mesh%cface(:,i)), i, gmesh%cneighbor(:,i)), &
      !       mesh%area(mesh%cface(:,i)), mesh%volume(i), i, gmesh)
      !   write(*,*) 'visco', tmp
      !   stop
      ! end if
      
      ! update the velocity explicitly if there is no need for the implicit solve
      ! note that updating velocity_cc inside this loop necessitates the use of
      ! the copy velocity_cc_n above when accessing neighbors
      ! TODO: * expand the "mass matrix" with the other physics that gets added,
      !         possibly warranting breaking this entire section into its own subroutine.
      !       * do the implicit stuff after this loop
      if (inviscid .or. viscous_implicitness == 0.0_r8) then
        mass = fluidRho(i)*fluidVof(i)
        velocity_cc(:,i) = merge(cell_rhs / mass, 0.0_r8, mask=mass>0.0_r8)
      else
        rhs((i-1)*ndim+1:(i-1)*ndim+ndim) = rhs((i-1)*ndim+1:(i-1)*ndim+ndim) + cell_rhs
      end if
    end do

    ! solve the implicit system for viscous flows
    if (.not.inviscid .and. viscous_implicitness > 0.0_r8) then
      ! initial guess is current velocity (TODO: better initial guess?)
      vel_cc_new = reshape(velocity_cc, [ndim*mesh%ncell])

      call solver%init (lhs, params)
      call solver%setup ()
      call solver%solve (rhs, vel_cc_new, stat)
      if (stat /= 0) call LS_fatal ("projection solver failed")

      ! copy result into velocity_cc (TODO: can I do all this directly without copying?)
      velocity_cc = reshape(vel_cc_new, [ndim, mesh%ncell])
    end if

  contains

    subroutine init_disc (matrix,ncell,cneighbor)

      use csr_matrix_type

      type(csr_matrix), intent(out) :: matrix
      integer, intent(in) :: ncell,cneighbor(:,:)

      type(csr_graph), pointer :: g
      integer :: i,n,m, index

      allocate(g)
      call g%init (ncell*ndim)
      do i = 1,ncell
        do n = 1,ndim
          index = (i-1)*ndim+n
          call g%add_edge (index,index)
          do m = 1,ndim
            call g%add_edge (index, (pack(cneighbor(:,i), mask=cneighbor(:,i)>0) - 1) * ndim + m)
          end do
        end do
      end do
      call g%add_complete ()

      call matrix%init (g, take_graph=.true.)

    end subroutine init_disc

    subroutine apply_velocity_bcs (grad_vel_face_bndry, viscosity_face_bndry, rhs, velocity_bc, &
        velocity_cc, dt, viscous_implicitness, inviscid, mesh, gmesh)

      use bndry_func_class
      use mesh_geom_type
      use differential_operators, only: faceGradient

      real(r8),          intent(inout) :: grad_vel_face_bndry(:,:,:), viscosity_face_bndry(:), rhs(:)
      class(bndry_func), intent(in)    :: velocity_bc
      real(r8),          intent(in)    :: velocity_cc(:,:), dt, viscous_implicitness
      logical,           intent(in)    :: inviscid
      type(unstr_mesh),  intent(in)    :: mesh
      type(mesh_geom),   intent(in)    :: gmesh

      real(r8) :: velocity_bc_value, dx(ndim)
      integer  :: bndry_f, fid, i, f, n


      ! update values on boundary faces
      do bndry_f = 1,size(velocity_bc%index)
        fid = velocity_bc%index(bndry_f)
        velocity_bc_value = velocity_bc%value(bndry_f)
        
        i = gmesh%fcell(1,fid) ! id of cell attached to this face
        f = gmesh%flid(1,fid)  ! local id of the face

        if (.not.inviscid) then
          ! viscosity
          viscosity_face_bndry(fid) = viscosityCell(mprop, vof(:,i), fluidVof(i))

          ! explicit
          if (viscous_implicitness < 1.0_r8) then
            do n = 1,ndim
              grad_vel_face_bndry(:,n,fid) = faceGradient (&
                  [velocity_cc(n,i), velocity_bc_value*gmesh%outnorm(n,f,i)], &
                  reshape([gmesh%xc(:,i), gmesh%fc(:,fid)], [ndim,2]))
            end do
          end if

          ! implicit
          ! TODO: allow tangential components to dirichlet BCs
          if (viscous_implicitness > 0.0_r8) then
            dx = gmesh%xc(:,i) - gmesh%fc(:,fid)
            rhs((i-1)*ndim+1:(i-1)*ndim+ndim) = rhs((i-1)*ndim+1:(i-1)*ndim+ndim) &
                + dt*viscous_implicitness &
                *viscosity_face_bndry(fid)*velocity_bc_value*gmesh%outnorm(:,f,i) &
                *dot_product(gmesh%outnorm(:,f,i), dx) / sum(dx**2) / mesh%volume(i)
          end if
        end if
      end do

    end subroutine apply_velocity_bcs

    subroutine apply_pressure_bcs (viscosity_face_bndry, rhs, pressure_bc, &
        velocity_cc, dt, viscous_implicitness, inviscid, mesh, gmesh)

      use bndry_func_class
      use mesh_geom_type

      real(r8),          intent(inout) :: viscosity_face_bndry(:), rhs(:)
      class(bndry_func), intent(in)    :: pressure_bc
      real(r8),          intent(in)    :: velocity_cc(:,:), dt, viscous_implicitness
      logical,           intent(in)    :: inviscid
      type(unstr_mesh),  intent(in)    :: mesh
      type(mesh_geom),   intent(in)    :: gmesh

      real(r8) :: pressure_bc_value, dx(ndim)
      integer  :: bndry_f, fid, i, f

      ! update values on boundary faces
      do bndry_f = 1,size(pressure_bc%index)
        fid = pressure_bc%index(bndry_f)
        pressure_bc_value = pressure_bc%value(bndry_f)

        i = gmesh%fcell(1,fid) ! id of cell attached to this face
        f = gmesh%flid(1,fid)  ! local id of the face

        if (.not.inviscid) then
          viscosity_face_bndry(fid) = viscosityCell(mprop, vof(:,i), fluidVof(i))

          ! zero velocity gradient at Dirichlet pressure boundaries
                    
          ! explicit
          ! no additional component to the rhs
          !face_gradient ( reshape([gmesh%xc(:,i), gmesh%fc(:,fid)], [ndim,2]))

          ! implicit
          if (viscous_implicitness > 0.0_r8) then
            dx = gmesh%xc(:,i) - gmesh%fc(:,fid)
            rhs((i-1)*ndim+1:(i-1)*ndim+ndim) = rhs((i-1)*ndim+1:(i-1)*ndim+ndim) - &
                dt*viscous_implicitness*viscosity_face_bndry(fid)*velocity_cc(:,i) &
                *mesh%area(fid)*dot_product(gmesh%outnorm(:,f,i), dx) / sum(dx**2) / mesh%volume(i)
          end if
        end if
      end do

    end subroutine apply_pressure_bcs

  end subroutine predictor

  subroutine apply_viscosity (rhs, lhs, mprop, dt, viscous_implicitness, velocity_cc, vof, &
      fluidVof, fluidRho, grad_vel_face_bndry, viscosity_face_bndry, face_area, cell_vol, &
      solid_face, is_pure_immobile, i, cface, gmesh)
    
    use csr_matrix_type
    use mesh_geom_type
    use matl_props_type

    real(r8),         intent(inout) :: rhs(:)
    type(csr_matrix), intent(inout) :: lhs
    type(matl_props), intent(in)    :: mprop
    real(r8),         intent(in)    :: dt, viscous_implicitness, velocity_cc(:,:), vof(:,:), &
        fluidVof(:), fluidRho, grad_vel_face_bndry(:,:,:), viscosity_face_bndry(:), face_area(:), cell_vol
    logical,          intent(in)    :: solid_face(:), is_pure_immobile
    integer,          intent(in)    :: i, cface(:)
    type(mesh_geom),  intent(in)    :: gmesh

    real(r8) :: dx(ndim), viscosity_face(nfc), tmp
    integer  :: f,n, index, index_ngbr, i_ngbr

    ! calculate face viscosities
    viscosity_face = viscosityFaces (mprop, vof, fluidVof, viscosity_face_bndry, solid_face, i, &
        gmesh%cneighbor(:,i))

    ! explicit part
    if (viscous_implicitness < 1.0_r8) &
        rhs = rhs + viscousExplicit (dt, viscous_implicitness, velocity_cc, grad_vel_face_bndry, &
        viscosity_face, face_area, cell_vol, i, gmesh)

    ! implicit part
    if (viscous_implicitness > 0.0_r8) then
      if (is_pure_immobile .or. fluidRho == 0.0_r8) then
        ! TODO: set this row such that Ax = x
      else
        do n = 1,ndim
          index = (i-1)*ndim+n

          call lhs%increment (index,index, fluidRho)

          do f = 1,nfc
            i_ngbr = gmesh%cneighbor(f,i)
            index_ngbr = (i_ngbr-1)*ndim+n

            if (i_ngbr > 0) then
              dx = gmesh%xc(:,i) - gmesh%xc(:,i_ngbr)
            else
              dx = gmesh%xc(:,i) - gmesh%fc(:,cface(f))
            end if

            tmp = dt * viscous_implicitness &
                * viscosity_face(f) * face_area(f) * dot_product(gmesh%outnorm(:,f,i), dx) &
                / sum(dx**2) / cell_vol

            if (i_ngbr > 0) call lhs%increment (index,index_ngbr, tmp)
            call lhs%increment (index,index, -tmp)
          end do
        end do
      end if
    end if

  end subroutine apply_viscosity

  function viscousExplicit (dt, viscous_implicitness, velocity_cc, &
      grad_vel_face_bndry, viscosity_face, face_area, cell_vol, i, gmesh) result(dMomentum)

    use mesh_geom_type

    real(r8),        intent(in) :: dt, viscous_implicitness, velocity_cc(:,:), &
        grad_vel_face_bndry(:,:,:), viscosity_face(:), face_area(:), cell_vol
    integer,         intent(in) :: i
    type(mesh_geom), intent(in) :: gmesh
    real(r8)                    :: dMomentum(ndim)

    dMomentum = dt * (1.0_r8 - viscous_implicitness) &
        * divStress(velocity_cc, grad_vel_face_bndry, viscosity_face, face_area, cell_vol, i, gmesh)

  end function viscousExplicit

  function divStress (velocity_cc, grad_vel_face_bndry, viscosity_face, face_area, cell_vol, i, gmesh)

    use mesh_geom_type
    use differential_operators, only: faceGradient, divergence

    real(r8),        intent(in) :: velocity_cc(:,:), grad_vel_face_bndry(:,:,:), viscosity_face(:), &
        face_area(:), cell_vol
    integer,         intent(in) :: i
    type(mesh_geom), intent(in) :: gmesh
    real(r8)                    :: divStress(ndim)

    real(r8) :: grad_vel_face_out(nfc)
    integer  :: f, n, i_ngbr

    do n = 1,ndim
      ! calculate the outward component of the gradient of the nth velocity component on every face
      do f = 1,nfc
        i_ngbr = gmesh%cneighbor(f,i)
        if (i_ngbr > 0) then
          grad_vel_face_out(f) = dot_product( &
              faceGradient (velocity_cc(n,[i, i_ngbr]), gmesh%xc(:,[i, i_ngbr])), &
              gmesh%outnorm(:,f,i))
        else ! boundary face gradient already calculated
          grad_vel_face_out(f) = dot_product(grad_vel_face_bndry(:,n,f), gmesh%outnorm(:,f,i))
        end if
      end do

      ! calculate the divergence of the viscosity * the gradient of the nth velocity component
      divStress(n) = divergence (viscosity_face * grad_vel_face_out, face_area, cell_vol)

      ! if (n==2 .and. i==1) then ! .and. divStress(2) > 1e-12_r8) then
      !   write(*,*) 'divstress cell', i
      !   write(*,*) divstress(2)
      !   write(*,*) gmesh%cneighbor(:,i) < 1
      !   write(*,'(6es14.4)') grad_vel_face_out
      !   write(*,'(a,6es14.4)') 'mu', viscosity_face

      !   write(*,*) "myvel      ",velocity_cc(2,i)
      !   do f = 1,nfc
      !     i_ngbr = gmesh%cneighbor(f,i)
      !     if (i_ngbr > 0) then
      !       write(*,'(i,2es14.4,a,3es14.4)') f, velocity_cc(2,i_ngbr), grad_vel_face_out(f), "              ",gmesh%outnorm(:,f,i)
      !     else
      !       write(*,'(i,6es14.4)') f,grad_vel_face_bndry(:,n,f), gmesh%outnorm(:,f,i)
      !     end if
      !   end do

      !   write(*,*)
      !   !stop
      ! end if
    end do

  end function divStress

  ! calculates the viscosity on all faces for a given cell
  pure function viscosityFaces (mprop, vof, fluidVof, viscosity_face_bndry, solid_face, i, cneighbor)

    use matl_props_type

    type(matl_props), intent(in) :: mprop
    real(r8),         intent(in) :: vof(:,:), fluidVof(:), viscosity_face_bndry(:)
    logical,          intent(in) :: solid_face(:)
    integer,          intent(in) :: i, cneighbor(:)
    real(r8)                     :: viscosityFaces(nfc)

    real(r8) :: viscosity_cc, viscosity_ngbr
    integer :: f, i_ngbr

    viscosity_cc = viscosityCell(mprop, vof(:,i), fluidVof(i))
    do f = 1,nfc
      i_ngbr = cneighbor(f)
      if (i_ngbr > 0) then
        viscosity_ngbr = viscosityCell(mprop, vof(:,i_ngbr), fluidVof(i_ngbr))
        viscosityFaces(f) = viscosityFace([viscosity_cc, viscosity_ngbr], solid_face(f))
      else
        viscosityFaces(f) = viscosity_face_bndry(f)
      end if
    end do

  end function viscosityFaces

  ! calculates viscosity on the face from the two neighboring cell-centered viscosities
  real(r8) pure function viscosityFace (viscosity_cc, solid_face)

    use array_utils, only: meanHarmonic

    real(r8), intent(in) :: viscosity_cc(:)
    logical,  intent(in) :: solid_face

    viscosityFace = merge(meanHarmonic(viscosity_cc), viscosity_cc(1), mask=.not.solid_face)

  end function viscosityFace

  ! calculates the cell-centered viscosity from the material viscosities and the vof in the cell
  real(r8) pure function viscosityCell (mprop, vof, fluidVof)

    use matl_props_type

    type(matl_props), intent(in) :: mprop
    real(r8),         intent(in) :: vof(:), fluidVof

    viscosityCell = sum(mprop%viscosity*vof, mask=.not.mprop%is_immobile) / fluidVof

  end function viscosityCell

end module predictor_module
