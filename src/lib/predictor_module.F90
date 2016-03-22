!
! Define the Navier-Stokes predictor step
!

module predictor_module  

  use kinds,  only: r8
  use consts, only: ndim,nfc
  use matl_props_type
  use bndry_func_class
  use unstr_mesh_type
  use mesh_geom_type
  use parameter_list_type
  use csr_matrix_type
  use logging_services
  use timer_tree_type
  implicit none
  private

  type, public :: predictor_solver
    private
    
    !! intent inout state variables
    real(r8), pointer :: velocity_cc(:,:)

    !! intent in state variables
    real(r8), pointer :: gradP_dynamic_over_rho_cc(:,:), dt, volume_flux(:,:,:), fluidRho(:), &
        fluidRho_n(:), vof(:,:), fluidVof(:), fluidVof_n(:), velocity_cc_n(:,:), &
        fluxing_velocity(:,:), viscous_implicitness
    type(matl_props),  pointer :: mprop
    logical,           pointer :: inviscid, solid_face(:), is_pure_immobile(:)
    class(bndry_func), pointer :: velocity_bc, pressure_bc
    type(unstr_mesh),  pointer :: mesh
    type(mesh_geom),   pointer :: gmesh
    type(parameter_list), pointer    :: params

    !! locally owned persistent variables
    type(csr_graph), pointer :: g
    type(csr_matrix) :: lhs
  contains
    procedure :: init => init_predictor_solver
    procedure :: solve
    procedure :: update_matrix
    final :: predictor_solver_delete
  end type predictor_solver

contains

  subroutine init_predictor_solver (this, velocity_cc, gradP_dynamic_over_rho_cc, volume_flux, &
      fluidRho, vof, fluidVof, fluxing_velocity, viscous_implicitness, mprop, inviscid, &
      velocity_bc, pressure_bc, mesh, gmesh, params)

    class(predictor_solver), intent(out) :: this
    real(r8),             intent(in), target :: velocity_cc(:,:), gradP_dynamic_over_rho_cc(:,:), &
        volume_flux(:,:,:), fluidRho(:), vof(:,:), fluidVof(:), fluxing_velocity(:,:), &
        viscous_implicitness
    type(matl_props),     intent(in), target :: mprop
    logical,              intent(in), target :: inviscid
    class(bndry_func),    intent(in), target :: velocity_bc, pressure_bc
    type(unstr_mesh),     intent(in), target :: mesh
    type(mesh_geom),      intent(in), target :: gmesh
    type(parameter_list), intent(in), target :: params

    this%velocity_cc => velocity_cc
    this%gradP_dynamic_over_rho_cc => gradP_dynamic_over_rho_cc
    this%volume_flux => volume_flux
    this%fluidRho => fluidRho
    this%vof => vof
    this%fluidVof => fluidVof
    this%fluxing_velocity => fluxing_velocity
    this%viscous_implicitness => viscous_implicitness
    this%mprop => mprop
    this%inviscid => inviscid
    this%velocity_bc => velocity_bc
    this%pressure_bc => pressure_bc
    this%mesh => mesh
    this%gmesh => gmesh
    this%params => params

    call this%update_matrix ()

  end subroutine init_predictor_solver

  subroutine update_matrix (this)

    class(predictor_solver), intent(inout) :: this

    integer :: i,n,m, index

    if (associated(this%g)) deallocate(this%g)
    allocate(this%g)
    call this%g%init (this%mesh%ncell*ndim)
    do i = 1,this%mesh%ncell
      do n = 1,ndim
        index = (i-1)*ndim+n
        call this%g%add_edge (index,index)
        do m = 1,ndim
          call this%g%add_edge (index, &
              (pack(this%gmesh%cneighbor(:,i), mask=this%gmesh%cneighbor(:,i)>0) - 1) * ndim + m)
        end do
      end do
    end do
    call this%g%add_complete ()

    call this%lhs%init (this%g, take_graph=.false.)

  end subroutine update_matrix

  subroutine predictor_solver_delete (this)
    type(predictor_solver) :: this
    if (associated(this%g)) deallocate(this%g)
  end subroutine predictor_solver_delete
  
  ! note 1: This is a crude way to account for solid material within the cell.
  !         In SOLVE_FOR_VELOCITY we also divide by FluidVof to account for
  !         the mass of fluid in the cell, more or less canceling out this
  !         term. Momentum advection is specifically excluded because the VOF
  !         is already accounting for the solid material.

  ! velocity_cc is inout
  subroutine solve (this, dt, fluidRho_n, fluidVof_n, velocity_cc_n, solid_face, is_pure_immobile)

    use hypre_hybrid_type

    class(predictor_solver), intent(inout) :: this
    real(r8),          intent(in)    :: dt, fluidRho_n(:), fluidVof_n(:), velocity_cc_n(:,:)
    logical,           intent(in)    :: solid_face(:), is_pure_immobile(:)

    real(r8) :: rhs(ndim*this%mesh%ncell), mass(ndim), mass_flux, cell_rhs(ndim), &
        vel_cc_new(ndim*this%mesh%ncell), &
        grad_vel_face_out_bndry(ndim,this%mesh%nface), viscosity_face_bndry(this%mesh%nface)
    integer :: i, i_ngbr, f, stat
    type(hypre_hybrid) :: solver

    call start_timer ('prediction')

    ! initialize the sparse left hand side matrix
    if (.not.this%inviscid .and. this%viscous_implicitness > 0.0_r8) call this%lhs%set_all (0.0_r8)
    
    ! update boundary face values prior to the rest of the domain
    call start_timer ('boundaries')
    grad_vel_face_out_bndry = 0.0_r8
    viscosity_face_bndry = 0.0_r8
    rhs = 0.0_r8
    call apply_velocity_bcs (grad_vel_face_out_bndry, viscosity_face_bndry, rhs, this%velocity_bc, &
        this%velocity_cc, this%vof, this%fluidVof, dt, this%viscous_implicitness, this%inviscid, &
        this%mprop, this%mesh, this%gmesh)
    call apply_pressure_bcs (viscosity_face_bndry, rhs, this%pressure_bc, this%velocity_cc, &
        this%vof, this%fluidVof, dt, this%viscous_implicitness, this%inviscid, this%mprop, &
        this%mesh, this%gmesh)
    call stop_timer ('boundaries')

    call start_timer ('internal')
    do i = 1,this%mesh%ncell
      cell_rhs = 0.0_r8

      ! pressure
      cell_rhs = cell_rhs - dt*this%fluidRho(i)*this%gradP_dynamic_over_rho_cc(:,i)

      ! viscosity
      if (.not.this%inviscid) &
          call apply_viscosity (cell_rhs, this%lhs, this%mprop, dt, this%viscous_implicitness, &
          velocity_cc_n, this%vof, this%fluidVof, this%fluidRho(i), &
          grad_vel_face_out_bndry(:,this%mesh%cface(:,i)), &
          viscosity_face_bndry(this%mesh%cface(:,i)), this%mesh%area(this%mesh%cface(:,i)), &
          this%mesh%volume(i), solid_face(this%mesh%cface(:,i)), is_pure_immobile(i), i, &
          this%mesh%cface(:,i), this%gmesh)

      ! TODO: turbulence
      
      ! account for solid material (see note 1)
      cell_rhs = cell_rhs * this%fluidVof(i)

      ! advect momentum
      do f = 1,nfc
        mass_flux = sum(this%mprop%density*this%volume_flux(:,f,i))
        ! upwind to determine which velocity to grab
        if (mass_flux > 0.0_r8) then ! outflow
          cell_rhs = cell_rhs - mass_flux*velocity_cc_n(:,i) / this%mesh%volume(i)
        else ! inflow
          i_ngbr = this%gmesh%cneighbor(f,i)
          if (i_ngbr > 0) then ! not a boundary
            cell_rhs = cell_rhs - mass_flux*velocity_cc_n(:,i_ngbr) / this%mesh%volume(i)
          else
            ! TODO: get this from boundary condition, including tangential components
            !       Truchas grabs the dirichlet velocity BC when available, and the
            !       cell center velocity where not (zero velocity gradient there)
            cell_rhs = cell_rhs &
                - mass_flux*this%fluxing_velocity(f,i)*this%gmesh%outnorm(:,f,i)/this%mesh%volume(i)
          end if
        end if
      end do
      
      ! TODO: solidifying flow

      ! TODO: surface tension

      ! TODO: porous drag
      
      ! TODO: body forces

      ! previous momentum value
      cell_rhs = cell_rhs + fluidRho_n(i)*fluidVof_n(i)*velocity_cc_n(:,i)

      ! update the velocity explicitly if there is no need for the implicit solve
      ! note that updating velocity_cc inside this loop necessitates the use of
      ! the copy velocity_cc_n above when accessing neighbors
      ! TODO: * expand the "mass matrix" with the other physics that gets added,
      !         possibly warranting breaking this entire section into its own subroutine.
      !       * do the implicit stuff after this loop
      if (this%inviscid .or. this%viscous_implicitness == 0.0_r8) then
        mass = this%fluidRho(i)*this%fluidVof(i)
        this%velocity_cc(:,i) = merge(cell_rhs / mass, 0.0_r8, mask=mass>0.0_r8)
      else
        rhs((i-1)*ndim+1:(i-1)*ndim+ndim) = rhs((i-1)*ndim+1:(i-1)*ndim+ndim) + cell_rhs
      end if
    end do
    call stop_timer ('internal')

    ! solve the implicit system for viscous flows
    call start_timer ('solver')
    if (.not.this%inviscid .and. this%viscous_implicitness > 0.0_r8) then
      ! initial guess is current velocity (TODO: better initial guess?)
      vel_cc_new = reshape(velocity_cc_n, [ndim*this%mesh%ncell])

      call solver%init (this%lhs, this%params)
      call solver%setup ()
      call solver%solve (rhs, vel_cc_new, stat)
      if (stat /= 0) call LS_fatal ("projection solver failed")

      ! copy result into velocity_cc (TODO: can I do all this directly without copying?)
      this%velocity_cc = reshape(vel_cc_new, [ndim, this%mesh%ncell])
    end if
    call stop_timer ('solver')

    call stop_timer ('prediction')

  contains

    subroutine apply_velocity_bcs (grad_vel_face_out_bndry, viscosity_face_bndry, rhs, velocity_bc, &
        velocity_cc, vof, fluidVof, dt, viscous_implicitness, inviscid, mprop, mesh, gmesh)

      use bndry_func_class
      use mesh_geom_type
      use differential_operators, only: faceGradient
      use array_utils,            only: magnitude2

      real(r8),          intent(inout) :: grad_vel_face_out_bndry(:,:), viscosity_face_bndry(:), rhs(:)
      class(bndry_func), intent(in)    :: velocity_bc
      real(r8),          intent(in)    :: velocity_cc(:,:), vof(:,:), fluidVof(:), dt, &
          viscous_implicitness
      logical,           intent(in)    :: inviscid
      type(matl_props),  intent(in)    :: mprop
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
              grad_vel_face_out_bndry(n,fid) = dot_product( &
                  faceGradient ([velocity_cc(n,i), velocity_bc_value*gmesh%outnorm(n,f,i)], &
                  reshape([gmesh%xc(:,i), gmesh%fc(:,fid)], [ndim,2])), &
                  gmesh%outnorm(:,f,i))
            end do
          end if

          ! implicit
          ! TODO: allow tangential components to dirichlet BCs
          if (viscous_implicitness > 0.0_r8) then
            dx = gmesh%xc(:,i) - gmesh%fc(:,fid)
            rhs((i-1)*ndim+1:(i-1)*ndim+ndim) = rhs((i-1)*ndim+1:(i-1)*ndim+ndim) &
                + dt*viscous_implicitness &
                *viscosity_face_bndry(fid)*velocity_bc_value*gmesh%outnorm(:,f,i) &
                *dot_product(gmesh%outnorm(:,f,i), dx) / magnitude2(dx) / mesh%volume(i)
          end if
        end if
      end do

    end subroutine apply_velocity_bcs

    subroutine apply_pressure_bcs (viscosity_face_bndry, rhs, pressure_bc, &
        velocity_cc, vof, fluidVof, dt, viscous_implicitness, inviscid, mprop, mesh, gmesh)

      use bndry_func_class
      use mesh_geom_type
      use array_utils, only: magnitude2

      real(r8),          intent(inout) :: viscosity_face_bndry(:), rhs(:)
      class(bndry_func), intent(in)    :: pressure_bc
      real(r8),          intent(in)    :: velocity_cc(:,:), vof(:,:), fluidVof(:), dt, &
          viscous_implicitness
      logical,           intent(in)    :: inviscid
      type(matl_props),  intent(in)    :: mprop
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
                *mesh%area(fid)*dot_product(gmesh%outnorm(:,f,i), dx) / magnitude2(dx) / mesh%volume(i)
          end if
        end if
      end do

    end subroutine apply_pressure_bcs

  end subroutine solve

  subroutine apply_viscosity (rhs, lhs, mprop, dt, viscous_implicitness, velocity_cc, vof, &
      fluidVof, fluidRho, grad_vel_face_out_bndry, viscosity_face_bndry, face_area, cell_vol, &
      solid_face, is_pure_immobile, i, cface, gmesh)
    
    use csr_matrix_type
    use mesh_geom_type
    use matl_props_type
    use array_utils, only: magnitude2

    real(r8),         intent(inout) :: rhs(:)
    type(csr_matrix), intent(inout) :: lhs
    type(matl_props), intent(in)    :: mprop
    real(r8),         intent(in)    :: dt, viscous_implicitness, velocity_cc(:,:), vof(:,:), &
        fluidVof(:), fluidRho, grad_vel_face_out_bndry(:,:), viscosity_face_bndry(:), &
        face_area(:), cell_vol
    logical,          intent(in)    :: solid_face(:), is_pure_immobile
    integer,          intent(in)    :: i, cface(:)
    type(mesh_geom),  intent(in)    :: gmesh

    real(r8) :: dx(ndim), viscosity_face(nfc), tmp
    integer  :: f,n, index, index_ngbr, i_ngbr

    call start_timer ('viscosity')

    ! calculate face viscosities
    call start_timer ('face viscosity')
    viscosity_face = viscosityFaces (mprop, vof, fluidVof, viscosity_face_bndry, solid_face, i, &
        gmesh%cneighbor(:,i))
    call stop_timer ('face viscosity')

    ! explicit part
    call start_timer ('explicit')
    if (viscous_implicitness < 1.0_r8) &
        rhs = rhs + viscousExplicit (dt, viscous_implicitness, velocity_cc, &
        grad_vel_face_out_bndry, viscosity_face, face_area, cell_vol, i, gmesh)
    call stop_timer ('explicit')

    ! implicit part
    call start_timer ('implicit')
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
                / magnitude2(dx) / cell_vol

            if (i_ngbr > 0) call lhs%increment (index,index_ngbr, tmp)
            call lhs%increment (index,index, -tmp)
          end do
        end do
      end if
    end if
    call stop_timer ('implicit')

    call stop_timer ('viscosity')

  end subroutine apply_viscosity

  pure function viscousExplicit (dt, viscous_implicitness, velocity_cc, &
      grad_vel_face_out_bndry, viscosity_face, face_area, cell_vol, i, gmesh) result(dMomentum)

    use mesh_geom_type

    real(r8),        intent(in) :: dt, viscous_implicitness, velocity_cc(:,:), &
        grad_vel_face_out_bndry(:,:), viscosity_face(:), face_area(:), cell_vol
    integer,         intent(in) :: i
    type(mesh_geom), intent(in) :: gmesh
    real(r8)                    :: dMomentum(ndim)

    dMomentum = dt * (1.0_r8 - viscous_implicitness) * divStress(velocity_cc, &
        grad_vel_face_out_bndry, viscosity_face, face_area, cell_vol, i, gmesh)

  end function viscousExplicit

  pure function divStress (velocity_cc, grad_vel_face_out_bndry, viscosity_face, face_area, &
      cell_vol, i, gmesh)

    use mesh_geom_type
    use differential_operators, only: faceGradient, divergence

    real(r8),        intent(in) :: velocity_cc(:,:), grad_vel_face_out_bndry(:,:), &
        viscosity_face(:), face_area(:), cell_vol
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
          grad_vel_face_out(f) = grad_vel_face_out_bndry(n,f)
        end if
      end do

      ! calculate the divergence of the viscosity * the gradient of the nth velocity component
      divStress(n) = divergence (viscosity_face * grad_vel_face_out, face_area, cell_vol)
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
