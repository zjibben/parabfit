!
!   Define procedures cell-face and cell-centered velocity 
!   projection operations.
!

module projection_module

  use kinds,  only: r8
  use consts, only: ndim,nfc
  use csr_matrix_type
  use unstr_mesh_type
  use mesh_geom_type
  use parameter_list_type
  use bndry_func_class
  use fischer_guess_type
  use logging_services
  use timer_tree_type
  implicit none
  private

  ! TODO: make these a user-set parameters
  real(r8), parameter :: void_pressure = 0.0_r8
  real(r8), parameter :: min_face_fraction = 1e-3_r8
  
  type, public :: projection_solver
    private

    type(unstr_mesh),     pointer :: mesh => null()
    type(mesh_geom),      pointer :: gmesh => null()
    type(parameter_list), pointer :: params

    type(csr_graph), pointer :: g => null()
    type(csr_matrix), pointer :: lhs => null()
    type(fischer_guess) :: fischer
  contains
    procedure :: init => init_projection_solver
    procedure :: solve
    procedure :: update_matrix
    procedure, private :: pressure_poisson
    procedure, private :: velocity_projection
    final :: projection_solver_delete
  end type projection_solver

contains

  ! the projection step involves two calculations:
  !   1. solving the pressure Poisson system for the change in pressure
  !   2. updating the cell center and fluxing velocities
  subroutine solve (this, velocity_cc, pressure_cc, gradP_dynamic_over_rho_cc, fluxing_velocity, &
      dt, fluidRho, fluidRho_n, min_fluidRho, fluidVof, vof, body_force, solid_face, &
      is_pure_immobile, mprop, velocity_bc, pressure_bc)

    use matl_props_type
    use bndry_func_class

    class(projection_solver), intent(inout) :: this
    real(r8), intent(inout) :: velocity_cc(:,:), pressure_cc(:), gradP_dynamic_over_rho_cc(:,:)
    real(r8), intent(out) :: fluxing_velocity(:,:)
    real(r8), intent(in) :: dt, fluidRho(:), fluidRho_n(:), min_fluidRho, fluidVof(:), vof(:,:), &
        body_force(:)
    logical, intent(in) :: solid_face(:), is_pure_immobile(:)
    type(matl_props), intent(in) :: mprop
    class(bndry_func), intent(in) :: velocity_bc, pressure_bc
    
    real(r8) :: dP(this%mesh%ncell), rho_face(nfc,this%mesh%ncell)

    call start_timer ('projection')

    call this%pressure_poisson (pressure_cc, dP, rho_face, fluxing_velocity, dt, velocity_cc, &
      gradP_dynamic_over_rho_cc, fluidRho, fluidRho_n, min_fluidRho, fluidVof, vof, body_force, &
      solid_face, is_pure_immobile, mprop, velocity_bc, pressure_bc)
    call this%velocity_projection (velocity_cc, pressure_cc, fluxing_velocity, &
      gradP_dynamic_over_rho_cc, rho_face, dP, dt, fluidRho, body_force, solid_face, &
      is_pure_immobile, velocity_bc, pressure_bc)

    call stop_timer ('projection')

  end subroutine solve

  ! set up and solve the pressure Poisson system
  ! here, pressure_cc is inout, and fluxing_velocity is out
  subroutine pressure_poisson (this, pressure_cc, dP, rho_face, fluxing_velocity, dt, velocity_cc, &
      gradP_dynamic_over_rho_cc, fluidRho, fluidRho_n, min_fluidRho, fluidVof, vof, body_force, &
      solid_face, is_pure_immobile, mprop, velocity_bc, pressure_bc)

    use hypre_hybrid_type
    use matl_props_type
    use array_utils, only: meanArithmetic

    class(projection_solver), intent(inout) :: this
    real(r8), intent(inout) :: pressure_cc(:)
    real(r8), intent(out)   :: dP(:), rho_face(:,:), fluxing_velocity(:,:)
    real(r8), intent(in)    :: dt, velocity_cc(:,:), gradP_dynamic_over_rho_cc(:,:), &
        fluidRho(:), fluidRho_n(:), min_fluidRho, fluidVof(:), vof(:,:), body_force(:)
    logical,  intent(in)    :: solid_face(:), is_pure_immobile(:)
    type(matl_props), intent(in) :: mprop
    class(bndry_func), intent(in) :: velocity_bc, pressure_bc

    real(r8) :: rhs(this%mesh%ncell)
    integer :: i,f,i_ngbr,fid

    call start_timer ('pressure poisson')

    ! initialize the sparse left hand side matrix
    call this%lhs%set_all (0.0_r8)

    ! update boundary face values prior to the rest of the domain
    call start_timer ('boundaries')
    call apply_velocity_bcs (fluxing_velocity, rho_face, rhs, velocity_bc, pressure_cc, fluidRho, &
        min_face_fraction, min_fluidRho, dt, solid_face, is_pure_immobile, this%mesh, this%gmesh)
    call apply_pressure_bcs (fluxing_velocity, rho_face, rhs, pressure_bc, velocity_cc, &
        pressure_cc, gradP_dynamic_over_rho_cc, fluidRho, body_force, min_face_fraction, &
        min_fluidRho, dt, solid_face, is_pure_immobile, this%mesh, this%gmesh)
    call stop_timer ('boundaries')

    ! WARNING: what about boundaries without a boundary condition set?
    
    ! update the fluxing velocity, and set the values for the lhs and rhs of the Poisson system
    call start_timer ('internal')
    do i = 1,this%mesh%ncell
      ! calculate quantities on this cell's faces
      do f = 1,nfc
        i_ngbr = this%gmesh%cneighbor(f,i) ! neighbor index
        fid = this%mesh%cface(f,i) ! face id
        if (i_ngbr < 1) cycle ! skip boundary faces

        ! calculate face quantities
        rho_face(f,i) = fluid_density_face (this%mesh%volume([i, i_ngbr]), fluidVof([i, i_ngbr]), &
            fluidRho([i, i_ngbr]), min_face_fraction, min_fluidRho, solid_face(fid))

        fluxing_velocity(f,i) = calc_fluxing_velocity (dt, rho_face(f,i), fluidRho([i, i_ngbr]), &
            fluidRho_n([i, i_ngbr]), velocity_cc(:,[i, i_ngbr]), pressure_cc([i, i_ngbr]), &
            gradP_dynamic_over_rho_cc(:,[i, i_ngbr]), body_force, solid_face(fid), &
            this%gmesh%outnorm(:,f,i), this%gmesh%xc(:,[i, i_ngbr]), this%gmesh%fc(:,fid))
      end do
      
      ! calculate the base Poisson system left and right hand sides
      rhs(i) = cell_rhs (fluxing_velocity(:,i), this%mesh%area(this%mesh%cface(:,i)), &
          this%mesh%volume(i), fluidRho(i), pressure_cc(i), this%gmesh%outnorm(:,:,i), &
          min_face_fraction, min_fluidRho, dt, solid_face(this%mesh%cface(:,i)), is_pure_immobile(i))
      
      call cell_lhs (this%lhs, this%mesh%volume(i), dt, rho_face(:,i), fluidRho(i), vof(:,i), &
          this%gmesh%fc(:,this%mesh%cface(:,i)), i, is_pure_immobile(i), &
          solid_face(this%mesh%cface(:,i)), mprop, this%mesh, this%gmesh)

      ! TODO: add different physics (e.g., surface tension)
    end do
    call stop_timer ('internal')

    call linear_solve (dP, this%lhs, rhs, this%fischer, this%params)

    ! re-normalize pressure when not using dirichlet pressure BCs
    if ((.not.allocated(pressure_bc%index) .or. size(pressure_bc%index)==0) &
        .and. .not.any(fluidRho==0.0_r8 .or. is_pure_immobile)) &
        dP = dP - meanArithmetic(dP)

    ! increment the cell centered pressure
    pressure_cc = pressure_cc + dP

    call stop_timer ('pressure poisson')
    
  contains

    subroutine linear_solve (dP, lhs, rhs, fischer, params)

      use consts, only: alittle

      real(r8),             intent(out)   :: dP(:)
      type(csr_matrix),     intent(in), target :: lhs
      real(r8),             intent(in)    :: rhs(:)
      type(fischer_guess),  intent(inout) :: fischer
      type(parameter_list), pointer       :: params

      type(hypre_hybrid) :: solver
      integer :: stat

      call start_timer ('solver')

      call solver%init (lhs, params)
      call solver%setup ()

      ! initial guess
      if (sum(abs(rhs)) > alittle) then
        call fischer%guess (dP, rhs)
      else
        dP = 1.0_r8
      end if

      call solver%solve (rhs, dP, stat)
      if (stat /= 0) call LS_fatal ("projection solver failed")

      call fischer%update (dP, rhs, this%lhs)
      
      call stop_timer ('solver')

    end subroutine linear_solve

    subroutine apply_velocity_bcs (fluxing_velocity, rho_face, rhs, velocity_bc, &
        pressure_cc, fluidRho, min_face_fraction, min_fluidRho, dt, solid_face, is_pure_immobile, &
        mesh, gmesh)

      use bndry_func_class

      real(r8),          intent(inout) :: fluxing_velocity(:,:), rho_face(:,:), rhs(:)
      class(bndry_func), intent(in)    :: velocity_bc
      real(r8),          intent(in)    :: pressure_cc(:), fluidRho(:), min_face_fraction, &
          min_fluidRho, dt
      logical,           intent(in)    :: solid_face(:), is_pure_immobile(:)
      type(unstr_mesh),  intent(in)    :: mesh
      type(mesh_geom),   intent(in)    :: gmesh

      real(r8) :: fluxing_velocity_bc
      integer  :: bndry_f, fid, i, f

      ! update values on boundary faces
      do bndry_f = 1,size(velocity_bc%index)
        fid = velocity_bc%index(bndry_f)
        
        i = gmesh%fcell(1,fid) ! id of cell attached to this face
        f = gmesh%flid(1,fid)  ! local id of the face
        
        ! calculate face boundary quantities we missed earlier
        rho_face(f,i) = merge(max(min_face_fraction*min_fluidRho, fluidRho(i)), 0.0_r8, &
            fluidRho(i) > 0.0_r8)
        fluxing_velocity(f,i) = merge(velocity_bc%value(bndry_f), 0.0_r8, .not.solid_face(fid))
      end do

    end subroutine apply_velocity_bcs

    subroutine apply_pressure_bcs (fluxing_velocity, rho_face, rhs, pressure_bc, &
        velocity_cc, pressure_cc, gradP_dynamic_over_rho_cc, fluidRho, body_force, &
        min_face_fraction, min_fluidRho, dt, solid_face, is_pure_immobile, mesh, gmesh)

      use body_force_module, only: calc_gravityhead_bndry
      use bndry_func_class
      use differential_operators, only: faceGradient

      real(r8),          intent(inout) :: fluxing_velocity(:,:), rho_face(:,:), rhs(:)
      class(bndry_func), intent(in)    :: pressure_bc
      real(r8),          intent(in)    :: velocity_cc(:,:), pressure_cc(:), &
          gradP_dynamic_over_rho_cc(:,:), fluidRho(:), body_force(:), min_face_fraction, &
          min_fluidRho, dt
      logical,           intent(in)    :: solid_face(:), is_pure_immobile(:)
      type(unstr_mesh),  intent(in)    :: mesh
      type(mesh_geom),   intent(in)    :: gmesh

      real(r8) :: pressure_bc_value, gradP_dynamic_face(ndim) !, face_area_over_rho
      integer  :: bndry_f, fid, i, f

      ! update values on boundary faces
      do bndry_f = 1,size(pressure_bc%index)
        fid = pressure_bc%index(bndry_f)
        pressure_bc_value = pressure_bc%value(bndry_f)

        i = gmesh%fcell(1,fid) ! id of cell attached to this face
        f = gmesh%flid(1,fid)  ! local id of the face
        
        ! rho face
        rho_face(f,i) = merge(max(min_face_fraction*min_fluidRho, fluidRho(i)), 0.0_r8, &
            fluidRho(i) > 0.0_r8)

        ! fluxing velocity
        if (solid_face(fid)) then
          fluxing_velocity(f,i) = 0.0_r8
        else
          fluxing_velocity(f,i) = dot_product(gmesh%outnorm(:,f,i), &
              velocity_cc(:,i) + dt*gradP_dynamic_over_rho_cc(:,i))
          if (rho_face(f,i) /= 0.0_r8) then
            gradP_dynamic_face = faceGradient (&
                [pressure_cc(i), pressure_bc_value] + calc_gravityhead_bndry (fluidRho(i), &
                velocity_cc(:,i), gmesh%xc(:,i), gmesh%fc(:,fid), body_force), &
                ! TODO: confirm this is correct/make it prettier
                reshape([gmesh%xc(:,i), gmesh%fc(:,fid)], [ndim,2])) 

            fluxing_velocity(f,i) = fluxing_velocity(f,i) &
                - dt/rho_face(f,i) * dot_product(gradP_dynamic_face, gmesh%outnorm(:,f,i))
          end if
        end if
      end do

      ! ! calculate the base Poisson system left and right hand sides
      !   rhs(i) = cell_rhs (fluxing_velocity(:,i), mesh%area(mesh%cface(:,i)), mesh%volume(i), &
      !       fluidRho(i), pressure_cc(i), gmesh%outnorm(:,:,i), min_face_fraction, min_fluidRho, &
      !       dt, solid_face(mesh%cface(:,i)), is_pure_immobile(i)) !&
      !   !+ pressure_bc_value*face_area_over_rho / mesh%volume(i)
      !   ! NOTE: looking closer at Truchas, it seems like this line is executed, but with
      !   !       pressure_bc_value set to zero, because we need to enforce zero pressure
      !   !       difference, not some the pressure value in this case. I'm not sure why they evaluate
      !   !       this line at all. Also note that in general the pressure BC could change in time,
      !   !       so it isn't necessarily zero, but I don't think Truchas accounts for this.
      !   ! TODO: modify this for time-varying pressure boundary conditions.

    end subroutine apply_pressure_bcs

    real(r8) function cell_rhs (fluxing_velocity, face_area, cell_vol, fluidRho, pressure, &
        outnorm, min_face_fraction, min_fluidRho, dt, solid_face, is_pure_immobile)

      use differential_operators, only: divergence

      real(r8), intent(in) :: fluxing_velocity(:), face_area(:), cell_vol, fluidRho, pressure, &
          outnorm(:,:), min_face_fraction, min_fluidRho, dt
      logical,  intent(in) :: solid_face(:), is_pure_immobile

      ! Compute the RHS for the MAC projection, which is the
      ! cell-centered divergence of the face fluxing velocity.
      ! Compute the cell-centered divergence of the fluxing velocity as a
      ! discrete sum over faces (Green-Gauss).  Recall that the fluxing
      ! velocity has already been projected along the face normal. This will
      ! be the RHS for MAC projection.
      if (FluidRho == 0.0_r8 .or. is_pure_immobile) then
        ! set RHS to zero in void cells because we are solving for the change in pressure
        ! (RHS divided by dt in NONORTHO_PROJECTION)
        ! Scaling of the RHS  requires multiplying by dt here for global scaling of RHS by dt/Volume.
        cell_rhs = void_pressure - pressure &
            / (min_face_fraction * min_fluidRho * cell_vol**(2.0_r8/3.0_r8))
      else
        ! the RHS must be divided by 'dt' so that we're solving for the pressure increment,
        ! and not pressure*dt. Scaling of the RHS here by dt and the cell volume is in
        ! conjunction with the scaled convergence criteria below.

        ! note the fluxing velocity is set to zero on solid faces, so no need to mask by solid_face
        cell_rhs = divergence (fluxing_velocity, face_area, cell_vol) / dt
      end if

      ! we are multiplying the entire Poisson equation by -1,
      ! since the -laplacian operator is positive definite,
      ! allowing the use of CG
      cell_rhs = -cell_rhs

    end function cell_rhs

    ! this calculates a row of the left hand matrix of the pressure poisson system
    ! note: * we are multiplying the entire Poisson equation by -1,
    !         since the -laplacian operator is positive definite,
    !         allowing the use of CG
    subroutine cell_lhs (lhs, cell_vol, dt, rho_face, fluidRho, vof, xf, i, is_pure_immobile, &
        solid_face, mprop, mesh, gmesh)

      use array_utils, only: magnitude2

      type(csr_matrix), intent(inout) :: lhs
      real(r8),         intent(in)    :: cell_vol, dt, fluidRho, vof(:), rho_face(:), xf(:,:)
      integer,          intent(in)    :: i ! local cell id (matrix row)
      logical,          intent(in)    :: is_pure_immobile, solid_face(:)
      type(matl_props), intent(in)    :: mprop
      type(unstr_mesh), intent(in)    :: mesh
      type(mesh_geom),  intent(in)    :: gmesh

      integer  :: f,i_ngbr
      real(r8) :: dx(3), tmp, face_area_over_rho

      ! TODO: set row to zero here rather than outside the threaded loop prior to this subroutine

      ! void and solid cells
      if (fluidRho == 0.0_r8 .or. is_pure_immobile) then
        call lhs%increment (i,i, -1.0_r8/(min_face_fraction*min_fluidRho*cell_vol**(2.0_r8/3.0_r8)))
        return
      end if
      
      ! divergence of pressure gradient term
      ! TODO: modify for nonorthogonal meshes
      do f = 1,nfc
        i_ngbr = gmesh%cneighbor(f,i)
        !if (solid_face(f)) cycle
        ! if (i_ngbr < 1) cycle ! any boundary conditions are applied to the rhs
        
        ! dx  = gmesh%xc(:,i) - gmesh%xc(:,i_ngbr)
        
        if (i_ngbr > 0) then
          dx = gmesh%xc(:,i) - gmesh%xc(:,i_ngbr)
        else
          dx = gmesh%xc(:,i) - gmesh%fc(:,mesh%cface(f,i))
        end if

        face_area_over_rho = merge(0.0_r8, mesh%area(mesh%cface(f,i)) / rho_face(f), &
            rho_face(f) == 0.0_r8 .or. fluidRho == 0.0_r8 .or. is_pure_immobile .or. solid_face(f))
        tmp = face_area_over_rho / cell_vol * dot_product(gmesh%outnorm(:,f,i), dx) / magnitude2(dx)
        
        if (i_ngbr > 0) call lhs%increment (i,i_ngbr, tmp)
        call lhs%increment (i,i, -tmp)
      end do
      
      ! void collapse term
      call lhs%increment (i,i, -cell_compressibility(vof, fluidRho, mprop) / dt**2)

    end subroutine cell_lhs

    ! calculates the cell compressibility for the void collapse model
    real(r8) function cell_compressibility (vof, fluidRho, mprop)

      real(r8),         intent(in) :: vof(:), fluidRho
      type(matl_props), intent(in) :: mprop

      real(r8) :: tmp
      integer, parameter :: navg = 1 ! power for averaging of sound speed by material

      cell_compressibility = - sum(vof**navg / mprop%sound_speed**2, mask=mprop%sound_speed > 0.0_r8)

      ! weight by the vofs in the cell
      tmp = sum(vof**navg, mask=.not.mprop%is_immobile)
      if (tmp      > 0.0_r8) cell_compressibility = cell_compressibility / tmp

      if (fluidRho > 0.0_r8) cell_compressibility = cell_compressibility / fluidRho

    end function cell_compressibility

  end subroutine pressure_poisson

  subroutine init_projection_solver (this, mesh, gmesh, params)

    use parameter_list_type ! intel seems to need this despite it being at the top. might be a bug.

    class(projection_solver), intent(out)        :: this
    type(unstr_mesh),         intent(in), target :: mesh
    type(mesh_geom),          intent(in), target :: gmesh
    type(parameter_list),     intent(in), target, optional :: params

    this%mesh => mesh
    this%gmesh => gmesh
    if (present(params)) this%params => params
    
    call this%update_matrix ()

    call this%fischer%init (this%mesh%ncell)

  end subroutine init_projection_solver

  subroutine update_matrix (this)

    class(projection_solver), intent(inout) :: this

    integer :: i

    if (associated(this%g)) deallocate(this%g)
    if (associated(this%lhs)) deallocate(this%lhs)
    allocate(this%g, this%lhs)
    call this%g%init (this%mesh%ncell)
    do i = 1,this%mesh%ncell
      call this%g%add_edge (i,i)
      call this%g%add_edge (i, pack(this%gmesh%cneighbor(:,i), mask=this%gmesh%cneighbor(:,i)>0))
    end do
    call this%g%add_complete ()

    call this%lhs%init (this%g, take_graph=.false.)

  end subroutine update_matrix

  subroutine projection_solver_delete (this)
    type(projection_solver) :: this
    if (associated(this%g)) deallocate(this%g)
    if (associated(this%lhs)) deallocate(this%lhs)
  end subroutine projection_solver_delete
  
  ! update the cell centered and fluxing velocities
  subroutine velocity_projection (this, velocity_cc, pressure_cc, fluxing_velocity, &
      gradP_dynamic_over_rho_cc, rho_face, dP, dt, fluidRho, body_force, solid_face, &
      is_pure_immobile, velocity_bc, pressure_bc)

    use body_force_module,      only: calc_gravityhead
    use differential_operators, only: faceGradient

    class(projection_solver), intent(inout) :: this
    real(r8), intent(inout) :: velocity_cc(:,:), pressure_cc(:), fluxing_velocity(:,:), &
        gradP_dynamic_over_rho_cc(:,:), rho_face(:,:)
    real(r8), intent(in) :: dP(:), dt, fluidRho(:), body_force(:)
    logical,  intent(in) :: solid_face(:), is_pure_immobile(:)
    class(bndry_func), intent(in) :: velocity_bc, pressure_bc
    
    real(r8) :: grad_dP(ndim), gradP_dynamic_over_rho_cc_n(ndim), gradP_face(ndim,nfc,this%mesh%ncell)
    integer  :: i,f,nid,nfid

    call start_timer ('velocity projection')

    ! apply boundary conditions
    call apply_velocity_bcs (gradP_dynamic_over_rho_cc, gradP_face, rho_face, velocity_bc, &
        velocity_cc, fluidRho, body_force, dt, is_pure_immobile, this%gmesh)
    call apply_pressure_bcs (fluxing_velocity, gradP_face, pressure_bc, velocity_cc, &
        pressure_cc, dP, rho_face, fluidRho, body_force, dt, is_pure_immobile, &
        this%gmesh)

    ! update interior face values and cell centered values
    do i = 1,this%mesh%ncell
      do f = 1,nfc
        nid = this%gmesh%cneighbor(f,i) ! neighbor cell id
        if (nid < 1) cycle ! skip boundary faces
        
        ! correct the fluxing velocities with grad(dP)
        if (rho_face(f,i) > 0.0_r8 .and. .not.solid_face(this%mesh%cface(f,i))) then
          grad_dP = faceGradient (dP([i, nid]), this%gmesh%xc(:,[i, nid]))
          fluxing_velocity(f,i) = fluxing_velocity(f,i) &
              - dt/rho_face(f,i) * dot_product(grad_dP, this%gmesh%outnorm(:,f,i))
        end if

        ! calculate the pressure (not dynamic) gradient on the face
        gradP_face(:,f,i) = faceGradient (pressure_cc([i, nid]) &
            + calc_gravityhead (fluidRho([i, nid]), velocity_cc(:,[i, nid]), &
            this%gmesh%xc(:,[i, nid]), this%gmesh%fc(:,this%mesh%cface(f,i)), body_force), &
            this%gmesh%xc(:,[i, nid]))
      end do

      ! calculate new pressure gradient
      gradP_dynamic_over_rho_cc_n = gradP_dynamic_over_rho_cc(:,i)
      gradP_dynamic_over_rho_cc(:,i) = calc_gradP_dynamic_over_rho_cc (gradP_face(:,:,i), &
          rho_face(:,i), fluidRho(i), this%mesh%area(this%mesh%cface(:,i)), &
          this%gmesh%outnorm(:,:,i), solid_face(this%mesh%cface(:,i)))

      ! correct cell centered velocities
      if (fluidRho(i) == 0.0_r8 .or. is_pure_immobile(i)) then
        velocity_cc(:,i) = 0.0_r8
      else
        velocity_cc(:,i) = velocity_cc(:,i) &
            - dt*(gradP_dynamic_over_rho_cc(:,i) - gradP_dynamic_over_rho_cc_n)
      end if
    end do

    ! Void cells don't have physical (nor solenoidal) Fluxing_Velocity's
    ! Set Fluxing_Velocity on faces that adjoin 'real' fluid cells to the 
    ! corresponding velocity on the other side of the face (corrected for direction)
    do i = 1,this%mesh%ncell
      do f = 1,nfc
        nid = this%gmesh%cneighbor(f,i) ! neighbor cell id
        if (fluidRho(i)==0.0_r8 .and. nid > 0) then
          nfid = this%gmesh%fneighbor(f,i) ! neighbor's local id for face f
          if (fluidRho(nid) > 0.0_r8) &
              fluxing_velocity(f,i) = -fluxing_velocity(nfid,nid)
        end if
      end do
    end do

    call stop_timer ('velocity projection')

  contains
    
    subroutine apply_pressure_bcs (fluxing_velocity, gradP_face, pressure_bc, velocity_cc, &
        pressure_cc, dP, rho_face, fluidRho, body_force, dt, is_pure_immobile, gmesh)

      use body_force_module,      only: calc_gravityhead_bndry
      use differential_operators, only: faceGradient

      real(r8),          intent(inout) :: fluxing_velocity(:,:), gradP_face(:,:,:)
      class(bndry_func), intent(in)    :: pressure_bc
      real(r8),          intent(in)    :: velocity_cc(:,:), pressure_cc(:), dP(:), rho_face(:,:), &
          fluidRho(:), body_force(:), dt
      logical,           intent(in)    :: is_pure_immobile(:)
      type(mesh_geom),   intent(in)    :: gmesh

      real(r8) :: pressure_bc_value
      integer  :: bndry_f, fid, i, f, if

      ! update values on boundary faces
      do bndry_f = 1,size(pressure_bc%index)
        fid = pressure_bc%index(bndry_f)
        pressure_bc_value = pressure_bc%value(bndry_f)

        i = gmesh%fcell(1,fid) ! id of cell attached to this face
        f = gmesh%flid(1,fid)  ! local id of the face

        ! calculate the dynamic pressure gradient for this cell
        gradP_face(:,f,i) = faceGradient ([pressure_cc(i),pressure_bc_value] &
            + calc_gravityhead_bndry (fluidRho(i), velocity_cc(:,i), gmesh%xc(:,i), &
            gmesh%fc(:,fid), body_force), &
            reshape([gmesh%xc(:,i),gmesh%fc(:,fid)], [ndim,2]))

        ! update the fluxing velocity
        if (rho_face(f,i) > 0.0_r8 .and. .not.solid_face(fid)) &
            fluxing_velocity(f,i) = fluxing_velocity(f,i) - dt/rho_face(f,i) * dot_product(&
            faceGradient ([dP(i), 0.0_r8], reshape([gmesh%xc(:,i),gmesh%fc(:,fid)], [ndim,2])), &
            gmesh%outnorm(:,f,i))
      end do

    end subroutine apply_pressure_bcs

    subroutine apply_velocity_bcs (gradP_dynamic_over_rho_cc, gradP_face, rho_face, velocity_bc, &
        velocity_cc, fluidRho, body_force, dt, is_pure_immobile, gmesh)

      use body_force_module,      only: calc_gravityhead_bndry
      use differential_operators, only: faceGradient

      real(r8),          intent(inout) :: gradP_dynamic_over_rho_cc(:,:), gradP_face(:,:,:), &
          rho_face(:,:)
      class(bndry_func), intent(in)    :: velocity_bc
      real(r8),          intent(in)    :: velocity_cc(:,:), fluidRho(:), body_force(:), dt
      logical,           intent(in)    :: is_pure_immobile(:)
      type(mesh_geom),   intent(in)    :: gmesh

      real(r8) :: velocity_bc_value, gradP_dynamic_over_rho_cc_n(ndim)
      integer  :: bndry_f, fid, i, f, if

      do bndry_f = 1,size(velocity_bc%index)
        fid = velocity_bc%index(bndry_f)
        velocity_bc_value = velocity_bc%value(bndry_f) ! TODO: remove this?

        i = gmesh%fcell(1,fid) ! id of cell attached to this face
        f = gmesh%flid(1,fid)  ! local id of the face

        ! pressure gradient on this face is forced to zero ! 0.0_r8
        gradP_face(:,f,i) = 0.0_r8 ! faceGradient (calc_gravityhead_bndry (fluidRho(i), velocity_cc(:,i), &
            ! gmesh%xc(:,i), gmesh%fc(:,fid), body_force), &
        ! reshape([gmesh%xc(:,i),gmesh%fc(:,fid)], [ndim,2]))

        ! this is just a cheap, tricky way to deactivate this face's contribution to the
        ! cell-centered dynamic pressure gradient in calc_gradP_dynamic_over_rho_cc
        rho_face(f,i) = 0.0_r8
      end do

    end subroutine apply_velocity_bcs
    
    ! the cell centered dynamic pressure gradient is calculated as
    ! the average of face centered dynamic pressure gradients,
    ! weighted by the face area and orientation
    ! 
    ! Note: * Truchas's original implementation skips solid and boundary faces.
    !         Since those are already set to zero by the face gradient routine,
    !         there is no need to explicitly skip them again here. However,
    !         we still need to weight by only the cells which are allowed.
    function calc_gradP_dynamic_over_rho_cc (gradP_face, rho_face, fluidRho, face_area, outnorm, &
        solid_face) result(gradP_d)

      use array_utils, only: meanArithmetic

      real(r8), intent(in) :: gradP_face(:,:), rho_face(:), fluidRho, face_area(:), outnorm(:,:)
      logical,  intent(in) :: solid_face(:)
      real(r8)             :: gradP_d(ndim)

      integer :: n

      gradP_d = 0.0_r8
      if (fluidRho /= 0.0_r8) then
        do n = 1,ndim
          gradP_d(n) = meanArithmetic(&
              pack(gradP_face(n,:)/rho_face,    mask=rho_face>0.0_r8 .and. .not.solid_face), &
              pack(abs(face_area*outnorm(n,:)), mask=rho_face>0.0_r8 .and. .not.solid_face))
        end do
      end if

    end function calc_gradP_dynamic_over_rho_cc

  end subroutine velocity_projection


  ! calculate the face-centered density, given cell and flow properties from each neighboring cell
  function fluid_density_face (cell_vol, fluidVof, fluidRho, min_face_fraction, min_fluidRho, &
      solid_face) result(rho_face)

    use array_utils, only: meanArithmetic

    real(r8), intent(in) :: cell_vol(:), fluidVof(:), fluidRho(:), min_face_fraction, min_fluidRho
    logical,  intent(in) :: solid_face
    real(r8)             :: rho_face
    
    ! the face density is the average of the cell centered densities,
    ! weighted by the volume of fluid in each cell
    rho_face = meanArithmetic(fluidRho, cell_vol*fluidVof)

    if (rho_face > 0.0_r8) rho_face = max(min_face_fraction*min_fluidRho, rho_face)

  end function fluid_density_face

  real(r8) function calc_fluxing_velocity (dt, rho_face, fluidRho, fluidRho_n, velocity_cc, &
      pressure_cc, gradP_dynamic_over_rho_cc, body_force, solid_face, outnorm, xc, xf)

    use body_force_module,      only: calc_gravityhead
    use differential_operators, only: faceGradient

    real(r8), intent(in) :: dt, rho_face, velocity_cc(:,:), pressure_cc(:), &
        fluidRho(:), fluidRho_n(:), gradP_dynamic_over_rho_cc(:,:), outnorm(:), xc(:,:), xf(:), &
        body_force(:)
    logical,  intent(in) :: solid_face

    real(r8) :: gradP_dynamic_face(ndim), gh(2)

    if (solid_face) then
      ! there may be sufficient checks in the later lines to remove this if-check altogether
      calc_fluxing_velocity = 0.0_r8
    else
      ! interpolate augmented velocities
      calc_fluxing_velocity = interpolate_velocity_to_faces (velocity_cc, fluidRho_n, &
          gradP_dynamic_over_rho_cc, outnorm, xc, xf, dt)

      if (.not.solid_face .and. rho_face /= 0.0_r8) then
        gh = calc_gravityhead (fluidRho, velocity_cc, xc, xf, body_force)
        gradP_dynamic_face = faceGradient (pressure_cc + gh, xc)
        calc_fluxing_velocity = calc_fluxing_velocity &
            - dt/rho_face * dot_product(gradP_dynamic_face, outnorm)
      end if
    end if

  end function calc_fluxing_velocity

  ! This function interpolates the augmented cell centered velocities to the faces
  ! I say the velocities are augmented because they are really the velocity + dt/rho * gradP_d
  function interpolate_velocity_to_faces (velocity_cc, fluidRho_n, gradP_dynamic_over_rho_cc, &
      outnorm, xc, xf, dt) result(interp_vel)

    use array_utils, only: magnitude2

    real(r8), intent(in) :: velocity_cc(:,:), fluidRho_n(:), gradP_dynamic_over_rho_cc(:,:), &
        outnorm(:), xc(:,:), xf(:), dt
    real(r8)             :: interp_vel

    real(r8)             :: interp_factor, aug_vel_cc(ndim,2)

    aug_vel_cc(:,1) = velocity_cc(:,1) + dt*gradP_dynamic_over_rho_cc(:,1)
    aug_vel_cc(:,2) = velocity_cc(:,2) + dt*gradP_dynamic_over_rho_cc(:,2)

    if (fluidRho_n(2)==0.0_r8) then
      ! WARNING: is this right??
      ! use cell center value if adjacent cell is void
      interp_vel = dot_product(outnorm, aug_vel_cc(:,1))
    else
      ! TODO: note the interpolation factor here is constant with the mesh, so it might be worth
      ! calculating these elsewhere (over in mesh_geom_type?) and saving them.
      interp_factor = dot_product(xf-xc(:,1), xc(:,2)-xc(:,1))/magnitude2(xc(:,2)-xc(:,1))
      interp_vel = dot_product(outnorm, &
          (1.0_r8 - interp_factor)*aug_vel_cc(:,1) + interp_factor*aug_vel_cc(:,2))
    end if

    ! WARNING: need to make sure void cell face velocities match the adjacent non-void cell

  end function interpolate_velocity_to_faces
  
end module projection_module
