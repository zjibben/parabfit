program driver

  !use exo_c_binding
  use unstr_mesh_type
  use unstr_mesh_factory
  use logging_services
  use parameter_list_type
  use parameter_list_json
  use,intrinsic :: iso_fortran_env, only: output_unit, error_unit
  use,intrinsic :: iso_c_binding, only: C_NEW_LINE
  use kinds, only: r8
  implicit none

  !type(exo_info) :: exo
  integer :: n, num_arg, inlun, j
  character(256) :: arg
  character(:), allocatable :: prog, infile, mesh_file, errmsg
  type(unstr_mesh), pointer :: mesh
  type(parameter_list), pointer :: params

  !! Get the program name from the command line.
  call get_command (arg)
  n = scan(arg, '/', back=.true.)
  prog = trim(arg(n+1:))  ! remove the leading path component, if any

  !! Get the input file name from the command line.
  num_arg = command_argument_count()
  if (num_arg == 1) then
    call get_command_argument (1, arg)
    infile = trim(arg)
  else
    write(error_unit,'(a)') 'usage: ' // prog // ' infile'
    stop
  end if

  !! Initialize the logging service routines; output goes to the stdout.
  call LS_initialize ([output_unit], LS_VERB_NOISY)

  !! Read the parameter list from the input file.
  open(newunit=inlun,file=infile,action='read',access='stream')
  call parameter_list_from_json_stream (inlun, params, errmsg)
  if (.not.associated(params)) call LS_fatal ("error reading input file:" // C_NEW_LINE // errmsg)
  close(inlun)

  !! Do something with the mesh
  call params%get ('mesh file', mesh_file)
  mesh => new_unstr_mesh(mesh_file)
  !call mesh%dump

  !call test_func_factories
  !call test_BC
  !call test_graph_module
  !call test_csr_matrix
  !call test_mfd_diff_precon_type
  !call test_HT_model_type
  !call test_face_solve
  !call test_initial_state
  !call test_mesh_plist
  call test_HT_sim

  call LS_exit

contains

  subroutine test_HT_sim
    use HT_sim_type
    type(HT_sim) :: sim
    integer :: stat
    character(:), allocatable :: errmsg
    call sim%init (params)
    call sim%run (stat, errmsg)
    if (stat /= 0) call LS_fatal ('HT_SIM: '//errmsg)
  end subroutine

  subroutine test_mesh_plist
    use unstr_mesh_factory
    use parameter_list_type
    type(parameter_list), pointer :: plist
    type(unstr_mesh), pointer :: my_mesh
    plist => params%sublist('mesh')
    my_mesh => new_unstr_mesh(plist)
  end subroutine

  subroutine test_initial_state

    use mfd_disc_type
    use HT_model_type
    use HT_solver_type
    use unstr_mesh_func
    use unstr_mesh_gmv

    type(parameter_list), pointer :: plist
    type(mfd_disc), target :: disc
    type(HT_model), target :: model
    type(HT_solver) :: solver

    real(r8), allocatable, target :: u(:), udot(:)
    real(r8), pointer :: hcell(:), tcell(:), hdot(:), tdot(:), tface(:)
    real(r8) :: temp(mesh%ncell)

    !! Initialize the discretization
    call disc%init (mesh)

    !! Create the heat transfer model
    plist => params%sublist('HT-model')
    call model%init (disc, plist)

    !! Create the heat transfer solver
    plist => params%sublist('HT-solver')
    call solver%init (model, plist)

    !! Generate the initial temperature field
    plist => params%sublist('initial-temperature')
    call compute_mesh_func (mesh, plist, temp)

    !! Solve for the initial state and its derivative
    allocate(u(model%num_dof()), udot(model%num_dof()))
    call solver%test_initial_state (0.0_r8, temp, 1.0d-3, u, udot)

    call model%get_cell_temp_view (u, tcell)
    call model%get_cell_heat_view (u, hcell)

    call model%get_cell_temp_view (udot, tdot)
    call model%get_cell_heat_view (udot, hdot)

    call gmv_open ('output.gmv')
    call gmv_write_unstr_mesh (mesh)
    call gmv_begin_variables
    call gmv_write_cell_var (mesh, tcell, 'T')
    call gmv_write_cell_var (mesh, hcell, 'H')
    call gmv_write_cell_var (mesh, tdot, 'dT/dt')
    call gmv_write_cell_var (mesh, hdot, 'dH/dt')
    call gmv_end_variables
    call gmv_write_unstr_mesh_surf (mesh)
    call gmv_begin_surfvars
    call model%get_face_temp_view (u, tface)
    call gmv_write_surf_var (mesh, tface, 'T')
    call model%get_face_temp_view (udot, tface)
    call gmv_write_surf_var (mesh, tface, 'dT/dt')
    call gmv_end_surfvars
    call gmv_close

  end subroutine test_initial_state

  subroutine test_face_solve

    use mfd_disc_type
    use HT_model_type
    use HT_AE_solver_type

    type(parameter_list), pointer :: plist
    type(mfd_disc), target :: disc
    type(HT_model), target :: model
    type(HT_AE_solver) :: solver
    real(r8) :: tcell(mesh%ncell), tface(mesh%nface)
    integer :: i, stat
    character(:), allocatable :: errmsg

    !! Initialize the discretization
    call disc%init (mesh)

    plist => params%sublist('HT-model')
    call model%init (disc, plist)

    plist => params%sublist('PCG-params')
    call solver%init (model, plist)

    !tcell = 1.0_r8
    do i = 1, mesh%ncell
      tcell(i) = sum(mesh%x(1,mesh%cnode(:,i))) / size(mesh%cnode,1)
    end do
    tface = 0.0_r8
    call solver%solve (0.0_r8, tcell, tface, stat, errmsg)
    if (stat /= 0) then
      print *, 'ERROR: HT_AE_solver%solve returned stat=', stat, ':', errmsg
    end if
    !print *, tface

    do i = 1, size(mesh%face_set_id)
      print *, 'Boundary surface ', mesh%face_set_id(i), ' temps:'
      print *, pack(tface, mask=btest(mesh%face_set_mask,pos=i))
    end do

  end subroutine

  subroutine test_HT_model_type

    use mfd_disc_type
    use HT_model_type

    type(parameter_list), pointer :: plist
    type(mfd_disc), target :: disc
    type(HT_model) :: model
    integer :: j

    !! Initialize the discretization
    call disc%init (mesh)

    plist => params%sublist('HT-test')
    call model%init (disc, plist)

    call model%temp_bc%compute (0.0_r8)
    associate (faces => model%temp_bc%index, values => model%temp_bc%value)
      print '(i0,":",es12.4)', (faces(j), values(j), j = 1, size(faces))
    end associate

    call model%flux_bc%compute (0.0_r8)
    associate (faces => model%flux_bc%index, values => model%flux_bc%value)
      print '(i0,":",es12.4)', (faces(j), values(j), j = 1, size(faces))
    end associate

  end subroutine

  subroutine test_mfd_diff_precon_type

    use mfd_disc_type
    use mfd_diff_matrix_type
    use mfd_diff_precon_type
    use parameter_list_type

    use bc_factory_type
    use bndry_func_class

    type(mfd_diff_matrix), allocatable :: dm
    type(mfd_diff_matrix), pointer :: dm_ptr
    type(mfd_diff_precon) :: pc
    type(mfd_disc), target :: disc
    type(parameter_list), pointer :: plist
    real(r8) :: coef(mesh%ncell), f1(mesh%ncell), f2(mesh%nface)

    type(bc_factory) :: bcfac
    class(bndry_func), allocatable :: bc_dir, bc_flux

    !! Initialize the discretization
    call disc%init (mesh)

    !! Initialize the diffusion matrix
    allocate(dm)
    call dm%init (disc)

    !! Initialize the preconditioner
    plist => params%sublist ('PC-test-precon')
    call pc%init (dm, plist)

    !! Now lets define the matrix values.
    dm_ptr => pc%matrix()
    coef = 1.0_r8
    call dm_ptr%compute (coef)

    call dm_ptr%incr_cell_diag (coef)

    !! And apply some BC fixups.
    plist => params%sublist ('PC-test-BC')
    call bcfac%init(mesh, plist)
    call bcfac%alloc_bc ('dirichlet', bc_dir)
    call bc_dir%compute(0.0_r8)
    call dm_ptr%set_dir_faces (bc_dir%index)

    !! Now compute the preconditioner
    call pc%compute

    f1 = 1.0_r8
    f2 = 0.0_r8
    f2(bc_dir%index) = bc_dir%value

    call pc%apply (f1, f2)
    print '(/,(8es10.2))', f1
    print '(/,(8es10.2))', f2

  end subroutine

  subroutine test_csr_matrix

    use csr_matrix_type

    type(csr_matrix) :: A
    type(csr_graph), pointer :: g
    integer :: j, k1, k2

    allocate(g)
    call g%init (mesh%nface)
    do j = 1, mesh%ncell
      call g%add_clique (mesh%cface(:,j))
    end do
    call g%add_complete

    call A%init (g, take_graph=.true.)
    call A%set_all (0.0_r8)
    do j = 1, mesh%ncell
      do k2 = 1, size(mesh%cface,dim=1)
        do k1 = 1, size(mesh%cface,dim=1)
          call A%increment (mesh%cface(k1,j), mesh%cface(k2,j), 1.0_r8)
        end do
      end do
    end do

  end subroutine test_csr_matrix

!  subroutine test_graph_module
!
!    use GraphModule
!    use graph_type
!
!    type(NGraphType) :: g1
!    type(graph) :: g2
!    integer, allocatable :: xadj1(:), adjncy1(:)
!    integer, allocatable :: xadj2(:), adjncy2(:)
!    real :: cpu0, cpu
!
!    call cpu_time (cpu0)
!    g1 = CreateGraph(mesh%nface)
!    do j = 1, mesh%ncell
!      call AddClique (g1, mesh%cface(:,j))
!    end do
!    call GetNeighborStructure (g1, xadj1, adjncy1)
!    call cpu_time (cpu)
!    print *, 'Old=', cpu-cpu0
!
!    cpu0 = cpu
!    call g2%init (mesh%nface)
!    do j = 1, mesh%ncell
!      call g2%add_clique (mesh%cface(:,j))
!    end do
!    call g2%get_adjacency (xadj2, adjncy2)
!    call cpu_time (cpu)
!    print *, 'New=', cpu-cpu0
!
!    print *, 'same xadj?', all(xadj1==xadj2)
!    print *, 'same adjncy?', all(adjncy1==adjncy2)
!
!  end subroutine

  subroutine test_func_factories

    use scalar_func_factories
    use parameter_list_type

    type(parameter_list), pointer :: plist
    class(scalar_func), allocatable :: foo

    plist => params%sublist ('f1')
    call alloc_scalar_func (foo, plist)
    print *, foo%eval([0.0_r8])

    plist => params%sublist ('f2')
    call alloc_scalar_func (foo, plist)
    print *, foo%eval([1.5_r8])

    plist => params%sublist ('f3')
    call alloc_scalar_func (foo, plist)
    print *, foo%eval([1.5_r8])

    plist => params%sublist ('f4')
    call alloc_scalar_func (foo, plist)
    print *, foo%eval([0.0_r8])

  end subroutine

  subroutine test_bc

    use bc_factory_type
    use bndry_func_class
    use parameter_list_type

    type(parameter_list), pointer :: plist
    type(bc_factory) :: bcfac
    class(bndry_func), allocatable :: bc_dir, bc_flux

    plist => params%sublist ('BC')
    call bcfac%init(mesh, plist)
    !bcfac = bc_factory(mesh, plist)

    call bcfac%alloc_bc ('dirichlet', bc_dir)
    call bc_dir%compute(0.0_r8)
    print '("[",i0,"] ",i0,es10.2)', (j, bc_dir%index(j), bc_dir%value(j), j=1, size(bc_dir%index))

    call bcfac%alloc_bc ('flux', bc_flux)
    call bc_flux%compute(0.5_r8)
    print '("[",i0,"] ",i0,es10.2)', (j, bc_flux%index(j), bc_flux%value(j), j=1, size(bc_flux%index))

    call bc_flux%compute(1.5_r8)
    print '("[",i0,"] ",i0,es10.2)', (j, bc_flux%index(j), bc_flux%value(j), j=1, size(bc_flux%index))

  end subroutine

end program
