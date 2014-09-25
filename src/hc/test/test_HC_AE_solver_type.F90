program test_HC_AE_solver_type

  use kinds, only: r8
  use unstr_mesh_type
  use unstr_mesh_factory
  use unstr_mesh_gmv
  use mfd_disc_type
  use HC_model_type
  use HC_AE_solver_type
  use parameter_list_type
  use logging_services
  use,intrinsic :: iso_fortran_env, only: output_unit
#ifdef NAGFOR
  use,intrinsic :: f90_unix, only: exit
#endif
  implicit none

  integer :: status = 0

  call LS_initialize ([output_unit], LS_VERB_NOISY)

  call test1
  call test2
  call exit (status)

contains

  !! Regular grid.  Should recover test function to machine epsilon, but
  !! only a single CG iteration.

  subroutine test1

    type(unstr_mesh), pointer :: mesh
    type(mfd_disc),   target :: disc
    type(HC_model),   target :: model
    type(HC_AE_solver) :: solver
    type(parameter_list), pointer :: HC_params, PCG_params, plist
    real(r8), allocatable :: ucell(:), uface(:), uface_ref(:)
    real(r8) :: xc(3), error
    integer :: n, j, stat
    character(:), allocatable :: errmsg

    !! 2D mesh on [0,1]^2 (one cell thick in z)
    mesh => new_unstr_mesh([0.0_r8, 0.0_r8, 0.0_r8], [1.0_r8, 1.0_r8, 0.1_r8], [5,8,1])
    call disc%init (mesh)

    !! Instantiate the HC model...
    !! First define the HC parameter list.  Only conductivity and the BC are
    !! significant.  Our test function is u(x,y) = 2x + y.  Cell values will
    !! be set using this function and we aim to recover matching face values
    !! and so need to set BC consistent with that.
    allocate(HC_params)
    call HC_params%set ('density', 2.0_r8)
    call HC_params%set ('specific-heat', 0.25_r8)
    call HC_params%set ('conductivity', 0.5_r8)
    !! Dirichlet BC on the left and bottom
    plist => HC_params%sublist ('bc')
    plist => plist%sublist ('left-bottom')
    call plist%set ('condition', 'dirichlet')
    call plist%set ('face-sets', [1,3])
    plist => plist%sublist('data-function')
    call plist%set ('type', 'polynomial')
    call plist%set ('poly-coef', [2.0_r8, 1.0_r8])
    call plist%set ('poly-powers', reshape([0,1,0, 0,0,1],shape=[3,2]))
    !! Flux BC on the right
    plist => HC_params%sublist ('bc')
    plist => plist%sublist ('right')
    call plist%set ('condition', 'flux')
    call plist%set ('face-sets', [2])
    call plist%set ('data-constant', -1.0_r8)
    !! Flux BC on the top
    plist => HC_params%sublist ('bc')
    plist => plist%sublist ('top')
    call plist%set ('condition', 'flux')
    call plist%set ('face-sets', [4])
    call plist%set ('data-constant', -0.5_r8)
    !! Symmetry BC on the z=const boundaries
    plist => HC_params%sublist ('bc')
    plist => plist%sublist ('symmetry')
    call plist%set ('condition', 'flux')
    call plist%set ('face-sets', [5,6])
    call plist%set ('data-constant', 0.0_r8)
    !! Now create the HC model.
    call model%init (disc, HC_params)

    !! Instantiate the face solver.
    allocate(PCG_params)
    call PCG_params%set ('rel-tol', 1.0d-12)
    call PCG_params%set ('max-iter', 50)
    !call PCG_params%set ('error tolerance', 1.0d-12)
    !call PCG_params%set ('max iterations', 100)
    !call PCG_params%set ('num cycles', 1)
    !call PCG_params%set ('print level', 2)
    !call PCG_params%set ('debug level', 0)
    !call PCG_params%set ('logging level', 1)
    call solver%init (model, PCG_params)

    !! Define cell values and initial guess for face values
    allocate(ucell(mesh%ncell), uface(mesh%nface), uface_ref(mesh%nface))
    do j = 1, mesh%ncell
      xc = sum(mesh%x(:,mesh%cnode(:,j)),dim=2) / size(mesh%cnode,1)
      ucell(j) = 2*xc(1) + xc(2)
    end do
    call random_number (uface)  ! initial guess

    call solver%solve (0.0_r8, ucell, uface, stat, errmsg)
    if (stat /= 0) then
      status = 1
      write(*,'(a,es9.2)') 'test1: solver returned an error: ', errmsg
    end if

    !! Compare solution to expected values.
    do j = 1, mesh%nface
      xc = sum(mesh%x(:,mesh%fnode(:,j)),dim=2) / size(mesh%fnode,1)
      uface_ref(j) = 2*xc(1) + xc(2)
    end do
    error = maxval(abs(uface-uface_ref))
    write(*,'(a,es9.2)') 'test1: ||uface-uface_ref||_max =', error
    if (error > 1.0e-15) then
      status = 1
      write(*,'(a)') 'test1: error exceeds tolerance'
    end if

  end subroutine

  !! Similar problem but with a randomized mesh to make CG do some work.
  !! Can no longer recover test function exactly because of nonorthog mesh.

  subroutine test2

    type(unstr_mesh), pointer :: mesh
    type(mfd_disc),   target :: disc
    type(HC_model),   target :: model
    type(HC_AE_solver) :: solver
    type(parameter_list), pointer :: HC_params, PCG_params, plist
    real(r8), allocatable :: ucell(:), uface(:), uface_ref(:)
    real(r8) :: xc(3), error
    integer :: n, j, stat
    character(:), allocatable :: errmsg

    !! 3D mesh on [0,1]^3
    mesh => new_unstr_mesh([real(r8)::0,0,0], [real(r8)::1,1,1], [17,15,16], eps=0.001_r8)
    call disc%init (mesh)

    !! Instantiate the HC model...
    !! First define the HC parameter list.  Only conductivity and the BC are
    !! significant.  Our test function is u(x,y,z) = x + 2y - 4z.  Cell values
    !! will be set using this function and we aim to recover matching face
    !! values and so need to set BC consistent with that.
    allocate(HC_params)
    call HC_params%set ('density', 2.0_r8)
    call HC_params%set ('specific-heat', 0.25_r8)
    call HC_params%set ('conductivity', 0.5_r8)
    !! Dirichlet BC on the left, front, and bottom
    plist => HC_params%sublist ('bc')
    plist => plist%sublist ('left-bottom')
    call plist%set ('condition', 'dirichlet')
    call plist%set ('face-sets', [1,3,5])
    plist => plist%sublist('data-function')
    call plist%set ('type', 'polynomial')
    call plist%set ('poly-coef', [1.0_r8, 2.0_r8, -4.0_r8])
    call plist%set ('poly-powers', reshape([0,1,0,0, 0,0,1,0, 0,0,0,1],shape=[4,3]))
    !! Flux BC on the right
    plist => HC_params%sublist ('bc')
    plist => plist%sublist ('right')
    call plist%set ('condition', 'flux')
    call plist%set ('face-sets', [2])
    call plist%set ('data-constant', -0.5_r8)
    !! Flux BC on the back
    plist => HC_params%sublist ('bc')
    plist => plist%sublist ('back')
    call plist%set ('condition', 'flux')
    call plist%set ('face-sets', [4])
    call plist%set ('data-constant', -1.0_r8)
    !! Flux BC on the top
    plist => HC_params%sublist ('bc')
    plist => plist%sublist ('symmetry')
    call plist%set ('condition', 'flux')
    call plist%set ('face-sets', [6])
    call plist%set ('data-constant', 2.0_r8)
    !! Now create the HC model.
    call model%init (disc, HC_params)

    !! Instantiate the face solver.
    allocate(PCG_params)
    call PCG_params%set ('rel-tol', 1.0d-12)
    call PCG_params%set ('max-iter', 50)
    call solver%init (model, PCG_params)

    !! Define cell values and initial guess for face values
    allocate(ucell(mesh%ncell), uface(mesh%nface), uface_ref(mesh%nface))
    do j = 1, mesh%ncell
      xc = sum(mesh%x(:,mesh%cnode(:,j)),dim=2) / size(mesh%cnode,1)
      ucell(j) = xc(1) + 2*xc(2) - 4*xc(3)
    end do
    call random_number (uface)  ! initial guess

    call solver%solve (0.0_r8, ucell, uface, stat, errmsg)
    if (stat /= 0) then
      status = 1
      write(*,'(a,es9.2)') 'test2: solver returned an error: ', errmsg
    end if

    !! Compare solution to expected values.
    do j = 1, mesh%nface
      xc = sum(mesh%x(:,mesh%fnode(:,j)),dim=2) / size(mesh%fnode,1)
      uface_ref(j) = xc(1) + 2*xc(2) - 4*xc(3)
    end do
    error = maxval(abs(uface-uface_ref))
    write(*,'(a,es9.2)') 'test2: ||uface-uface_ref||_max =', error
    if (error > 1.0e-4) then
      status = 1
      write(*,'(a)') 'test2: error exceeds tolerance'
    end if

  end subroutine

end program test_HC_AE_solver_type
