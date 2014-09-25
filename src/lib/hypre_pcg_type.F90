!!
!! HYPRE_PCG_TYPE
!!
!! This module defines a preconditioned CG solver class built on the
!! ParCSR PCG solver from Hypre that uses BoomerAMG as the preconditioner.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!!

#include "f90_assert.fpp"

module hypre_pcg_type

  use kinds, only: r8
  use fhypre
  use csr_matrix_type
  use parameter_list_type
  implicit none
  private

  type, public :: hypre_pcg
    private
    type(csr_matrix), pointer :: A => null()
    integer :: nrows = 0, ilower = 0, iupper = 0
    type(hypre_obj) :: solver = hypre_null_obj    ! HYPRE_Solver object handle
    type(hypre_obj) :: precon = hypre_null_obj    ! HYPRE_Solver object handle
    type(hypre_obj) :: Ah = hypre_null_obj        ! HYPRE_IJMatrix object handle
    type(hypre_obj) :: bh = hypre_null_obj, xh = hypre_null_obj ! HYPRE_IJVector object handles
    !! PCG parameters -- these are set at initialization
    integer  :: print_level       ! ??? OFF=0, SETUP=1, SOLVE=2, SETUP+SOLVE=3
    integer  :: max_iter          ! maximum number of iterations
    real(r8) :: err_tol           ! error tolerance (||r||_2 / ||b||_2)
    !! BoomerAMG preconditioner parameters -- these are set at initialization
    integer  :: num_cycles        ! number of cycles
    integer  :: debug_level       ! OFF=0, ON=1
    integer  :: logging_level     ! OFF=0, ON=1, >1=residual available from hypre
    !! BoomerAMG preconditioner parameters -- these are currently hardwired
    integer  :: coarsen_type = 6  ! Falgout coarsening
    integer  :: relax_type = 3    ! hybrid Gauss-Seidel smoothing
    integer  :: num_sweeps = 1    ! number of smoother sweeps
    integer  :: max_levels = 25   ! max number of multigrid levels
    real(r8) :: strong_threshold = 0.5_r8 ! should be 0.5 for 3D problems and 0.25 for 2D
  contains
    procedure :: init
    procedure :: compute
    procedure :: solve
    procedure :: matrix
    procedure :: get_metrics
    final :: hypre_pcg_delete
  end type hypre_pcg

contains

  !! Final subroutine for HYPRE_PCG objects.
  subroutine hypre_pcg_delete (this)
    type(hypre_pcg) :: this
    integer :: ierr
    ierr = 0
    call fHYPRE_ClearAllErrors
    if (hypre_associated(this%Ah)) call fHYPRE_IJMatrixDestroy (this%Ah, ierr)
    if (hypre_associated(this%bh)) call fHYPRE_IJVectorDestroy (this%bh, ierr)
    if (hypre_associated(this%xh)) call fHYPRE_IJVectorDestroy (this%xh, ierr)
    if (hypre_associated(this%solver)) call fHYPRE_PCGDestroy (this%solver, ierr)
    if (hypre_associated(this%precon)) call fHYPRE_BoomerAMGDestroy (this%precon, ierr)
    INSIST(ierr == 0)
  end subroutine hypre_pcg_delete

  subroutine init (this, A, params)

    class(hypre_pcg), intent(out) :: this
    type(csr_matrix), intent(in), target :: A
    type(parameter_list) :: params

    integer :: ierr

    this%A => A

    this%nrows  = A%nrow    ! number of on-process rows (if parallel)
    this%ilower = 1         ! global index of first on-process row (if parallel)
    this%iupper = A%nrow    ! global index of last on-process row (if parallel)

    call fHYPRE_ClearAllErrors

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%bh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%bh, 0, ierr)
    INSIST(ierr == 0)

    call fHYPRE_IJVectorCreate (this%ilower, this%iupper, this%xh, ierr)
    call fHYPRE_IJVectorSetMaxOffProcElmts (this%xh, 0, ierr)
    INSIST(ierr == 0)

    call params%get ('error tolerance', this%err_tol)
    INSIST(this%err_tol > 0.0_r8)
    call params%get ('max iterations', this%max_iter)
    INSIST(this%max_iter > 0)
    call params%get ('num cycles', this%num_cycles)
    INSIST(this%num_cycles > 0)
    call params%get ('print level', this%print_level, default=0)
    INSIST(this%print_level >=0 .and. this%print_level <= 3)
    call params%get ('debug level', this%debug_level, default=0)
    INSIST(this%debug_level >= 0)
    call params%get ('logging level', this%logging_level, default=0)
    INSIST(this%logging_level >= 0)

  end subroutine init

  function matrix (this)
    class(hypre_pcg), intent(in) :: this
    type(csr_matrix), pointer :: matrix
    matrix => this%A
  end function matrix

  subroutine get_metrics (this, num_itr)
    class(hypre_pcg), intent(in) :: this
    integer, intent(out), optional :: num_itr
    integer :: ierr
    INSIST(hypre_associated(this%solver))
    if (present(num_itr)) then
      call fHYPRE_PCGGetNumIterations (this%solver, num_itr, ierr)
      INSIST(ierr == 0)
    end if
  end subroutine get_metrics

  subroutine compute (this)

    class(hypre_pcg), intent(inout) :: this

    integer :: ierr, i
    real(r8) :: dummy(this%nrows)
    integer :: rows(this%nrows)

    ASSERT(associated(this%A))
    call copy_to_ijmatrix (this%A, this%Ah)

    !! Create the Hypre BoomerAMG solver object (preconditioner).  Note that
    !! once the solver has been setup, it is not possible to change the matrix
    !! values without completely destroying the solver and recreating it.
    if (hypre_associated(this%precon)) call fHYPRE_BoomerAMGDestroy (this%precon, ierr)
    call fHYPRE_BoomerAMGCreate (this%precon, ierr)
    INSIST(ierr == 0)

    !! Set some debugging/diagnostic output options.
    call fHYPRE_BoomerAMGSetDebugFlag (this%precon, this%debug_level, ierr)
    call fHYPRE_BoomerAMGSetLogging (this%precon, this%logging_level, ierr)
    INSIST(ierr == 0)

    !! Set the Boomer AMG parameters.
    call fHYPRE_BoomerAMGSetPrintLevel  (this%precon, this%print_level, ierr)
    call fHYPRE_BoomerAMGSetCoarsenType (this%precon, this%coarsen_type, ierr)
    call fHYPRE_BoomerAMGSetRelaxType   (this%precon, this%relax_type, ierr)
    call fHYPRE_BoomerAMGSetNumSweeps   (this%precon, this%num_sweeps, ierr)
    call fHYPRE_BoomerAMGSetMaxLevels   (this%precon, this%max_levels, ierr)
    call fHYPRE_BoomerAMGSetMaxIter     (this%precon, this%num_cycles, ierr)
    call fHYPRE_BoomerAMGSetTol         (this%precon, 0.0_r8, ierr)
    call fHYPRE_BoomerAMGSetStrongThreshold  (this%precon, this%strong_threshold, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%precon, 3, 1, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%precon, 4, 2, ierr)
    call fHYPRE_BoomerAMGSetCycleRelaxType (this%precon, 9, 3, ierr)
    INSIST(ierr == 0)

    !! Create the Hypre PCG solver object. This supposes that once the solver
    !! has been setup, it is not possible to change the matrix values without
    !! completely destroying the solver and recreating it.  VERIFY THIS.
    if (hypre_associated(this%solver)) call fHYPRE_PCGDestroy (this%solver, ierr)
    call fHYPRE_PCGCreate (this%solver, ierr)
    INSIST(ierr == 0)

    !! Set the PCG parameters.
    call fHYPRE_PCGSetPrintLevel (this%solver, this%print_level, ierr)
    call fHYPRE_PCGSetTwoNorm (this%solver, 1, ierr)
    call fHYPRE_PCGSetTol (this%solver, this%err_tol, ierr)
    call fHYPRE_PCGSetAbsoluteTol (this%solver, 0.0_r8, ierr)
    call fHYPRE_PCGSetMaxIter (this%solver, this%max_iter, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre solution and RHS vectors;
    !! structure (not values) needed to setup PCG solver.
    rows = [ (i, i = this%ilower, this%iupper) ]  ! global row indices for this process.
    dummy = 0.0_r8
    call fHYPRE_IJVectorInitialize (this%bh, ierr)
    call fHYPRE_IJVectorSetValues  (this%bh, this%nrows, rows, dummy, ierr)
    call fHYPRE_IJVectorAssemble   (this%bh, ierr)
    INSIST(ierr == 0)
    call fHYPRE_IJVectorInitialize (this%xh, ierr)
    call fHYPRE_IJVectorSetValues  (this%xh, this%nrows, rows, dummy, ierr)
    call fHYPRE_IJVectorAssemble   (this%xh, ierr)
    INSIST(ierr == 0)

    !! After setup the solver is ready to go.
    call fHYPRE_PCGSetPrecond (this%solver, this%precon, ierr)
    call fHYPRE_PCGSetup (this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)

  end subroutine compute

  subroutine solve (this, b, x)

    class(hypre_pcg), intent(inout) :: this
    real(r8), intent(in) :: b(:)
    real(r8), intent(inout) :: x(:)

    integer :: i, ierr, rows(this%nrows)

    call fHYPRE_ClearAllErrors

    !! Global row indices for this process.
    rows = [ (i, i = this%ilower, this%iupper) ]

    !! Initialize the Hypre RHS vector.
    call fHYPRE_IJVectorInitialize (this%bh, ierr)
    call fHYPRE_IJVectorSetValues  (this%bh, this%nrows, rows, b, ierr)
    call fHYPRE_IJVectorAssemble   (this%bh, ierr)
    INSIST(ierr == 0)

    !! Initialize the Hypre initial guess vector.
    call fHYPRE_IJVectorInitialize (this%xh, ierr)
    call fHYPRE_IJVectorSetValues  (this%xh, this%nrows, rows, x, ierr)
    call fHYPRE_IJVectorAssemble   (this%xh, ierr)
    INSIST(ierr == 0)

    !! Solve the system.
    call fHYPRE_PCGSolve (this%solver, this%Ah, this%bh, this%xh, ierr)
    INSIST(ierr == 0)

    !! Retrieve the solution vector from HYPRE.
    call fHYPRE_IJVectorGetValues (this%xh, this%nrows, rows, x, ierr)
    INSIST(ierr == 0)

  end subroutine solve

 !!
 !! This auxillary routine copies a CSR_MATRIX object SRC to an equivalent
 !! HYPRE_IJMatrix object.  The HYPRE matrix is created if it does not exist.
 !! Otherwise the elements of the existing HYPRE matrix are overwritten with
 !! the values from SRC.  In the latter case the sparsity pattern of the two
 !! matrices must be identical.
 !!

  subroutine copy_to_ijmatrix (src, matrix)

    type(csr_matrix), intent(in) :: src
    type(hypre_obj), intent(inout) :: matrix

    integer :: j, ierr, ilower, iupper, nrows, nnz
    integer, allocatable :: ncols_onP(:), ncols_offP(:), ncols(:), rows(:), cols(:)

    nrows  = src%nrow
    ilower = 1
    iupper = src%nrow

    call fHYPRE_ClearAllErrors

    if (.not.hypre_associated(matrix)) then
      call fHYPRE_IJMatrixCreate (ilower, iupper, ilower, iupper, matrix, ierr)
      !! For each row we know how many column entries are on-process and how many
      !! are off-process.  HYPRE is allegedly much faster at forming its CSR matrix
      !! if it knows this info up front.
      allocate(ncols_onP(nrows), ncols_offP(nrows))
      do j = 1, nrows
        ncols_offP(j) = count(src%graph%adjncy(src%graph%xadj(j):src%graph%xadj(j+1)-1) > nrows)
        ncols_onP(j)  = src%graph%xadj(j+1) - src%graph%xadj(j) - ncols_offP(j)
      end do
      call fHYPRE_IJMatrixSetDiagOffdSizes (matrix, ncols_onP, ncols_offP, ierr)
      deallocate(ncols_onP, ncols_offP)
      !! Let HYPRE know that we won't be setting any off-process matrix values.
      call fHYPRE_IJMatrixSetMaxOffProcElmts (matrix, 0, ierr)
      INSIST(ierr == 0)
    end if

    !! After initialization the HYPRE matrix elements can be set.
    call fHYPRE_IJMatrixInitialize (matrix, ierr)
    INSIST(ierr == 0)

    !! Copy the matrix elements into the HYPRE matrix.  This defines both the
    !! nonzero structure of the matrix and the values of those elements. HYPRE
    !! expects global row and column indices.
    nnz = src%graph%xadj(nrows+1) - src%graph%xadj(1)
    allocate(ncols(nrows), rows(nrows), cols(nnz))
    rows = [ (j, j = ilower, iupper) ]
    ncols = src%graph%xadj(2:nrows+1) - src%graph%xadj(1:nrows)
    !cols = global_index(src%graph%row_ip, src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1))
    cols = src%graph%adjncy(src%graph%xadj(1):src%graph%xadj(nrows+1)-1)
    call fHYPRE_IJMatrixSetValues (matrix, nrows, ncols, rows, cols, src%values, ierr)
    deallocate(ncols, rows, cols)
    INSIST(ierr == 0)

    !! After assembly the HYPRE matrix is ready to use.
    call fHYPRE_IJMatrixAssemble (matrix, ierr)
    INSIST(ierr == 0)

  end subroutine copy_to_ijmatrix

end module hypre_pcg_type
