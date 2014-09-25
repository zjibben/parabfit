!!
!! MFD_DIFF_MATRIX
!!
!! This module defines a derived type that encapsulates the mimetic finite
!! difference (MFD) diffusion matrix.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! Adapted for F2008, March 2014
!!
!! PROGRAMMING INTERFACE
!!
!! The module defines the derived type MFD_DIFF_MATRIX.  Objects of this type
!! are initialized using one of the following type bound subroutines:
!!
!!  INIT (DISC) initializes the matrix object.  DISC is a MFD_DISC object that
!!    defines the discrete context for the matrix. The object holds a reference
!!    to DISC, and so the actual argument must either be a pointer or have the
!!    target attribute, and must persist for the lifetime of the object.
!!    The data components of the object are allocated.
!!
!!  INIT (MOLD) initializes the matrix object using another MFD_DIFF_MATRIX
!!    object MOLD as a template.  The matrix will share the same underlying
!!    MFD_DISC object as MOLD.
!!
!! NB: A matrix object initialized using the first call should not be finalized
!! until all other matrix objects initialized using it as a mold, directly or
!! indirectly, have been finalized.  All such objects share a common CRS_GRAPH
!! object that the first matrix owns and will finalize.
!!
!! The following type bound procedures are defined.  COMPUTE evaluates the
!! basic diffusion matrix elements, and the remaining procedures perform
!! specific fix-ups to the basic matrix associated with boundary conditions
!! and other discretization details.  If used, COMPUTE_FACE_SCHUR_MATRIX
!! should be called after all the necessary fix-ups.
!!
!!  COMPUTE (COEF) computes (or re-computes) the diffusion matrix using the
!!    given cell-based array COEF of diffusion coefficients.  This resets the
!!    internal list of Dirichlet faces, if any, to empty.
!!
!!  SET_DIR_FACES (FACES) sets the faces listed in the rank-1 array FACES to
!!    be Dirichlet faces.  These are added to the existing set of Dirichlet
!!    faces, if any.  The implication is that the unknowns associated with
!!    these faces are not actually unknowns and should not be included in the
!!    diffusion matrix.  This is handled by projecting out the corresponding
!!    face rows and columns and setting unit diagonal values for those faces,
!!    effectively replacing equations for those unknowns with dummy equations
!!    that are decoupled from all other unknowns.
!!
!!  INCR_CELL_DIAG (VALUES) increments the diagonal cell-cell diffusion
!!    submatrix with the given cell-based array VALUES.  This is used to
!!    incorporate time derivative terms into the diffusion matrix.
!!
!!  INCR_FACE_DIAG (INDICES, VALUES) increments selected diagonal elements
!!    of the face-face submatrix.  The rank-1 array INDICES is a list of face
!!    indices and the rank-1 array VALUES is the the corresponding list of
!!    values to be added to the diagonal.  This is used to implement certain
!!    flux-type boundary conditions.
!!
!!  COMPUTE_FACE_SCHUR_MATRIX (SFF) computes the face-face Schur complement
!!    matrix.  SFF is an intent(inout) CSR_MATRIX object that has the same
!!    non-zero structure as the A22 submatrix component of the diffusion
!!    matrix object.  It must be defined on entry (values are ignored), and
!!    its values are overwritten with the face-face Schur complement of the
!!    diffusion matrix.
!!

#include "f90_assert.fpp"

module mfd_diff_matrix_type

  use kinds, only: r8
  use unstr_mesh_type
  use mfd_disc_type
  use csr_matrix_type
  implicit none
  private

  type, public :: mfd_diff_matrix
    type(mfd_disc),   pointer :: disc => null() ! reference only - do not own
    type(unstr_mesh), pointer :: mesh => null() ! reference only - do not own
    real(r8), allocatable :: a11(:)     ! the cell-cell submatrix
    real(r8), allocatable :: a12(:,:)   ! the cell-face submatrix
    type(csr_matrix)      :: a22        ! the face-face submatrix
    integer, allocatable :: dir_faces(:)
  contains
    procedure, private :: init_disc
    procedure, private :: init_mold
    generic   :: init => init_disc, init_mold
    procedure :: compute
    procedure :: set_dir_faces
    procedure :: incr_cell_diag
    procedure :: incr_face_diag
    procedure :: compute_face_schur_matrix
    procedure :: forward_elimination
    procedure :: backward_substitution
  end type mfd_diff_matrix

contains

  !! Initialize starting from a MFD_DISC object.
  subroutine init_disc (this, disc)
    class(mfd_diff_matrix), intent(out) :: this
    type(mfd_disc), intent(in), target :: disc
    integer :: j
    type(csr_graph), pointer :: g
    this%disc => disc
    this%mesh => disc%mesh
    allocate(this%a11(this%mesh%ncell))
    allocate(this%a12(size(this%mesh%cface,1),this%mesh%ncell))
    !! Create a CSR matrix graph for the A22 submatrix.
    allocate(g)
    call g%init (this%mesh%nface)
    do j = 1, this%mesh%ncell
      call g%add_clique (this%mesh%cface(:,j))
    end do
    call g%add_complete
    !! Create the A22 submatrix.
    call this%a22%init (g, take_graph=.true.)
  end subroutine init_disc

  !! Initialize using another MFD_DIFF_MATRIX object as a template.
  subroutine init_mold (this, mold)
    class(mfd_diff_matrix), intent(out) :: this
    class(mfd_diff_matrix), intent(in)  :: mold
    type(csr_graph), pointer :: g
    this%disc => mold%disc
    this%mesh => mold%mesh
    allocate(this%a11(size(mold%a11)))
    allocate(this%a12(size(mold%a12,1),size(mold%a12,2)))
    g => mold%a22%graph
    call this%a22%init (g, take_graph=.false.)
  end subroutine init_mold

  !! Compute the diffusion matrix given the cell-based array of diffusion
  !! coefficients.  Any existing dirichlet faces are dropped
  subroutine compute (this, coef)

    use upper_packed_matrix, only: upm_col_sum

    class(mfd_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: coef(:)

    integer :: j, n, ir, ic
    real(r8) :: w(size(this%mesh%cface,dim=1))
    real(r8) :: minv(size(w)*(size(w)+1)/2)

    ASSERT(size(coef) == this%mesh%ncell)

    call this%a22%set_all (0.0_r8)
    do j = 1, this%mesh%ncell
      if (coef(j) == 0.0_r8) then
        this%a11(j) = 0.0_r8
        this%a12(:,j) = 0.0_r8
        cycle
      end if
      minv = coef(j) * this%disc%minv(:,j)
      !! Fill the A11 and A12 submatrices
      call upm_col_sum (minv, w)
      this%a11(j) = sum(w)
      this%a12(:,j) = -w
      !! Assemble the A22 CSR submatrix.
      associate (index => this%mesh%cface(:,j))
        n = 1
        do ic = 1, size(index)
          do ir = 1, ic-1
            call this%a22%increment (index(ir), index(ic), minv(n))
            call this%a22%increment (index(ic), index(ir), minv(n))
            n = n + 1
          end do
          call this%a22%increment (index(ic), index(ic), minv(n))
          n = n + 1
        end do
      end associate
    end do

    if (allocated(this%dir_faces)) deallocate(this%dir_faces)

  end subroutine compute

  !! Set the specified faces to be Dirichlet faces.  These are added to the
  !! existing set of Dirichlet faces, if any.  Unknowns associated with these
  !! faces are not actually unknowns and should no be included in the diffusion
  !! matrix.  This is handled by projecting out the corresponding face rows and
  !! columns and setting unit diagonal values for those faces, effectively
  !! replacing equations for those unknowns with dummy equations that are
  !! decoupled from all other unknowns.  The face-face submatrix is directly
  !! modified.  The modified cell-face submatrix is kept in factored form as
  !! the product of the original cell-face submatrix and the projection matrix
  !! described by the list of Dirichlet faces; this is more efficient than
  !! directly modifying the cell-face submatrix.

  subroutine set_dir_faces (this, dir_faces)

    class(mfd_diff_matrix), intent(inout) :: this
    integer, intent(in) :: dir_faces(:)

    integer :: j, n
    integer, allocatable :: tmp(:)

    ASSERT(minval(dir_faces) >= 1)
    ASSERT(maxval(dir_faces) <= this%mesh%nface)

    if (allocated(this%dir_faces)) then
      if (size(dir_faces) > 0) then
        n = size(this%dir_faces) + size(dir_faces)
        allocate(tmp(n))
        n = size(this%dir_faces)
        tmp(:n) = this%dir_faces
        tmp(n+1:) = dir_faces
        call move_alloc (tmp, this%dir_faces)
      end if
    else
      this%dir_faces = dir_faces
    end if

    do j = 1, size(dir_faces)
      n = dir_faces(j)
      call this%a22%project_out (n)
      call this%a22%set (n, n, 1.0_r8)
    end do

  end subroutine set_dir_faces

  !! Increment the (entire) diagonal cell-cell diffusion submatrix with
  !! the specified values.  The intended use is the incorporation of time
  !! derivative terms into the base diffusion matrix.

  subroutine incr_cell_diag (this, values)
    class(mfd_diff_matrix), intent(inout) :: this
    real(r8), intent(in) :: values(:)
    ASSERT(size(values) == size(this%a11))
    this%a11 = this%a11 + values
  end subroutine incr_cell_diag

  !! Increment selected diagonal elements of the face-face diffusion submatrix.
  !! The array INDICES is a list of face indices and VALUES the corresponding
  !! array of values to be added to the diagonal.  The intended use is the
  !! implementation of certain types of flux boundary condtions.

  subroutine incr_face_diag (this, indices, values)
    class(mfd_diff_matrix), intent(inout) :: this
    integer,  intent(in) :: indices(:)
    real(r8), intent(in) :: values(:)
    integer :: j
    ASSERT(size(indices) == size(values))
    ASSERT(minval(indices) >= 1)
    ASSERT(maxval(indices) <= this%mesh%nface)
    do j = 1, size(indices)
      call this%a22%increment (indices(j), indices(j), values(j))
    end do
  end subroutine incr_face_diag

  !! Compute the face Schur complement matrix SFF.  SFF is intent inout.  It
  !! must be defined on input (structure only), and its values are overwritten
  !! with the values of the Schur complement.  Its structure must match the
  !! structure of the face-face diffusion submatrix (same CSR_GRAPH component).

  subroutine compute_face_schur_matrix (this, Sff)

    class(mfd_diff_matrix), intent(in) :: this
    type(csr_matrix), intent(inout) :: Sff

    integer :: j, n, ir, ic
    real(r8) :: tmp

    ASSERT(associated(this%a22%graph, Sff%graph))

    Sff%values = this%a22%values
    do j = 1, this%mesh%ncell
      associate (indices => this%mesh%cface(:,j))
        do ir = 1, size(indices)
          do ic = 1, size(indices)
            tmp = -this%a12(ir,j)*this%a12(ic,j)/this%a11(j)
            call Sff%increment (indices(ir), indices(ic), tmp)
          end do
        end do
      end associate
    end do

    !! Apply the Dirichlet projections.
    if (allocated(this%dir_faces)) then
      do j = 1, size(this%dir_faces)
        n = this%dir_faces(j)
        call Sff%project_out (n)
        call Sff%set (n, n, 1.0_r8)
      end do
    end if

  end subroutine compute_face_schur_matrix

  subroutine forward_elimination (this, b1, b2)

    class(mfd_diff_matrix), intent(in) :: this
    real(r8), intent(in) :: b1(:)
    real(r8), intent(inout) :: b2(:)

    integer :: j, k, n
    real(r8) :: s
    real(r8), allocatable :: b2_dir(:)

    ASSERT(size(b1) == this%mesh%ncell)
    ASSERT(size(b2) == this%mesh%nface)

    if (allocated(this%dir_faces)) b2_dir = b2(this%dir_faces)

    do j = 1, size(this%a12,dim=2)
      s = b1(j) / this%a11(j)
      do k = 1, size(this%a12,dim=1)
        n = this%mesh%cface(k,j)
        b2(n) = b2(n) - this%a12(k,j) * s
      end do
    end do

    if (allocated(this%dir_faces)) b2(this%dir_faces) = b2_dir

  end subroutine forward_elimination

  subroutine backward_substitution (this, b1, u2)

    class(mfd_diff_matrix), intent(in) :: this
    real(r8), intent(inout) :: b1(:), u2(:)

    integer :: j, k
    real(r8) :: s
    real(r8), allocatable :: u2_dir(:)

    if (allocated(this%dir_faces)) then
      u2_dir = u2(this%dir_faces)
      u2(this%dir_faces) = 0.0_r8
    end if

    do j = 1, this%mesh%ncell
      s = b1(j)
      do k = 1, size(this%a12,dim=1)
        s = s - this%a12(k,j) * u2(this%mesh%cface(k,j))
      end do
      b1(j) = s / this%a11(j)
    end do

    if (allocated(this%dir_faces)) u2(this%dir_faces) = u2_dir

  end subroutine backward_substitution

end module mfd_diff_matrix_type
