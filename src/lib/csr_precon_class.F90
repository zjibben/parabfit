!!
!! CSR_PRECON_CLASS
!!
!! A common interface to preconditioners for CSR format matrices.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! March 2014
!!
!! PROGRAMMING INTERFACE
!!
!!  This module defines the abstract base class CSR_PRECON which defines a
!!  common interface to preconditioners for CSR_MATRIX-type matrices.
!!  Application code is expected to use polymorphic variables of this type
!!  and not work directly with extensions of the type.  The base type defines
!!  The following type bound procedures
!!
!!  INIT(A, PARAMS) configures the object to be a preconditioner for the
!!    specified CSR_MATRIX A.  The object holds a reference to the matrix A,
!!    and so the matrix must not go out of scope during the lifetime of the
!!    object.  Moreover, the actual argument must be a pointer or have the
!!    target attribute.  Only the structure of the matrix needs to be defined
!!    at this point; the matrix values are not referenced.  Also note that
!!    matrix will not be modified in any way by this, or any other procedure.
!!    The PARAMETER_LIST object PARAMS gives the parameters associated with
!!    the preconditioner; its expected content will depend on the dynamic
!!    type of the CSR_PRECON object (see the comments for the extended types.)
!!
!!  COMPUTE() performs the final setup and configuration of the preconditioner.
!!    It is at this point that the values of the CSR_MATRIX A are referenced.
!!    This must be called before calling the APPLY procedure and after the
!!    matrix values are defined, and must be called again whenever the matrix
!!    values are modified.
!!
!!  APPLY(X) applies the preconditioner to the vector X.  The size of X must
!!    be at least A%NROW; only the initial A%NROW elements will be referenced
!!    or modified.
!!
!!  MATRIX() returns a pointer the the CSR_MATRIX object A with which the
!!    preconditioner was configured in the INIT call.
!!

module csr_precon_class

  use kinds, only: r8
  use csr_matrix_type
  use parameter_list_type
  implicit none
  private

  type, abstract, public :: csr_precon
    type(csr_matrix), pointer :: A => null()  ! reference only - do not own
  contains
    procedure(init),    deferred :: init
    procedure(compute), deferred :: compute
    procedure(apply),   deferred :: apply
    procedure :: matrix
  end type csr_precon

  abstract interface
    subroutine init (this, A, params)
      import csr_precon, csr_matrix, parameter_list
      class(csr_precon), intent(out) :: this
      type(csr_matrix), intent(in), target :: A
      type(parameter_list) :: params
    end subroutine
    subroutine compute (this)
      import csr_precon
      class(csr_precon), intent(inout) :: this
    end subroutine
    subroutine apply (this, x)
      import csr_precon, r8
      class(csr_precon), intent(in) :: this
      real(r8), intent(inout) :: x(:)
    end subroutine
  end interface

contains

  function matrix (this)
    class(csr_precon), intent(in) :: this
    type(csr_matrix), pointer :: matrix
    matrix => this%A
  end function matrix

end module csr_precon_class
