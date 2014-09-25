!!
!! CSR_PRECON_SSOR_TYPE
!!
!! A concrete implementation of the abstract base class CSR_PRECON that
!! uses SSOR preconditioning.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! March 2014, adapted for F2008
!!
!! PROGRAMMING INTERFACE
!!
!!  This module defines the derived type CSR_PRECON_SSOR which is an
!!  extension of the abstract base class CSR_PRECON that implements SSOR
!!  preconditioning.  See the base class comments for a description of the
!!  common type bound procedures.
!!
!!  The INIT procedure expects to find the following parameters in the
!!  TYPE(PARAMETER_LIST) argument PARAMS.  Parameters with a default value
!!  are optional; the others are required.
!!
!!    'number of sweeps'  - The number of SSOR sweeps; > 0.
!!    'relaxation factor' - The (over) relaxation factor (default 1.0); > 0.0
!!
!! IMPLEMENTATION NOTES
!!
!! TODO: proper error handling of parameter list values.
!!

#include "f90_assert.fpp"

module csr_precon_ssor_type

  use kinds, only: r8
  use csr_matrix_type
  use csr_precon_class
  use parameter_list_type
  implicit none
  private

  type, extends(csr_precon), public :: csr_precon_ssor
    private
    real(r8), allocatable :: diag(:)
    integer :: num_iter
    real(r8) :: omega
  contains
    procedure :: init
    procedure :: compute
    procedure :: apply
  end type csr_precon_ssor

contains

  subroutine init (this, A, params)

    use logging_services

    class(csr_precon_ssor), intent(out) :: this
    type(csr_matrix), intent(in), target :: A
    type(parameter_list) :: params

    integer :: stat
    character(:), allocatable :: context, errmsg

    this%A => A

    allocate(this%diag(A%nrow))

    context = 'processing ' // params%name() // ': '
    call params%get ('num-sweeps', this%num_iter, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (this%num_iter <= 0) call LS_fatal (context//'"num-sweeps" must be > 0')
    call params%get ('omega', this%omega, default=1.0_r8, stat=stat, errmsg=errmsg)
    if (stat /= 0) call LS_fatal (context//errmsg)
    if (this%omega <= 0.0_r8) call LS_fatal (context//'"omega" must be > 0.0')

  end subroutine init

  subroutine compute (this)
    class(csr_precon_ssor), intent(inout) :: this
    call this%A%get_diag_copy (this%diag)
  end subroutine compute

  subroutine apply (this, x)

    class(csr_precon_ssor), intent(in) :: this
    real(r8), intent(inout) :: x(:)

    integer :: n, j, k
    real(r8) :: s, u(this%A%nrow)

    ASSERT(size(x) >= this%A%nrow)

    u = 0.0_r8
    do n = 1, this%num_iter
      !! Forward sweep.
      do j = 1, this%A%nrow
        s = x(j)
        do k = this%A%graph%xadj(j), this%A%graph%xadj(j+1)-1
          s = s - this%A%values(k) * u(this%A%graph%adjncy(k))
        end do
        u(j) = u(j) + this%omega * (s / this%diag(j))
      end do
      !! Backward sweep.
      do j = this%A%nrow, 1, -1
        s = x(j)
        do k = this%A%graph%xadj(j), this%A%graph%xadj(j+1)-1
          s = s - this%A%values(k) * u(this%A%graph%adjncy(k))
        end do
        u(j) = u(j) + this%omega * (s / this%diag(j))
      end do
    end do

    !! Copy solution to return array.
    x(:this%A%nrow) = u

  end subroutine apply

end module csr_precon_ssor_type
