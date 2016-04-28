!!
!! FISCHER_TYPE
!!
!! Implements the Fischer intial guess technique from
!! 
!! P. F. Fischer, Projection Techniques for iterative solution of Ax=b with
!!  successive right-hand sides.  Comput. Methods Appl. Mech. Engrg. v. 163
!!  pp. 193--204, 1998
!!
!! TODO: * threading
!!
!! Zechariah J Jibben <zjibben@lanl.gov>
!! April 2016
!!

#include "f90_assert.fpp"

module fischer_guess_type

  use kinds, only: r8
  implicit none
  private

  type, public :: fischer_guess
    private
    real(r8), allocatable :: b_tilde(:,:), x_tilde(:,:), x_guess(:)
    integer :: Nvecs, max_Nvecs
  contains
    procedure :: init
    procedure :: guess
    procedure :: update
  end type fischer_guess

contains

  subroutine init (this,ncell)

    class(fischer_guess), intent(out) :: this
    integer, intent(in)  :: ncell

    ! This is currently hard-wired, as it is in Truchas.
    ! I'm not sure if it would be useful to make it a
    ! user-specified parameter or not. -zjibben
    this%max_Nvecs = 6

    allocate(this%b_tilde(ncell,this%max_Nvecs), this%x_tilde(ncell,this%max_Nvecs), &
        this%x_guess(ncell))
    
    this%Nvecs = 0
    this%x_guess = 0.0_r8
    this%b_tilde = 0.0_r8
    this%x_tilde = 0.0_r8

  end subroutine init

  ! Given the new right hand side b, this computes a new guess for
  ! the solution x based on previous stored solution values
  !
  ! CAUTION: The space of vectors is updated after the true solution
  !  to Ax=b is found.  That updated requires the guess we calculate
  !  here.  As such, the guess is stored in the module (x_guess).  
  !  The usual code flow will be
  !
  !    call Fischer_Initial_Guess( x', b)
  !    solve Ax=b using x'
  !    call Fischer_update_space( x, b )
  !
  !  If Fischer_Initial_Guess is called repeatedly before the update
  !  then the update will be wrong.
  subroutine guess (this, x, b)

    class(fischer_guess), intent(inout) :: this
    real(r8),       intent(out)   :: x(:)
    real(r8),       intent(in)    :: b(:)

    real(r8) :: alpha(this%Nvecs)
    integer :: i

    ASSERT(size(x)==size(b))
    ASSERT(size(x)==size(this%x_guess))

    do i = 1,this%Nvecs
      alpha(i) = sum(b * this%b_tilde(:,i))
    end do

    do i = 1,size(x)
      this%x_guess(i) = sum(alpha*this%x_tilde(i,:))
    end do

    x = this%x_guess

  end subroutine guess

  ! After a PPE solution, the set of projection vectors
  ! (i.e. previous solutions and right hand sides) is updated.
  subroutine update (this, x_soln, b, lhs)

    use consts, only: alittle
    use array_utils, only: magnitude
    use csr_matrix_type

    class(fischer_guess),   intent(inout) :: this
    real(r8),         intent(in)    :: x_soln(:), b(:)
    type(csr_matrix), intent(in)    :: lhs

    integer :: i
    real(r8) :: b_norm, alpha(this%max_Nvecs)

    if (this%Nvecs == this%max_Nvecs) then

      this%x_tilde(:,1) = x_soln
      this%b_tilde(:,1) = b !lhs%matvec(this%x_tilde(:,1))

      b_norm = magnitude(this%b_tilde(:,1))
      if (b_norm < alittle) return ! disregard
      this%b_tilde(:,1) = this%b_tilde(:,1) / b_norm
      this%x_tilde(:,1) = this%x_tilde(:,1) / b_norm
      this%Nvecs = 1
    else
      this%Nvecs = this%Nvecs + 1

      this%x_tilde(:,this%Nvecs) = x_soln - this%x_guess
      this%b_tilde(:,this%Nvecs) = lhs%matvec(this%x_tilde(:,this%Nvecs))
      
      do i = 1,this%Nvecs - 1
        alpha(i) = sum(this%b_tilde(:,this%Nvecs) * this%b_tilde(:,i))
      end do

      do i = 1,size(x_soln)
        this%b_tilde(i,this%Nvecs) = this%b_tilde(i,this%Nvecs) - sum(alpha*this%b_tilde(i,:))
        this%x_tilde(i,this%Nvecs) = this%x_tilde(i,this%Nvecs) - sum(alpha*this%x_tilde(i,:))
      end do

      b_norm = magnitude(this%b_tilde(:,this%Nvecs))
      if (b_norm < alittle) then ! disregard the vector
        this%Nvecs = this%Nvecs - 1
        return
      end if
      this%b_tilde(:,this%Nvecs) = this%b_tilde(:,this%Nvecs) / b_norm
      this%x_tilde(:,this%Nvecs) = this%x_tilde(:,this%Nvecs) / b_norm
    end if

  end subroutine update

end module fischer_guess_type
