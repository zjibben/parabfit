!!
!! DIFFERENTIAL_OPERATORS
!!
!! This module contains a collection of small, public functions for calculation
!! various derivatives.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! March 2016
!!

#include "f90_assert.fpp"

module differential_operators

  use kinds, only: r8
  implicit none
  private

  public :: faceGradient, divergence

contains

  ! currently, this calculates the gradient of a cell-centered quantity
  ! on the face for orthogonal meshes
  ! example: dynamic pressure face gradients
  pure function faceGradient (y, xc)

    use consts, only: ndim

    real(r8), intent(in) :: y(:), xc(:,:)
    real(r8)             :: faceGradient(ndim)

    real(r8)             :: dx(ndim)

    ! ASSERT(size(y)==2)
    ! ASSERT(size(xc, dim=2)==2)

    dx = xc(:,1) - xc(:,2)
    faceGradient = (y(1) - y(2)) * dx / sum(dx**2)

  end function faceGradient

  ! calculates the divergence of a function y, which is given on the faces
  ! it is assumed that y has already been dotted with the outward normal
  ! on every face. It may be worth adding another function which does not
  ! make this assumption
  real(r8) pure function divergence (y, face_area, cell_vol)
    real(r8), intent(in) :: y(:), face_area(:), cell_vol
    divergence = sum(y*face_area) / cell_vol
  end function divergence

end module differential_operators
