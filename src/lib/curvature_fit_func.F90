!!
!! curvature_fit
!!
!! This module defines a function for calculating curvature from a
!! collection of points using a quadratic fitting procedure.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! December 2016
!!

#include "f90_assert.fpp"

module curvature_fit_func

  use kinds, only: r8
  use logging_services
  implicit none
  private

  public :: curvature_fit

contains

  real(r8) function curvature_fit ()

    curvature_fit = 0.0_r8

  end function curvature_fit

end module curvature_fit_func
