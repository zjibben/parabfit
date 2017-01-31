!!
!! interface_patch_type
!!
!! This module defines a class for calculating curvature from a "patch" of
!! VOF interface reconstructions. It takes a collection of interface reconstructions
!! from a cell and its neighbors, finds a set of points that represents them,
!! and hands those points off to an analytic_surface_type for fitting and
!! curvature calculation.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! January 2017
!!

#include "f90_assert.fpp"

module interface_patch_type

  use kinds, only: r8
  use consts, only: ndim
  use logging_services
  use polygon_type
  implicit none
  private

  public :: curvature_from_patch

contains
  
  real(r8) function curvature_from_patch (interface_reconstructions)

    use analytic_surface_type

    type(polygon), intent(in) :: interface_reconstructions(:)

    integer :: i
    real(r8) :: pts(ndim,3*size(interface_reconstructions))
    type(analytic_surface) :: surf

    ! get points representing each polygon
    do i = 1,size(interface_reconstructions)
      pts(:,3*(i-1)+1:3*(i-1)+3) = polygon_points(interface_reconstructions(i))
    end do

    ! calculate the analytic surface fit and curvature, calculated at the center point
    ! by convention, the first element of interface_reconstructions is the center polygon
    call surf%bestFit (pts)
    curvature_from_patch = surf%curvature(interface_reconstructions(1)%centroid())

    !print '(a,3es14.4,a)', 'c: ',interface_reconstructions(1)%centroid()
    do i = 1,size(pts,2)
      print '(a,3es14.4)', 'x: ',pts(:,i)
    end do
    print *, 'curvature ', curvature_from_patch
    print '(dt)', surf

  end function curvature_from_patch

  ! return a collection of 3 points which fit the polygon
  ! these points should surround the centroid, but not touch the edge of the polygon
  function polygon_points (interface_reconstruction)
    
    type(polygon), intent(in) :: interface_reconstruction
    real(r8) :: polygon_points(ndim,3)
    
    real(r8) :: xcent(ndim), q(ndim,2), d
    
    xcent = interface_reconstruction%centroid()
    !polygon_points(:,1) = 0.5_r8 * (xcent + interface_reconstruction%x(:,1))
    
    ! form an equilateral triangle in the scaled basis
    q = interface_reconstruction%basis()
    d = 0.5_r8
    polygon_points(:,1) = xcent + d*q(:,1)
    polygon_points(:,2) = xcent + d*0.5_r8*(-q(:,1) + sqrt(3.0_r8)*q(:,2))
    polygon_points(:,3) = xcent + d*0.5_r8*(-q(:,1) - sqrt(3.0_r8)*q(:,2))
    
  end function polygon_points

end module interface_patch_type
