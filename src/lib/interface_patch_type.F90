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

  public :: curvature_from_patch, normal_from_patch

contains

  real(r8) function curvature_from_patch(interface_reconstructions, weight_scale, normal, vof, &
      mesh, gmesh, cid, centroid, verbose)

    use analytic_surface_type
    use paraboloid_type
    use timer_tree_type
    use locate_plane_nd_module
    use polyhedron_type
    use unstr_mesh_type
    use mesh_geom_type
    !use hex_types, only: hex_f, hex_e

    type(polygon_box), intent(in) :: interface_reconstructions(:)
    real(r8), intent(in) :: weight_scale
    real(r8), intent(in) :: normal(:), vof(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    integer, intent(in) :: cid
    real(r8), intent(out), optional :: centroid(:)
    logical, intent(in), optional :: verbose

    integer :: i
    real(r8) :: pts(ndim,3*size(interface_reconstructions)), wgt(3*size(interface_reconstructions))
    real(r8) :: xc(ndim), area, a, d(3)
    !type(analytic_surface) :: surf
    type(paraboloid) :: surf
    logical :: verboseh

    integer :: in, npts, ierr, j
    type(polygon_box) :: patch_polygons
    type(polygon) :: interface_polygon
    type(polyhedron) :: cell

    if (present(verbose)) then
      verboseh = verbose
    else
      verboseh = .false.
    end if

    call start_timer ("fit curvature")

    ASSERT(size(interface_reconstructions(1)%elements) > 0)

    ! npts = 1
    ! pts(:,1) = interface_reconstructions(1)%centroid2()
    ! do in = 1,gmesh%caneighbor(cid)%n_elements
    !   i = gmesh%caneighbor(cid)%elements(in)

    !   if (vof(i) < 1-1e-3_r8 .and. vof(i) > 1e-3_r8) then
    !     npts = npts + 1

    !     call cell%init (ierr, mesh%x(:,mesh%cnode(:,i)), hex_f, hex_e, mesh%volume(i), &
    !         gmesh%outnorm(:,:,i))
    !     if (ierr /= 0) call LS_fatal ('could not initialize cell')

    !     interface_polygon = cell%intersection_verts(&
    !         locate_plane_nd (cell, normal, vof(i)*mesh%volume(i), mesh%volume(i)))
    !     pts(:,npts) = interface_polygon%centroid2()
    !   end if
    ! end do
    ! wgt = 1

    ! npts = size(interface_reconstructions)
    ! do i = 1,npts
    !   pts(:,i) = interface_reconstructions(i)%centroid2()
    !   wgt(i) = interface_reconstructions(i)%area() ** weight_scale
    !   !wgt(i) = 1 / (1 + norm2(pts(:,i) - pts(:,1))) ** weight_scale
    ! end do
    ! call surf%bestFit(pts(:,1:npts), wgt(1:npts), normal)

    ! ! get points representing each polygon
    ! do i = 1,size(interface_reconstructions)
    !   pts(:,3*(i-1)+1:3*(i-1)+3) = polygon_points(interface_reconstructions(i))
    !   wgt(3*(i-1)+1:3*(i-1)+3) = interface_reconstructions(i)%area() ** weight_scale
    ! end do

    ! ! calculate the analytic surface fit and curvature, calculated at the center point
    ! ! by convention, the first element of interface_reconstructions is the center polygon
    ! !call surf%bestOneSheetFit (pts)
    ! !call surf%bestParaboloidFit (pts)
    ! call surf%bestFit (pts, wgt, normal)

    xc = interface_reconstructions(1)%centroid()

    ! patch_polygons = flat_polygon_box(interface_reconstructions)
    ! call surf%volumetricFit(patch_polygons%elements)
    call surf%volumetricFit(interface_reconstructions)
    curvature_from_patch = surf%curvature(xc)

    ! ! height function style curvature point ############
    ! d = 0
    ! d(maxloc(abs(normal),1)) = 1
    ! xc = surf%point_along_line(gmesh%xc(:,cid), d)
    ! !xc = gmesh%xc(:,cid)

    ! ! call cell%init(ierr, cid, mesh, gmesh)
    ! ! if (cell%is_inside(xc)) then
    !   !xc = interface_reconstructions(1)%centroid()
    !   curvature_from_patch = surf%curvature(xc)
    ! ! else
    ! !   curvature_from_patch = 0
    ! ! end if
    ! ! ############################################

    if (present(centroid)) centroid = surf%point_on_surface(xc)

    if (verboseh) then
      ! do i = 1,size(interface_reconstructions)
      !   !print '(a,3es20.10)', 'x: ', pts(:,i)
      !   print '(3(a,es20.10),a)', '[',pts(1,i),',',pts(2,i),',',pts(3,i),'],'
      ! end do
      ! print *
      ! do i = 1,size(interface_reconstructions)
      !   print '(a,es20.10)', 'd: ', abs(norm2(pts(:2,i)) - 0.35_r8)
      ! end do

      ! print *
      ! do i = 1,size(interface_reconstructions)
      !   !print '(a,3es20.10)', 'x: ', pts(:,i)
      !   print '(3es20.10)', interface_reconstructions(i)%centroid2()
      ! end do

      print *, 'RECONSTRUCTION POLYGONS'
      do i = 1,size(interface_reconstructions)
        do j = 1,interface_reconstructions(i)%n_elements
          do in = 1,size(interface_reconstructions(i)%elements(j)%x, dim=2)
            print '(3es20.10)', &
                matmul(surf%rot, interface_reconstructions(i)%elements(j)%x(:,in) - surf%offset)
          end do
          print *
        end do
      end do

      print '(dt)', surf

      print *, 'RECONSTRUCTION POLYGON NORMALS'
      do i = 1,size(interface_reconstructions)
        print *, i, interface_reconstructions(i)%n_elements
        do j = 1,interface_reconstructions(i)%n_elements
          print '(3es20.10)', interface_reconstructions(i)%elements(j)%norm
        end do
        print *
      end do

      print '(a,3es18.8)', 'n:  ', normal
      print '(a,3es18.8)', 'n2: ', surf%normal(xc)
      print '(a,2f10.4)', 'curvature0 ', curvature_from_patch, surf%curvatureQdrtc(xc)

      ! call interface_reconstructions(1)%print_data()

      ! print *
      ! call surf%volumetricFit(interface_reconstructions)
      ! print '(a,f10.4, es20.10)', 'curvature2 ', &
      !     surf%curvature(interface_reconstructions(1)%centroid()), &
      !     abs(surf%curvature(interface_reconstructions(1)%centroid()) + 1/0.35_r8) * 0.35_r8
    end if

    call stop_timer ("fit curvature")

  end function curvature_from_patch

  function normal_from_patch (interface_reconstructions, weight_scale, normal, verbose)

    use analytic_surface_type
    use paraboloid_type
    use timer_tree_type

    type(polygon), intent(in) :: interface_reconstructions(:)
    real(r8), intent(in) :: weight_scale
    real(r8), intent(in) :: normal(:)
    logical, intent(in), optional :: verbose
    real(r8) :: normal_from_patch(3)

    integer :: i
    real(r8) :: pts(ndim,3*size(interface_reconstructions)), wgt(3*size(interface_reconstructions))
    !type(analytic_surface) :: surf
    type(paraboloid) :: surf
    logical :: verboseh

    if (present(verbose)) then
      verboseh = verbose
    else
      verboseh = .false.
    end if

    call start_timer ("fit normals")
    call LS_fatal ("currently turned off normal_from_patch")

    ! do i = 1,size(interface_reconstructions)
    !   pts(:,i) = interface_reconstructions(i)%centroid2()
    !   wgt(i) = interface_reconstructions(i)%area() ** weight_scale
    ! end do
    ! call surf%bestFit (pts(:,1:size(interface_reconstructions)), &
    !     wgt(1:size(interface_reconstructions)), normal)

    !call surf%volumetricFit(interface_reconstructions)

    !normal_from_patch = surf%normal(interface_reconstructions(1)%centroid())
    normal_from_patch = surf%normal_average(interface_reconstructions(1))

    call stop_timer ("fit normals")

  end function normal_from_patch

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
