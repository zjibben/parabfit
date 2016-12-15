!!
!! curvature
!!
!! This module defines a function for calculating curvature from a
!! set of interface reconstructions.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! October 2016
!!

!! TODO: * Instead of taking a collection of polygons, this could
!!         be a method inside the surface type. Note this would require
!!         the surface type to work with threading, and eventually
!!         distributed memory.

#include "f90_assert.fpp"

module curvature

  use kinds, only: r8
  use logging_services
  implicit none
  private

  public :: meanCurvature,exnorm

  interface exnorm
    procedure exnorm_r8, exnorm_r8r8
  end interface exnorm

contains

  ! Return the curvature for a cell, given it and its neighbors' cell reconstructions
  ! The argument, by convention, should contain an array of 27 polygons, most of which
  ! will be null. The center polygon, or element 14, must not be null.

  ! This particular version tries to calculate the curvature directly from the local vof
  real(r8) function meanCurvature (intnorm, vof, face_area, cell_vol, i, mesh, gmesh)

    use consts, only: nfc,cutvof,ndim
    use unstr_mesh_type
    use mesh_geom_type
    use differential_operators, only: divergence !, faceGradient
    use array_utils, only: normalize, interpolate, magnitude, meanArithmetic
    !use polygon_type

    real(r8),         intent(in) :: intnorm(:,:), vof(:), face_area(:), cell_vol
    integer,          intent(in) :: i
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh

    real(r8) :: intnorm_out(nfc), tmp(ndim), n(ndim,2)
    integer :: i_ngbr,f

    ! ASSERT(size(reconstruction)==27)
    ! ASSERT(reconstruction(14)/=NULL)

    ! print '(a,3f10.4,a)', '    (',intnorm(:,i),')'
    ! print '(a,3f10.4,a)', '    (',-normalize(gmesh%xc(:,i)-[0.5_r8,0.5_r8,gmesh%xc(3,i)]),')'

    ! calculate the outward component of the gradient of the vof on every face
    intnorm_out = 0.0_r8
    do f = 1,nfc
      i_ngbr = gmesh%cneighbor(f,i)

      ! intnorm_out(f) = dot_product(&
      !     sum(exnorm(mesh%x(:,mesh%fnode(:,mesh%cface(f,i)))),dim=2)/real(4,r8),&
      !     gmesh%outnorm(:,f,i))
      ! if (i_ngbr > 0) then
      !   ! intnorm_out(f) = dot_product(&
      !   !     normalize(faceGradient (vof([i, i_ngbr]), gmesh%xc(:,[i, i_ngbr]))), &
      !   !     gmesh%outnorm(:,f,i))
      !   ! if (.not.all(intnorm(:,i_ngbr)==1.0_r8) &
      !   !     .and. vof(i_ngbr) > cutvof .and. vof(i_ngbr) < 1.0_r8-cutvof) then! .and. &
      !   !dot_product(intnorm(:,i),gmesh%outnorm(:,f,i)) < 0.8_r8) then
        
      !   n(:,1) = -normalize(gmesh%xc(:,i)-[0.5_r8,0.5_r8,gmesh%xc(3,i)])
      !   n(:,2) = -normalize(gmesh%xc(:,i_ngbr)-[0.5_r8,0.5_r8,gmesh%xc(3,i_ngbr)])
      !   tmp = normalize(interpolate(n, gmesh%xc(:,[i,i_ngbr]), gmesh%fc(:,mesh%cface(f,i))))

          intnorm_out(f) = dot_product( &
              !normalize(interpolate(intnorm(:,[i,i_ngbr]), gmesh%xc(:,[i,i_ngbr]), gmesh%fc(:,mesh%cface(f,i)))), &
              -normalize(gmesh%fc(:,mesh%cface(f,i))-[0.5_r8,0.5_r8,gmesh%fc(3,mesh%cface(f,i))]), &
              !tmp, &
              gmesh%outnorm(:,f,i))
      !     ! else ! boundary face gradient already calculated
      !     !   intnorm_out(f) = intnorm_out_bndry(n,f)
      !   ! else
      !   !   intnorm_out(f) = dot_product(intnorm(:,i), gmesh%outnorm(:,f,i))
      !   ! end if

      !     ! !print '(3(a,3es12.2,a))', &
      !     ! print '(2(a,3es12.2,a),2es12.2)', &
      !     !   '    (',normalize(interpolate(intnorm(:,[i,i_ngbr]), gmesh%xc(:,[i,i_ngbr]), gmesh%fc(:,mesh%cface(f,i)))),')',&
      !     !   '    (',-normalize(gmesh%fc(:,mesh%cface(f,i))-[0.5_r8,0.5_r8,gmesh%fc(3,mesh%cface(f,i))]),')',&
      !     !   ! '    (',normalize(interpolate(intnorm(:,[i,i_ngbr]), gmesh%xc(:,[i,i_ngbr]), gmesh%fc(:,mesh%cface(f,i)))) &
      !     !   ! + normalize(gmesh%fc(:,mesh%cface(f,i))-[0.5_r8,0.5_r8,gmesh%fc(3,mesh%cface(f,i))]),')'
      !     !   dot_product(normalize(interpolate(intnorm(:,[i,i_ngbr]), gmesh%xc(:,[i,i_ngbr]), gmesh%fc(:,mesh%cface(f,i)))),&
      !     !   gmesh%outnorm(:,f,i)),&
      !     !   dot_product(-normalize(gmesh%fc(:,mesh%cface(f,i))-[0.5_r8,0.5_r8,gmesh%fc(3,mesh%cface(f,i))]),&
      !     !   gmesh%outnorm(:,f,i))

      !     ! tmp = -normalize(gmesh%xc(:,i)-[0.5_r8,0.5_r8,gmesh%xc(3,i)])
      !     ! print '((a,3es12.2,a), 2es12.2)', '    (',intnorm(:,i),')', acos(intnorm(2,i)),acos(tmp(2))
      !     ! tmp = -normalize(gmesh%xc(:,i_ngbr)-[0.5_r8,0.5_r8,gmesh%xc(3,i_ngbr)])
      !     ! print '(2(a,3es12.2,a),2es12.2)', '    (',intnorm(:,i_ngbr),')', &
      !     !     '    (',-normalize(gmesh%xc(:,i_ngbr)-[0.5_r8,0.5_r8,gmesh%xc(3,i_ngbr)]),')', &
      !     !     acos(intnorm(2,i_ngbr)),acos(tmp(2))

      !     ! print *

      !   ! print '(i3,a,3f10.3,a,f10.3,a,3f10.3,a,2l3,2es14.3)', f, &
      !   !   '    (',intnorm(:,i_ngbr),')', intnorm_out(f),&
      !   !   '    (',gmesh%outnorm(:,f,i),')', &
      !   !   .not.all(intnorm(:,i_ngbr)==1.0_r8) .and. vof(i_ngbr) > cutvof .and. vof(i_ngbr) < 1.0_r8-cutvof, &
      !   !   vof(i_ngbr) > cutvof .and. vof(i_ngbr) < 1.0_r8-cutvof, &
      !   !   vof(i_ngbr), dot_product(intnorm(:,i),gmesh%outnorm(:,f,i))
      !   else
      !     intnorm_out(f) = dot_product( &
      !         -normalize(gmesh%fc(:,mesh%cface(f,i))-[0.5_r8,0.5_r8,gmesh%fc(3,mesh%cface(f,i))]), &
      !         gmesh%outnorm(:,f,i))
      !     !intnorm_out(f) = dot_product(intnorm(:,i), gmesh%outnorm(:,f,i))
      ! end if
    end do

    !print '(6es10.2)', intnorm_out
    meanCurvature = -divergence (intnorm_out, face_area, cell_vol)

  end function meanCurvature

  function exnorm_r8(x) result(n)

    use consts, only: ndim
    use array_utils, only: normalize

    real(r8), intent(in) :: x(:)
    real(r8) :: n(ndim)

    n = -normalize(x - [0.5_r8,0.5_r8,x(3)])

  end function exnorm_r8

  function exnorm_r8r8(x) result(n)
    
    real(r8), intent(in) :: x(:,:)
    real(r8) :: n(size(x,1),size(x,2))

    integer :: i

    do i = 1,size(x,2)
      n(:,i) = exnorm_r8(x(:,i))
    end do

  end function exnorm_r8r8

end module curvature
