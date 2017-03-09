!!
!! curvature_hf
!!
!! This module provides a function returning curvature, given volume fractions.
!! It uses the height function method, assuming a regular mesh.
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! February 2017
!!

module curvature_hf

  use kinds, only: r8
  use unstr_mesh_type
  use mesh_geom_type
  implicit none
  private

  public :: curvatureHF

contains

  function curvatureHF (volfrac, normal, mesh, gmesh)

    real(r8), intent(in) :: volfrac(:), normal(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8) :: curvatureHF(size(volfrac))

    integer :: i

    do i = 1,mesh%ncell
      if (isMixedCell(volfrac(i))) then
        curvatureHF(i) = curvatureHFCell (volfrac, normal(:,i), mesh, gmesh, i)
      else
        curvatureHF(i) = 0.0_r8
      end if
    end do

  end function curvatureHF

  logical function isMixedCell (volfrac)
    use consts, only: cutvof
    real(r8), intent(in) :: volfrac
    isMixedCell = volfrac > cutvof .and. volfrac < 1.0_r8 - cutvof
  end function isMixedCell

  real(r8) function curvatureHFCell (volfrac, normal, mesh, gmesh, cid)

    use array_utils, only: magnitude, normalize, isZero

    real(r8), intent(in) :: volfrac(:), normal(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    integer, intent(in) :: cid

    real(r8) :: volfrac_local(3,3,7), hf(3,3), Hx, Hy, Hxx, Hxy, Hyy, dx
    integer :: i_dir, j_dir, k_dir

    dx = 2.0_r8 * magnitude(gmesh%xc(:,cid) - gmesh%fc(:,mesh%cface(1,cid)))
    k_dir = maxloc(abs(normal), 1)
    i_dir = modulo(k_dir  ,3) + 1
    j_dir = modulo(k_dir+1,3) + 1

    ! get the stencil for the height function,
    ! then calculate the height function from the volume fractions
    volfrac_local = HFStencil(volfrac, mesh, gmesh, cid, i_dir, j_dir, k_dir)
    hf = sum(volfrac_local, dim=3) * dx

    ! calculate height function derivatives
    Hx = (hf(3,2) - hf(1,2)) / (2.0_r8*dx)
    Hy = (hf(2,3) - hf(2,1)) / (2.0_r8*dx)
    Hxx = (hf(3,2) - 2.0_r8*hf(2,2) + hf(1,2)) / dx**2
    Hyy = (hf(2,3) - 2.0_r8*hf(2,2) + hf(2,1)) / dx**2
    Hxy = (hf(3,3) - hf(1,3) - hf(3,1) + hf(1,1)) / (4.0_r8*dx*dx)

    ! calculate curvature
    if (hf(2,2) > 3.0_r8*dx .and. hf(2,2) < 4.0_r8*dx) then
      curvatureHFCell = - & !sign(1.0_r8,interface_normal(dir,c)) * &
          (Hxx + Hyy + Hxx*Hy**2 + Hyy*Hx**2 - 2.0_r8*Hxy*Hx*Hy) / (1.0_r8+Hx**2+Hy**2)**1.5_r8
    else
      curvatureHFCell = 0.0_r8
    end if

  end function curvatureHFCell

  function HFStencil (volfrac, mesh, gmesh, cid, i_dir, j_dir, k_dir)

    real(r8), intent(in) :: volfrac(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    integer, intent(in) :: cid, i_dir, j_dir, k_dir
    real(r8) :: HFStencil(3,3,7)

    integer :: i,j,k, ic,jc,kc

    HFStencil = 0.0_r8

    do i = 1,3
      ic = meshNeighborID(cid, i_dir, i-2, mesh, gmesh)
      if (ic < 1) cycle

      do j = 1,3
        jc = meshNeighborID(ic, j_dir, j-2, mesh, gmesh)
        if (jc < 1) cycle

        HFStencil(i,j,4) = volfrac(jc)

        kc = jc
        do k = 5,7
          kc = meshNeighborID(kc, k_dir, 1, mesh, gmesh)
          if (kc < 1) exit
          HFStencil(i,j,k) = volfrac(kc)
        end do

        kc = jc
        do k = 3,1,-1
          kc = meshNeighborID(kc, k_dir, -1, mesh, gmesh)
          if (kc < 1) exit
          HFStencil(i,j,k) = volfrac(kc)
        end do

      end do
    end do

  end function HFStencil

  ! return the neighbor cell id in direction dir on side side
  integer function meshNeighborID (cid, dir, side, mesh, gmesh)

    integer, intent(in) :: cid, dir, side
    type(unstr_mesh), intent(in) :: mesh ! TODO: may not be needed. remove.
    type(mesh_geom), intent(in) :: gmesh

    integer, parameter :: face_id(2,3) = reshape([&
        3,4, & ! x-+
        2,1, & ! y-+
        5,6],& ! z-+
        shape(face_id))

    if (side==-1) then
      meshNeighborID = gmesh%cneighbor(face_id(1,dir), cid)
    else if (side==0) then
      meshNeighborID = cid
    else if (side==1) then
      meshNeighborID = gmesh%cneighbor(face_id(2,dir), cid)
    end if

  end function meshNeighborID

end module curvature_hf
