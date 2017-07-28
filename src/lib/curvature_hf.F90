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
  use logging_services
  implicit none
  private

  public :: curvatureHF, heightFunction

contains

  subroutine heightFunction (curvature_hf, normal_hf, volfrac, normal, mesh, gmesh)

    use timer_tree_type

    real(r8), allocatable, intent(out) :: curvature_hf(:), normal_hf(:,:)
    real(r8), intent(in) :: volfrac(:), normal(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8) :: curvatureHF(size(volfrac))

    integer :: i

    call start_timer ("height function")

    if (allocated(curvature_hf)) deallocate(curvature_hf)
    if (allocated(normal_hf)) deallocate(normal_hf)

    allocate(curvature_hf(mesh%ncell), normal_hf(3,mesh%ncell))

    do i = 1,mesh%ncell
      if (isMixedCell(volfrac(i))) then
        call HFCell(curvature_hf(i), normal_hf(:,i), volfrac, normal(:,i), mesh, gmesh, i)

        ! if (.not.any(gmesh%cneighbor(:,i)<1)) then
        !   print '(3es13.3)', normal_hf(:,i)
        !   print '(3es13.3)', normal(:,i)
        !   print '(es13.3)', norm2(normal_hf(:,i) - normal(:,i))
        !   print *
        ! end if
      else
        curvature_hf(i) = 0
        normal_hf(:,i) = 0
      end if
    end do

    call stop_timer ("height function")

  end subroutine heightFunction

  function curvatureHF (volfrac, normal, mesh, gmesh)

    real(r8), intent(in) :: volfrac(:), normal(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    real(r8) :: curvatureHF(size(volfrac))

    integer :: i
    real(r8) :: throwaway(3)

    do i = 1,mesh%ncell
      if (isMixedCell(volfrac(i))) then
        call HFCell(curvatureHF(i), throwaway, volfrac, normal(:,i), mesh, gmesh, i)
      else
        curvatureHF(i) = 0
      end if
    end do

  end function curvatureHF

  logical function isMixedCell (volfrac)
    use consts, only: cutvof
    real(r8), intent(in) :: volfrac
    isMixedCell = volfrac > cutvof .and. volfrac < 1.0_r8 - cutvof
  end function isMixedCell

  subroutine HFCell (curvature_hf, normal_hf, volfrac, normal, mesh, gmesh, cid)

    use array_utils, only: magnitude, normalize, isZero

    real(r8), intent(out) :: curvature_hf, normal_hf(:)
    real(r8), intent(in) :: volfrac(:), normal(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom), intent(in) :: gmesh
    integer, intent(in) :: cid

    real(r8) :: volfrac_local(3,3,7), hf(3,3), Hx, Hy, Hxx, Hxy, Hyy, dx
    integer :: i_dir, j_dir, k_dir

    dx = 2 * magnitude(gmesh%xc(:,cid) - gmesh%fc(:,mesh%cface(1,cid)))
    k_dir = maxloc(abs(normal), 1)
    i_dir = modulo(k_dir  ,3) + 1
    j_dir = modulo(k_dir+1,3) + 1

    ! get the stencil for the height function,
    ! then calculate the height function from the volume fractions
    volfrac_local = HFStencil(volfrac, mesh, gmesh, cid, i_dir, j_dir, k_dir)
    hf = sum(volfrac_local, dim=3) * dx

    ! calculate height function derivatives
    Hx = (hf(3,2) - hf(1,2)) / (2*dx)
    Hy = (hf(2,3) - hf(2,1)) / (2*dx)
    Hxx = (hf(3,2) - 2*hf(2,2) + hf(1,2)) / dx**2
    Hyy = (hf(2,3) - 2*hf(2,2) + hf(2,1)) / dx**2
    Hxy = (hf(3,3) - hf(1,3) - hf(3,1) + hf(1,1)) / (4*dx*dx)

    ! calculate normal
    normal_hf = - [Hx, Hy,  - sign(1.0_r8, normal(k_dir))] / sqrt(1 + Hx**2 + Hy**2)
    normal_hf([i_dir, j_dir, k_dir]) = normal_hf

    ! calculate curvature
    if (hf(2,2) > 3*dx .and. hf(2,2) < 4*dx) then
      curvature_hf = - & !sign(1.0_r8,interface_normal(dir,c)) * &
          (Hxx + Hyy + Hxx*Hy**2 + Hyy*Hx**2 - 2*Hxy*Hx*Hy) / (1+Hx**2+Hy**2)**1.5_r8
    else
      curvature_hf = 0
    end if

    ! if (.not.any(gmesh%cneighbor(:,cid)<1) &
    !     !.and. i_dir == 2 &
    !     .and. norm2(normal_hf - normal) > 1e-1_r8 &
    !     !.and. norm2(normal_hf + [0.8_r8, 0.6_r8, 0.0_r8]) > 5e-5_r8 &
    !     ) then
    !   print *, i_dir, j_dir, k_dir
    !   print *, hf

    !   print '(3es13.3)', normal
    !   print '(3es13.3)', normal_hf

    !   !print '(es13.3)', norm2(normal_hf + [0.8_r8, 0.6_r8, 0.0_r8])
    !   !print '(3es13.3)', normal_hf + [0.8_r8, 0.6_r8, 0.0_r8]

    !   print *
    !   call LS_fatal ("DEBUGGING")
    ! end if

  end subroutine HFCell

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
