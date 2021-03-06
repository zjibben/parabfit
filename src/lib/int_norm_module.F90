!!
!! INT_NORM_MODULE
!!
!! This module provides routines for calculating
!! interface normals from vof values
!!
!! Zechariah J. Jibben <zjibben@lanl.gov>
!! Last revised April 2016.
!!

module int_norm_module

  use kinds,  only: r8
  use consts, only: alittle,nfc,ndim
  use unstr_mesh_type
  use mesh_geom_type
  implicit none
  private

  public :: interface_normal, interface_normal_cell

contains

  ! This routine calculates the interface normal for each
  ! material by computing the gradient of material volume fractions
  function interface_normal (vof, mesh, gmesh, do_onion_skin)

    use timer_tree_type
    use array_utils, only: first_true_loc,last_true_loc

    real(r8),         intent(in) :: vof(:,:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    logical,          intent(in) :: do_onion_skin
    real(r8)                     :: interface_normal(3,size(vof,dim=1),mesh%ncell)

    integer :: m, n

    call start_timer ("normal calculation")

    ! Calculate the volume fraction gradient for each material
    ! recall interface normal is in opposite direction of gradient
    do m = 1,size(vof, dim=1)
      interface_normal(:,m,:) = -gradient (vof(m,:), mesh, gmesh)
    end do

    !$omp parallel do
    do n = 1,mesh%ncell
      call cell_gradients_to_normals (interface_normal(:,:,n), vof(:,n), do_onion_skin)
    end do
    !$omp end parallel do

    call stop_timer ("normal calculation")

  end function interface_normal

  ! calculate the interface normal in a single cell
  function interface_normal_cell (vof, n, mesh, gmesh, do_onion_skin) result(interface_normal)

    use timer_tree_type
    use array_utils, only: first_true_loc,last_true_loc

    real(r8),         intent(in) :: vof(:,:)
    integer,          intent(in) :: n
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    logical,          intent(in) :: do_onion_skin
    real(r8)                     :: interface_normal(3,size(vof,dim=1))

    integer :: m

    ! Calculate the volume fraction gradient for each material
    ! recall interface normal is in opposite direction of gradient
    do m = 1,size(vof, dim=1)
      interface_normal(:,m) = -gradient_cell (vof(m,:), n, mesh, gmesh)
    end do

    call cell_gradients_to_normals (interface_normal, vof(:,n), do_onion_skin)

  end function interface_normal_cell

  ! takes the collection of vof gradients in a cell, and ensures they are consistent,
  ! performs onion skin (if requested), then normalizes and removes noise.
  subroutine cell_gradients_to_normals (interface_normal, vof, do_onion_skin)

    real(r8),         intent(inout) :: interface_normal(:,:)
    real(r8),         intent(in) :: vof(:)
    logical,          intent(in) :: do_onion_skin

    integer :: m

    call normal_matching_cell (interface_normal, vof, do_onion_skin)

    do m = 1,size(interface_normal, dim=2)
      call grad2norm (interface_normal(:,m))
    end do

  end subroutine cell_gradients_to_normals

  ! when onion skin is requested, performs onion skin on interface normals
  ! sets normals to materials which are not present
  ! when exactly two materials are in a cell, ensures they are consistent
  subroutine normal_matching_cell (interface_normal, vof, do_onion_skin)

    use array_utils, only: first_true_loc,last_true_loc

    real(r8),         intent(inout) :: interface_normal(:,:)
    real(r8),         intent(in) :: vof(:)
    logical,          intent(in) :: do_onion_skin

    integer :: m, nmat, nmat_in_cell
    logical, allocatable :: matl_in_cell(:)

    nmat = size(vof, dim=1)

    ! If there are more than 2 materials in the cell and we are doing onion skin,
    ! sum the gradients in priority order. This is equivalent to calculating the
    ! interface normal for a material composed of the first few materials
    matl_in_cell = vof > 0.0_r8
    nmat_in_cell = count(matl_in_cell)
    if (nmat_in_cell > 2) then
      do m = 2,nmat
        if (vof(m) > alittle .and. do_onion_skin) then
          interface_normal(:,m) = interface_normal(:,m) + interface_normal(:,1)
        else if (vof(m) <= alittle) then
          interface_normal(:,m) = interface_normal(:,1)
        end if
      end do
    else if (nmat_in_cell == 2) then
      ! if there are only two materials in the cell,
      ! ensure the normals are consistent identify the material IDs in the cells
      interface_normal(:,last_true_loc(matl_in_cell)) = &
          -interface_normal(:,first_true_loc(matl_in_cell))
    end if

  end subroutine normal_matching_cell

  ! convert a gradient to a normal by eliminating tiny components,
  ! normalizing, and handling 0-gradients
  !
  ! note 1: The threshold here is higher, following what was done in Truchas.
  !         With the threshold lower, we can get normals like [1.0, 1e-8, 1e-8].
  !         This problem, along with planes barely touching polyhedron
  !         vertices, sometimes results in invalid polyhedra in nested
  !         dissection. -zjibben
  pure subroutine grad2norm (n)

    use array_utils, only: magnitude, eliminateNoise

    real(r8), intent(inout) :: n(:)

    real(r8) :: mag

    call eliminateNoise (n)          ! set tiny components to zero
    mag = magnitude(n)               ! calculate the gradient magnitude
    if (mag > alittle) n = n/mag     ! normalize when it won't blow up
    call eliminateNoise (n, 1e-6_r8) ! set tiny components to zero (see note 1)
    if (all(n == 0.0_r8)) n = 1.0_r8 ! set all components to 1.0 in pure cells

  end subroutine grad2norm

  !   Compute the cell-centered gradient (dPhi_dx,dPhi_dy,dPhi_dz) of a
  !   cell-centered scalar quantity Phi with the green-gauss method:
  !       With a discrete approximation to Gauss's theorem, whereby
  !       the integral of (dPhi_dx,dPhi_dy,dPhi_dz) over the cell
  !       volume is converted to an integral of Phi over the cell
  !       surface area vector.  The area integral is approximated
  !       as a discrete sum over cell faces.  This control volume
  !       formulation is discretely conservative, i.e., adjacent
  !       face contributions will telescope, leaving only boundary
  !       contributions.
  !
  !  note 1: Truchas does this with linear prop. It could also possibly
  !          be done locally just through the 2 cells sharing the face,
  !          using some flux function.
  function gradient (phi, mesh, gmesh)

    use timer_tree_type
    use array_utils, only: eliminateNoise

    real(r8),         intent(in) :: phi(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    real(r8)                     :: gradient(3,size(phi))

    integer  :: f,i
    real(r8) :: Phi_f(mesh%nface)

    call start_timer("gradient")

    ! compute face average of phi
    ! this is calculated through the average of node values at a face
    ! the node values are calculated through the average of neighboring cell values
    ! this could be modified and moved inside the following loop,
    ! although some effort would be duplicated that way
    phi_f = face_avg_global (vertex_avg_global (phi, mesh, gmesh), mesh)

    ! green-gauss method
    ! Loop over faces, accumulating the product
    ! phi_f*Face_Normal for each area vector component

    !$omp parallel do private(f)
    do i = 1,mesh%ncell
      ! Accumulate the dot product
      gradient(:,i) = 0.0_r8
      do f = 1,nfc
        gradient(:,i) = gradient(:,i) + &
            gmesh%outnorm(:,f,i)*mesh%area(mesh%cface(f,i))*phi_f(mesh%cface(f,i))
      end do
      gradient(:,i) = gradient(:,i) / mesh%volume(i) ! normalize by cell volume

      call eliminateNoise (gradient(:,i))
    end do
    !$omp end parallel do

    call stop_timer("gradient")

  end function gradient

  function gradient_cell (phi, i, mesh, gmesh) !result(gradient)

    use timer_tree_type
    use array_utils, only: eliminateNoise

    real(r8),         intent(in) :: phi(:)
    integer,          intent(in) :: i
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    real(r8)                     :: gradient_cell(ndim)

    integer  :: f

    ! green-gauss method
    ! Loop over faces, accumulating the product
    ! Phi_f*Face_Normal for each area vector component

    ! Accumulate the dot product
    gradient_cell = 0.0_r8
    do f = 1,nfc
      gradient_cell = gradient_cell + &
          gmesh%outnorm(:,f,i)*mesh%area(mesh%cface(f,i))*face_avg (phi,mesh%cface(f,i), mesh, gmesh)
    end do

    gradient_cell = gradient_cell / mesh%volume(i) ! Normalize by cell volume
    call eliminateNoise (gradient_cell)

  end function gradient_cell

  real(r8) function face_avg (x_cell, f, mesh, gmesh)

    use consts, only: nvf

    real(r8),         intent(in) :: x_cell(:)
    integer,          intent(in) :: f
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh

    integer  :: v

    ! loop through each vertex on the face,
    ! adding up the average value for the node (found through neighboring cells)
    face_avg = 0.0_r8
    do v = 1,nvf
      face_avg = face_avg + vertex_avg (x_cell, mesh%fnode(v,f), mesh, gmesh)
    end do
    face_avg = face_avg / nvf

  end function face_avg

  function face_avg_global (x_vtx, mesh) result(x_face)

    use consts, only: nvf

    real(r8),         intent(in) :: x_vtx(:)
    type(unstr_mesh), intent(in) :: mesh
    real(r8)                     :: x_face(mesh%nface)

    integer :: f

    ! interpolate vertex values to each face (see note 1)
    !$omp parallel do
    do f = 1,mesh%nface
      x_face(f) = sum(x_vtx(mesh%fnode(:,f))) / nvf
    end do
    !$omp end parallel do

  end function face_avg_global

  ! Compute a vertex-averaged quantity X_vertex from a cell-centered
  ! quantity X_cell.  X_vertex is an inverse-volume-weighted average
  ! of X_cell: X_vertex = SUM(X_cell/Vol)/SUM(1/Vol).
  function vertex_avg_global (x_cell, mesh, gmesh) result(x_vertex)
    real(r8),         intent(in) :: x_cell(:)
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    real(r8)                     :: x_vertex(mesh%nnode)

    integer  :: n

    ! calculate the inverse sum of inverse volumes
    !$omp parallel do
    do n = 1,mesh%nnode
      x_vertex(n) = vertex_avg (x_cell, n, mesh, gmesh)
    end do
    !$omp end parallel do

  end function vertex_avg_global

  real(r8) function vertex_avg (x_cell, n, mesh, gmesh)
    real(r8),         intent(in) :: X_cell(:)
    integer,          intent(in) :: n
    type(unstr_mesh), intent(in) :: mesh
    type(mesh_geom),  intent(in) :: gmesh
    !real(r8)                     :: X_vertex(mesh%nnode)

    integer  :: i,cellid
    real(r8) :: tmp1,tmp2

    ! calculate the inverse sum of inverse volumes
    tmp1 = 0.0_r8; tmp2 = 0.0_r8
    do i = 1,size(gmesh%vcell(:,n))
      cellid = gmesh%vcell(i,n)
      if (cellid>0) then
        tmp1 = tmp1 + X_cell(cellid) / mesh%volume(cellid)
        tmp2 = tmp2 +         1.0_r8 / mesh%volume(cellid)
      end if
    end do
    vertex_avg = tmp1/tmp2
    ! X_vertex(n) = sum( X_cell(gmesh%vcell(:,n)) / mesh%volume(gmesh%vcell(:,n)), mask=(gmesh%vcell(:,n)>0)) &
    !  /        sum(                   1.0_r8 / mesh%volume(gmesh%vcell(:,n)), mask=(gmesh%vcell(:,n)>0))

  end function vertex_avg

end module int_norm_module
