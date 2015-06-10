!!
!! MFD_DISC_TYPE
!!
!! This module defines a derived type and associated procedures that provide
!! some core functionality for implementing mimetic finite difference (MFD)
!! discretizations of differential operators such as the diffusion operator.
!!
!! Neil Carlson <nnc@lanl.gov>
!! Adapted for Fortran 2008, February 2014
!!
!! PROGRAMMING INTERFACE
!!
!!  This module defines the derived type MFD_DISC which composes a reference
!!  to the computational mesh with additional persistent mesh-dependent data
!!  that is specfic to MFD.  This type is intended to be used by trusted code
!!  and so the data components are public.  However, the components should
!!  be regarded as read-only, just as the mesh, because it is also expected
!!  that an instance of this type will be referenced by multiple higher-level
!!  objects.  The available components are:
!!
!!    MESH -- a pointer to the underlying computational mesh
!!    MINV -- the inverse of the local flux mass matrices.  MINV(:,j) is the
!!      inverse of the (symmetric) local flux mass matrix on cell j in upper
!!      packed matrix format.
!!
!!  The MFD_DISC type has the following type bound procedures.
!!
!!  INIT(MESH) initializes the object with MESH as the underlying computational
!!    mesh.  The object holds a reference to MESH, and so the actual argument
!!    must either be a pointer or have the target attribute, and must persist
!!    for the lifetime of the object.  After this call the MINV component of
!!    the object is defined.
!!
!!  APPLY_DIFF(COEF, UCELL, UFACE, RCELL, RFACE) applies the local MFD diffusion
!!    operator to a vector.  The vector is partitioned into cell and face-based
!!    parts in the UCELL and UFACE vectors, and the like-partitioned result is
!!    returned in the RCELL and RFACE vectors.  COEF is the cell-based vector of
!!    diffusion coefficients.
!!
!!  COMPUTE_CELL_GRAD (UFACE, GRAD) computes an approximation to the average
!!    gradient of a field u on each cell given average values of u on the faces
!!    of the mesh.  The face values are passed in the array UFACE and the
!!    computed gradient in the array GRAD; GRAD(:,j) is the gradient on cell j.
!!    Note that this ignores the average cell values of u which are also
!!    available (typically).  The results are really only intended to be used
!!    for output, and may not be at all suitable for use in a discretization
!!    scheme.
!!

#include "f90_assert.fpp"

module mfd_disc_type

  use kinds, only: r8
  use unstr_mesh_type
  use unstr_mesh_partition_type
  implicit none
  private

  type, public :: mfd_disc
    ! The next two fields are reference only - we do not own them
    type(unstr_mesh), pointer :: mesh => null() 
    type(unstr_mesh_partition), pointer :: partition => null()
    real(r8), allocatable :: minv(:,:)
  contains
    procedure :: init => mfd_disc_init
    procedure :: apply_diff => mfd_disc_apply_diff
    procedure :: compute_cell_gradient => mfd_disc_compute_cell_grad
  end type mfd_disc

  !! Private type for internal use.
  type :: mfd_hex
    real(r8) :: volume
    real(r8) :: corner_volumes(8)
    real(r8) :: face_normals(3,6)
  contains
    procedure :: init => mfd_hex_init
    procedure :: compute_flux_matrix => mfd_hex_compute_flux_matrix
  end type mfd_hex

contains

  subroutine mfd_disc_init (this, mesh, partition)

    class(mfd_disc), intent(out) :: this
    type(unstr_mesh), intent(in), target :: mesh
    type(unstr_mesh_partition), intent(in), target :: partition

    integer :: j
    type(mfd_hex) :: tmp

    !! Keep a reference to the mesh
    this%mesh => mesh

    !! Keep a reference to the partition
    this%partition => partition

    allocate(this%minv(21,mesh%ncell))
    do j = 1, mesh%ncell
      call tmp%init (mesh%x(:,mesh%cnode(:,j)))
      call tmp%compute_flux_matrix (1.0_r8, this%minv(:,j), invert=.true.)
    end do

  end subroutine mfd_disc_init

  subroutine mfd_disc_apply_diff (this, coef, ucell, uface, rcell, rface)

    use upper_packed_matrix, only: upm_matmul
    use timer_tree_type

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in)  :: coef(:)
    real(r8), intent(in)  :: ucell(:), uface(:)
    real(r8), intent(out) :: rcell(:), rface(:)

    integer :: j,st,en

    integer :: p
    real(r8), allocatable :: rface_exp(:)

    ASSERT(size(coef) == this%mesh%ncell)
    ASSERT(size(ucell) == size(coef))
    ASSERT(size(rcell) == size(ucell))
    ASSERT(size(uface) == this%mesh%nface)
    ASSERT(size(rface) == size(uface))

    rface = 0.0_r8

    ! Scatter rface to the partitioned data in rface_expanded.
    call this%partition%scatter_face(rface, rface_exp)

    ! Dispatch a parallel loop over the partitions.
    call start_timer('apply_diff')
    !$OMP parallel private(p)
    !$OMP do
    do p = 1,this%partition%num_part()
       call update_section(p)
    end do
    !$OMP end do
    !$OMP end parallel
    call stop_timer('apply_diff')

    ! Gather the partitioned data in rface_expanded back to rface
    call this%partition%gather_sum_face(rface_exp, rface)

  contains 

    subroutine update_section(p)
      integer, intent(in) :: p
      integer :: st,en,j,offset
      real(r8) :: flux(size(this%mesh%cface,dim=1))

      st=this%partition%part_start(p)
      en=this%partition%part_end(p) 
      offset = st-1
      
      do j = st,en
         flux = coef(j) * upm_matmul(this%minv(:,j), ucell(j) - uface(this%mesh%cface(:,j)))
         rface_exp(this%partition%cface(:,j-offset,p)) = &
              rface_exp(this%partition%cface(:,j-offset,p)) - flux
         rcell(j) = sum(flux)
      end do
    end subroutine update_section

  end subroutine mfd_disc_apply_diff

 !! This procedure computes a cell-wise average gradient of a scalar field
 !! given its values at the mesh faces.  The average gradient for a cell C
 !! is computed using the following formula:
 !!
 !!    1/V \int_C \nabla u dx = 1/V \int_{\partial C} u n\dot dS
 !!                           \approx 1/V \sum_{f\in C} u_f n_f
 !!
 !! where V is the volume of C, u_f is the (average) value of u on face f,
 !! and n_f is the (average) outward normal to face f with magnitude equal
 !! to the area of the face.

  subroutine mfd_disc_compute_cell_grad (this, uface, grad)

    class(mfd_disc), intent(in) :: this
    real(r8), intent(in) :: uface(:)
    real(r8), intent(out) :: grad(:,:)

    integer :: j
    type(mfd_hex) :: tmp

    INSIST(size(grad,1) == 3)
    INSIST(size(uface) == this%mesh%nface)
    INSIST(size(grad,2) <= this%mesh%ncell)

    do j = 1, size(grad,2)
      call tmp%init (this%mesh%x(:,this%mesh%cnode(:,j)))
      grad(:,j) = matmul(tmp%face_normals, uface(this%mesh%cface(:,j))) / tmp%volume
    end do

  end subroutine mfd_disc_compute_cell_grad

  subroutine mfd_hex_init (this, vertices)
    use cell_geometry, only: eval_hex_volumes, hex_face_normals
    class(mfd_hex), intent(out) :: this
    real(r8), intent(in) :: vertices(:,:)
    ASSERT(size(vertices,1) == 3 .and. size(vertices,2) == 8)
    call eval_hex_volumes (vertices, this%volume, this%corner_volumes)
    this%face_normals = hex_face_normals(vertices)
  end subroutine mfd_hex_init

  subroutine mfd_hex_compute_flux_matrix (this, coef, matrix, invert)

    use cell_topology, only: HEX8_VERT_FACE
    use upper_packed_matrix, only: upm_invert

    class(mfd_hex), intent(in) :: this
    real(r8), intent(in) :: coef
    real(r8), intent(out) :: matrix(:)
    logical, intent(in), optional :: invert

    integer :: c, i, j, ii, jj, loc
    real(r8) :: s, cwt(8), Nc(3,3), Mc(3,3)

    ASSERT(size(matrix) == 21)

    matrix = 0.0_r8
    cwt = this%corner_volumes / sum(this%corner_volumes)
    do c = 1, 8
      Nc = this%face_normals(:,HEX8_VERT_FACE(:,c))
      call invert_sym_3x3(matmul(transpose(Nc),Nc), Mc)
      !! Scatter the corner matrix into the full cell flux matrix.
      !! It is essential that HEX8_VERT_FACE(:,c) is an increasing sequence of indices.
      s = this%volume * cwt(c) / coef
      do j = 1, 3
        jj = HEX8_VERT_FACE(j,c)
        do i = 1, j
          ii = HEX8_VERT_FACE(i,c)
          loc = ii + jj*(jj - 1)/2
          matrix(loc) = matrix(loc) + s*Mc(i,j)
        end do
      end do
    end do

    if (present(invert)) then
      if (invert) call upm_invert (matrix)
    end if

  end subroutine mfd_hex_compute_flux_matrix

 !! Direct inversion of a 3x3 symmetrix matrix using the formula that the
 !! inverse equals the transponse of the matrix of cofactors divided by
 !! the determinant.

  subroutine invert_sym_3x3 (A, Ainv)
    real(r8), intent(in)  :: A(3,3)
    real(r8), intent(out) :: Ainv(3,3)
    !! Transpose of the matrix of cofactors ...
    Ainv(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
    Ainv(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
    Ainv(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
    Ainv(1,2) = Ainv(2,1)
    Ainv(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
    Ainv(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
    Ainv(1,3) = Ainv(3,1)
    Ainv(2,3) = Ainv(3,2)
    Ainv(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    !! and scale by the determinant to get the inverse.
    Ainv = (1.0_r8/(A(1,1)*Ainv(1,1) + A(2,1)*Ainv(2,1) + A(3,1)*Ainv(3,1))) * Ainv
  end subroutine invert_sym_3x3

end module mfd_disc_type
