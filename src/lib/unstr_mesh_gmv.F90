!!
!! DISTRIBUTED_MESH_GMV
!!
!! Neil N. Carlson <nnc@lanl.gov> 3/30/2006
!! Last revised 11 Aug 2006.
!!
!! This is a quick-n-dirty high-level layer over GMV's gmvwrite C library
!! that provides some procedures for writing a distributed mesh and cell-
!! based fields over that mesh to a GMV-format graphics file.  Although
!! the output is performed only on the IO processor, the routines can be
!! called in parallel, and when distributed data is passed, must be called
!! in parallel.
!!
!! These are the available procedures, presented in the order they should
!! be called:
!!
!! CALL GMV_OPEN (FILE) establishes FILE as the file where the graphics data
!!    will be written.  Any previous contents of the file will be overwritten.
!!
!! CALL GMV_WRITE_DIST_MESH (MESH) writes the distributed mesh MESH to the
!!    graphics file.  The cell block ID information is written as the material
!!    component.  When executing in parallel, mesh partition information is
!!    also written as flag components: cellpart gives the cell partitioning
!!    nodepart the node partitioning.
!!
!! CALL GMV_BEGIN_VARIABLES ([time] [, seq]) prepares the graphics file to
!!    receive data for field variables.  If the optional argument TIME is
!!    present its value is used for the gmv problem time.  If the optional
!!    argument SEQ is present, its value is used for the gmv cycle number.
!!
!! CALL GMV_WRITE_DIST_CELL_VAR (MESH, U, NAME) writes the variable data U to
!!    the graphics file.  U should be a distributed, cell-based field over the
!!    distributed mesh MESH (same as in the earlier call), and NAME is an
!!    arbitrary string used to label the variable; only the first 8 characters
!!    are significant.  This may be called multiple times.
!!
!! CALL GMV_END_VARIABLES () signals that no more variable data will be written.
!!
!! CALL GMV_CLOSE () finalizes the graphics file; nothing more can be written.
!!

#include "f90_assert.fpp"

module unstr_mesh_gmv

  use kinds, only: r8
  use gmvwrite_fortran_binding
  implicit none
  private

  public :: gmv_open, gmv_close, gmv_write_unstr_mesh
  public :: gmv_begin_variables, gmv_end_variables, gmv_write_cell_var
  public :: gmv_write_unstr_mesh_surf, gmv_begin_surfvars, gmv_end_surfvars, gmv_write_surf_var

contains

  subroutine gmv_open (file)
    character(*) :: file
    !! 4-byte integer data and 8-byte real data.
    !call gmvwrite_openfile_ir (file, 4, 8) ! GMV bug with node ids
    call gmvwrite_openfile_ir_ascii (file//C_NULL_CHAR, 4, 8)
  end subroutine gmv_open

  subroutine gmv_close ()
    call gmvwrite_closefile ()
  end subroutine gmv_close

  subroutine gmv_write_unstr_mesh (mesh)

    use unstr_mesh_type
    use string_utilities, only: i_to_c

    type(unstr_mesh), intent(in) :: mesh

    integer :: j
    character(9) :: cell_type, tmpname

    !! Infer cell type from the CNODE array (MESH ought to have a type component!)
    select case (size(mesh%cnode,1))
    case (4)
      cell_type = 'tet'
    case (8)
      cell_type = 'hex'
    case default
      INSIST(.false.)
    end select
    cell_type(9:9) = C_NULL_CHAR
    
    call gmvwrite_node_data (size(mesh%x,dim=2), mesh%x(1,:), mesh%x(2,:), mesh%x(3,:)) ! segfault here
    call gmvwrite_cell_header (size(mesh%cnode,dim=2))
    do j = 1, size(mesh%cnode,dim=2)
      call gmvwrite_cell_type (cell_type, size(mesh%cnode,dim=1), mesh%cnode(:,j))
    end do

    !! Write external mesh node numbers as the nodeids -- GMV uses these for display.
    call gmvwrite_nodeids (mesh%xnode)

    !! Write external mesh cell numbers as the cellids -- GMV uses these for display.
    call gmvwrite_cellids (mesh%xcell)

    !! Write the cell block IDs as the cell material.
    call gmvwrite_material_header (size(mesh%block_id), CELLDATA)
    do j = 1, size(mesh%block_id)
      tmpname = 'block' // i_to_c(mesh%block_id(j))
      tmpname(9:9) = C_NULL_CHAR
      call gmvwrite_material_name (tmpname)
    end do
    call gmvwrite_material_ids (mesh%cblock, CELLDATA)

  end subroutine gmv_write_unstr_mesh

  subroutine gmv_begin_variables (time, seq)
    real(r8), intent(in), optional :: time
    integer, intent(in), optional :: seq
    if (present(time)) call gmvwrite_probtime (time)
    if (present(seq))  call gmvwrite_cycleno (seq)
    call gmvwrite_variable_header
  end subroutine gmv_begin_variables

  subroutine gmv_end_variables
    call gmvwrite_variable_endvars
  end subroutine gmv_end_variables

  subroutine gmv_write_cell_var (mesh, u, name)

    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: mesh
    real(r8),   intent(in) :: u(:)
    character(*), intent(in) :: name

    character(9) :: tmpname

    ASSERT(size(u) == mesh%ncell)

    tmpname = name
    tmpname(9:9) = C_NULL_CHAR
    call gmvwrite_variable_name_data (CELLDATA, tmpname, u)

  end subroutine gmv_write_cell_var

  subroutine gmv_write_unstr_mesh_surf (mesh)

    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: mesh

    integer :: nsurf, j
    integer, allocatable :: snode(:,:)

    nsurf = count(mesh%face_set_mask /= 0)
    if (nsurf == 0) return
    call gmvwrite_surface_header (nsurf)

    do j = 1, mesh%nface
      if (mesh%face_set_mask(j) /= 0) &
          call gmvwrite_surface_data (size(mesh%fnode,1), mesh%fnode(:,j))
    end do

  end subroutine gmv_write_unstr_mesh_surf

  subroutine gmv_begin_surfvars (time, seq)
    real(r8), intent(in), optional :: time
    integer, intent(in), optional :: seq
    if (present(time)) call gmvwrite_probtime (time)
    if (present(seq))  call gmvwrite_cycleno (seq)
    call gmvwrite_surfvars_header
  end subroutine gmv_begin_surfvars

  subroutine gmv_end_surfvars
    call gmvwrite_surfvars_endsvar
  end subroutine gmv_end_surfvars

  subroutine gmv_write_surf_var (mesh, u, name)

    use unstr_mesh_type

    type(unstr_mesh), intent(in) :: mesh
    real(r8), intent(in) :: u(:)
    character(*), intent(in) :: name

    character(9) :: tmpname

    ASSERT(size(u) == mesh%nface)

    tmpname = name
    tmpname(9:9) = C_NULL_CHAR
    call gmvwrite_surfvars_name_data (tmpname, pack(u,mask=(mesh%face_set_mask /= 0)))

  end subroutine gmv_write_surf_var

end module unstr_mesh_gmv
