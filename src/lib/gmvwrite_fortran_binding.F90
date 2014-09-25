!!
!! GMVWRITE_FORTRAN_BINDING
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! May 2014
!!
!! NOTE 1: For binary format file, function writes 8 char of the variable.
!! This means the actual character argument must have length at least 8, and
!! no null termination is needed.  For ascii format file, function writes the
!! passed string, which must be null terminated.  Uncertain whether the string
!! must have a specific length, or whether there is a maximum length.
!!
!! NOTE 2: Argument is really void* and is cast to the appropriate type
!! (int/long or float/double) by the function according to the file types
!! specified when the file was opened.  Here the type (int or double) is
!! hardwired for our usage; THIS NEEDS TO BE GENERALIZED.
!!

module gmvwrite_fortran_binding

  use,intrinsic :: iso_c_binding
  implicit none
  public

  integer(c_int), parameter :: CELLDATA = 0, NODEDATA = 1

  character(3), parameter :: TRICELLTYPE  = 'tri'
  character(4), parameter :: QUADCELLTYPE = 'quad'
  character(3), parameter :: TETCELLTYPE  = 'tet'
  character(3), parameter :: HEXCELLTYPE  = 'hex'

  interface
    !void gmvwrite_openfile_ir(char filenam[], int isize, int rsize)
    subroutine gmvwrite_openfile_ir (filenam, isize, rsize) bind(c)
      import c_char, c_int
      character(kind=c_char), intent(in) :: filenam(*)  ! null-terminated string
      integer(c_int), intent(in), value :: isize, rsize
    end subroutine

    !void gmvwrite_openfile_ir_ascii(char filenam[], int isize, int rsize)
    subroutine gmvwrite_openfile_ir_ascii (filenam, isize, rsize) bind(c)
      import c_char, c_int
      character(kind=c_char), intent(in) :: filenam(*)  ! null-terminated string
      integer(c_int), intent(in), value :: isize, rsize
    end subroutine

    !void gmvwrite_closefile(void)
    subroutine gmvwrite_closefile () bind(c)
    end subroutine

    !void gmvwrite_node_data(void *nndes, void *x, void *y, void *z)
    subroutine gmvwrite_node_data (nndes, x, y, z) bind(c)
      import c_int, c_double
      integer(c_int), intent(in) :: nndes  ! SEE NOTE 2
      real(c_double), intent(in) :: x(*), y(*), z(*)  ! SEE NOTE 2
    end subroutine

    !void gmvwrite_cell_header(void *ncells)
    subroutine gmvwrite_cell_header (ncells) bind(c)
      import c_int
      integer(c_int), intent(in) :: ncells ! SEE NOTE 2
    end subroutine

    !void gmvwrite_cell_type(char cell_type[], int nverts, void *nodes)
    subroutine gmvwrite_cell_type (cell_type, nverts, nodes) bind(c)
      import c_char, c_int
      character(kind=c_char), intent(in) :: cell_type(*)  ! SEE NOTE 1
      integer(c_int), intent(in), value :: nverts
      integer(c_int), intent(in) :: nodes(*)  ! SEE NOTE 2
    end subroutine

    !void gmvwrite_material_header(int nmats, int data_type)
    subroutine gmvwrite_material_header (nmats, data_type) bind(c)
      import c_int
      integer(c_int), intent(in), value :: nmats, data_type
    end subroutine

    !void gmvwrite_material_name(char matname[])
    subroutine gmvwrite_material_name (matname) bind(c)
      import c_char
      character(kind=c_char), intent(in) :: matname(*)  ! SEE NOTE 1
    end subroutine

    !void gmvwrite_material_ids(int matids[], int data_type)
    subroutine gmvwrite_material_ids (matids, data_type) bind(c)
      import c_int
      integer(c_int), intent(in) :: matids(*)
      integer(c_int), intent(in), value :: data_type
    end subroutine

    !void gmvwrite_variable_header(void)
    subroutine gmvwrite_variable_header () bind(c)
    end subroutine

    !void gmvwrite_variable_name_data(int data_type, char varname[], void *vids)
    subroutine gmvwrite_variable_name_data (data_type, varname, vids) bind(c)
      import c_int, c_char, c_double
      integer(c_int), intent(in), value :: data_type
      character(kind=c_char), intent(in) :: varname(*)  ! SEE NOTE 1
      real(c_double), intent(in) :: vids(*) ! SEE NOTE 2
    end subroutine

    !void gmvwrite_variable_endvars(void)
    subroutine gmvwrite_variable_endvars () bind(c)
    end subroutine

    !void gmvwrite_flag_header(void)
    subroutine gmvwrite_flag_header () bind(c)
    end subroutine

    !void gmvwrite_flag_name(char flagname[], int numtypes, int data_type)
    subroutine gmvwrite_flag_name (flagname, numtypes, data_type) bind(c)
      import c_char, c_int
      character(kind=c_char), intent(in) :: flagname(*) ! SEE NOTE 1
      integer(c_int), intent(in), value :: numtypes, data_type
    end subroutine

    !void gmvwrite_flag_subname(char subname[])
    subroutine gmvwrite_flag_subname (subname) bind(c)
      import c_char
      character(kind=c_char), intent(in) :: subname(*)  ! SEE NOTE 1
    end subroutine

    !void gmvwrite_flag_data(int data_type, int flag_data[])
    subroutine gmvwrite_flag_data (data_type, flag_data) bind(c)
      import c_int
      integer(c_int), intent(in), value :: data_type
      integer(c_int), intent(in) :: flag_data(*)
    end subroutine

    !void gmvwrite_flag_endflag(void)
    subroutine gmvwrite_flag_endflag () bind(c)
    end subroutine

    !void gmvwrite_probtime(double ptime)
    subroutine gmvwrite_probtime (ptime) bind(c)
      import c_double
      real(c_double), intent(in), value :: ptime
    end subroutine

    !void gmvwrite_cycleno(int cyclenum)
    subroutine gmvwrite_cycleno (cyclenum) bind(c)
      import c_int
      integer(c_int), intent(in), value :: cyclenum
    end subroutine

    !void gmvwrite_nodeids(void *nodeids)
    subroutine gmvwrite_nodeids (nodeids) bind(c)
      import c_int
      integer(c_int), intent(in) :: nodeids(*)  ! SEE NOTE 2
    end subroutine

    !void gmvwrite_cellids(void *cellids)
    subroutine gmvwrite_cellids (cellids) bind(c)
      import c_int
      integer(c_int), intent(in) :: cellids(*)  ! SEE NOTE 2
    end subroutine

    !void gmvwrite_surface_header(void *nsurf)
    subroutine gmvwrite_surface_header (nsurf) bind(c)
      import c_int
      integer(c_int), intent(in) :: nsurf  ! SEE NOTE 2
    end subroutine

    !void gmvwrite_surface_data(int nverts, void *nodeids)
    subroutine gmvwrite_surface_data (nverts, nodeids) bind(c)
      import c_int
      integer(c_int), intent(in), value :: nverts
      integer(c_int), intent(in) :: nodeids(*)  ! SEE NOTE 2
    end subroutine

    !void gmvwrite_surfvars_header()
    subroutine gmvwrite_surfvars_header () bind(c)
    end subroutine

    !void gmvwrite_surfvars_name_data(char varname[], void *vids)
    subroutine gmvwrite_surfvars_name_data (varname, vids) bind(c)
      import c_char, c_double
      character(kind=c_char), intent(in) :: varname(*)  ! SEE NOTE 1
      real(c_double), intent(in) :: vids(*) ! SEE NOTE 2
    end subroutine

    !void gmvwrite_surfvars_endsvar()
    subroutine gmvwrite_surfvars_endsvar () bind(c)
    end subroutine

  end interface

end module gmvwrite_fortran_binding
