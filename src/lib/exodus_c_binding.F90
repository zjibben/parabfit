!!
!! EXODUS_C_BINDING
!!
!! Bindings to a subset of the ExodusII C library.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! January 2014
!!
!! NOTES
!!
!!  Bindings to just a very small of the Exodus C library are defined; only
!!  what is required to read the mesh data used by Truchas/Pececillo.
!!
!!  I've used the Exodus II document (2009 update) as the reference, but a
!!  quick scan of the library source (version 6) suggests that it is out of
!!  date (though not necessarily incorrect).  In particular it appears there
!!  is the possibility of using 64-bit integers for various entity types.
!!  Only a 32-bit integer (C int type) interface is documented, and that's
!!  what I've assumed here. This is a likely place of breakage for this module.
!!
!!  The documented ex_open function is actually a preprocessor macro that
!!  delegates to the auxillary function ex_open_int, passing an additional
!!  constant EX_API_VERS_NODOT defined by the exodusII.h header file.  For
!!  this Fortran binding, exo_open is defined as a module procedure.  The
!!  value of the parameter EX_API_VERS_NODOT must match the header file value
!!  used in the Exodus II library.
!!

module exodus_c_binding

  use,intrinsic :: iso_c_binding, only: c_int, c_float, c_char, c_ptr
  implicit none
  private

  !! Functions (see Exodus II manual)
  public :: ex_open                 ! open Exodus II file
  public :: ex_close                ! close Exodus II file
  public :: ex_get_init             ! read initialization parameters
  public :: ex_get_coord            ! read nodal coordinates
  public :: ex_get_elem_blk_ids     ! read element block IDs
  public :: ex_get_elem_block       ! read element block parameters
  public :: ex_get_elem_conn        ! read element block connectivity
  public :: ex_get_side_set_ids     ! read side set IDs
  public :: ex_get_side_set_param   ! read side set parameteters
  public :: ex_get_side_set         ! read side set

  !! Parameters from exodusII.h
  integer, parameter, public :: MAX_STR_LENGTH  = 32
  integer, parameter, public :: MAX_LINE_LENGTH = 80
  integer(c_int),  parameter, public :: EX_READ = 0
  integer(c_int),  parameter         :: EX_API_VERS_NODOT = 602 ! MUST MATCH THE INSTALLED LIBRARY

  interface
    function ex_open_int(path, mode, comp_ws, io_ws, version, run_version) &
        result(exoid) bind(c)
      import c_int, c_float, c_char
      character(kind=c_char), intent(in) :: path(*)
      integer(c_int), value, intent(in) :: mode
      integer(c_int), intent(inout) :: comp_ws, io_ws
      real(c_float), intent(out) :: version
      integer(c_int), value, intent(in) :: run_version
      integer(c_int) :: exoid
    end function
  end interface

  interface
    function ex_close(exoid) result(error) bind(c)
      import c_int
      integer(c_int), value, intent(in) :: exoid
      integer(c_int) :: error
    end function
  end interface

  interface
    function ex_get_init (exoid, title, num_dim, num_nodes, num_elem, &
        num_elem_blk, num_node_sets, num_side_sets) result(error) bind(c)
      import c_int, c_char
      integer(c_int), value, intent(in) :: exoid
      character(kind=c_char), intent(out) :: title(*) ! length <= MAX_LINE_LENGTH
      integer(c_int), intent(out) :: num_dim, num_nodes, num_elem, num_elem_blk
      integer(c_int), intent(out) :: num_node_sets, num_side_sets
      integer(c_int) :: error
    end function
  end interface

  interface
    function ex_get_elem_blk_ids(exoid, elem_blk_ids) result(error) bind(c)
      import c_int
      integer(c_int), value, intent(in) :: exoid
      integer(c_int), intent(out) :: elem_blk_ids(*)
      integer(c_int) :: error
    end function
  end interface

  interface
    function ex_get_elem_block(exoid, elem_blk_id, elem_type, &
        num_elem_this_blk, num_nodes_per_elem, num_attr) result(error) bind(c)
      import c_int, c_char
      integer(c_int), value, intent(in) :: exoid, elem_blk_id
      character(kind=c_char), intent(out) :: elem_type(*) ! length <= MAX_STR_LENGTH
      integer(c_int), intent(out) :: num_elem_this_blk, num_nodes_per_elem, num_attr
      integer(c_int) :: error
    end function
  end interface

  interface
    function ex_get_side_set_ids(exoid, side_set_ids) result(error) bind(c)
      import c_int
      integer(c_int), value, intent(in) :: exoid
      integer(c_int), intent(out) :: side_set_ids(*)
      integer(c_int) :: error
    end function
  end interface

  interface
    function ex_get_side_set_param(exoid, side_set_id, num_side_in_set, &
        num_dist_fact_in_set) result(error) bind(c)
      import c_int
      integer(c_int), value, intent(in) :: exoid, side_set_id
      integer(c_int), intent(out) :: num_side_in_set, num_dist_fact_in_set
      integer(c_int) :: error
    end function
  end interface

  interface
    function ex_get_side_set(exoid, side_set_id, side_set_elem_list, &
        side_set_side_list) result(error) bind(c)
      import c_int, c_ptr
      integer(c_int), value, intent(in) :: exoid, side_set_id
      type(c_ptr), value :: side_set_elem_list, side_set_side_list
      integer(c_int) :: error
    end function
  end interface

  interface
    function ex_get_elem_conn(exoid, elem_blk_id, connect) result(error) bind(c)
      import c_int, c_ptr
      integer(c_int), value, intent(in) :: exoid, elem_blk_id
      type(c_ptr), value :: connect
      integer(c_int) :: error
    end function
  end interface

  interface
    function ex_get_coord(exoid, x_coor, y_coor, z_coor) result(error) bind(c)
      import c_int, c_ptr
      integer(c_int), value, intent(in) :: exoid
      type(c_ptr), value :: x_coor, y_coor, z_coor
      integer(c_int) :: error
    end function
  end interface

contains

  !! EX_OPEN is actually defined as a macro in exodusII.h.
  function ex_open (path, mode, comp_ws, io_ws, version) result(exoid) bind(c)
    character(kind=c_char), intent(in) :: path(*)
    integer(c_int), value, intent(in) :: mode
    integer(c_int), intent(inout) :: comp_ws, io_ws
    real(c_float), intent(out) :: version
    integer(c_int) :: exoid
    exoid = ex_open_int(path, mode, comp_ws, io_ws, version, EX_API_VERS_NODOT)
  end function ex_open

end module exodus_c_binding
