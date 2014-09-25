Here is a brief description of each of the source files contained in this
directory.  The hope is that this provides a good overview of how the code
and its functionality are organized.  These comprise the core pececillo
library upon which individual mini-app PDE solvers will be based.

Neil Carlson <nnc@lanl.gov>
September 2014

================================================================================

GENERAL UTILITIES MODULES

The following modules provide general utilities and are largely standalone
and independent.

f90_assert.F90/.fpp
  Implements assertions via preprocessor macros.

kinds.F90
  Provides kind parameter values for 4 and 8-byte integers and reals.  With
  the addition of these types of parameters to the intrinsic iso_fortran_env
  module in F2008, this is perhaps no longer needed; however these parameters
  have the virtue of much shorter names (i4, i8, r4, r8) than the new ones.
  
string_utilities.F90
  Provides a few handy string functions, like converting an integer to a
  string and raising/lowering the case of string.

bitfield_type.F90
  To reduce memory requirements, it is often useful to use an integer as a
  bitmask to record collections of boolean data.  When the 32 or 64 bits of
  data provided by an integer is not enough, the bitfield type is what you
  need.

permutations.F90
  Provides procedures for manipulating permutations and doing in-place
  reordering of arrays according to a provided permutation.

graph_type.F90
  Defines a derived type for representing an N-node graph.  The implemented
  functionality is limited to dynamically building a graph and then extracting
  a compact static description of the graph in a standard format.

logging_services.F90
  Provides a common set of procedures for writing log/trace messages -- the
  types of messages typically written to the terminal or simulation log file.
  Its use it intended to produce more consistent looking output and eliminate
  raw Fortran formatted output code.
  
data_layout_type.F90
  Defines a data type for managing the layout of subcomponents within a
  single index space.  Useful for multi-physics problems, for example,
  where multiple variables must be packed into a single vector of unknowns.
 
 
MESH MODULES
  
unstr_mesh_type.F90
  Defines a derived type that describes the topology and geometry of an
  unstructured FE-like mesh.  The data components are public and meant to
  be accessed directly by application code.

unstr_mesh_factory.F90
  Provides several routines for instantiating a mesh object.  The mesh
  can either be read from an Exodus/Genesis format file, or internally
  generated (a simple regular cartesian mesh useful for testing). 

exodus_c_binding.F90
  Defines bindings to the Exodus C library which is the official interface
  to the Exodus file.  The underlying storage format is NetCDF.

exodus_file_type.F90
  An application-specific layer over the Exodus library.  This is what
  the mesh code uses.
    
cell_topology.F90
  Provides data and procedures that establish the topology of different cell
  types; e.g., how the the vertices of a hexahedral cell are ordered.
  
cell_geometry.F90
  Provides procedures for computing primitive cell geometry data, like
  volumes, face areas, face normals, etc.

unstr_mesh_tools.F90
  Provides procedures for generating higher-order mesh topology information,
  such as a numbering of the mesh faces.  Much of the work in instantiating
  a mesh object is done here.

facet_hash_type.F90
  Defines a hash function specialized for hashing facets of mesh cells.

facet_table_type.F90
  Defines a facet hash table (using the facet hash function).  This is used
  to enumerate mesh faces and establish their orientations, for example.


FUNCTION MODULES

The need for general functions occurs in many places in a PDE code; the
specification of boundary data, PDE coefficients, and initial conditions are
just a few examples.  To avoid hardcoding specific types of functions at each
place where one is needed, it is necessary to abstract out the notion of a
function as an object.  These modules provide this capability.

scalar_func_class.F90
  Defines an abstract base class that provides an interface to a general
  scalar-valued function of a vector argument.  Application code only deals
  with polymorphic objects of this type.

scalar_func_factories.F90
  Provides procedures for instantiating a scalar_func class object of
  different dynamic types.

scalar_func_containers.F90
  Application code sometimes needs to deal with a collection (e.g., array)
  of scalar_func objects, each with a possibly different dynamic type.
  Unfortunately a simple Fortran array cannot handle this because arrays
  must have a homogenous dynamic type (every element is the same type).
  This module provides some containers that work around this limitation.

The following modules provide concrete implementations of the scalar_func
class.  Application code should never be using these directly.
  
const_scalar_func_type.F90
  Defines a constant-valued function.

poly_scalar_func_type.F90
  Defines a polynomial function.

mpoly_scalar_func_type.F90
  Defines a multivariable polynomial function.
  
tabular_scalar_func_type.F90
  Defines a tabular function.

fptr_scalar_func_type.F90
  Wraps an application-provided parametrized function.


MESH FUNCTION MODULES

As for general continuous functions above, there is also a need for discrete
functions over a mesh; for example a scalar function of time and mesh cells
or time and boundary faces.  Usually these arise as the discretization of a
general continuous function over the mesh or mesh boundary.  Here again, we
need to abstract out the essential characteristics as an object in order to
avoid hard-coding specific cases in every place where they are needed.

bndry_func_class.F90
  Defines an abstract base class that provides an interface to a time-
  dependent function of a subset of mesh facets (e.g., cells or faces).
  Application code deals with polymorphic objects of this type and not
  with any specific concrete implementation.

bc_factory_type.F90
  Defines an object for instantiating bndry_func objects of different
  dynamic types.

bndry_face_func_type.F90
  A concrete implementation of bndry_func that is based on an unstr_mesh
  object and discretizes scalar_func objects over portions of the mesh
  boundary.

unstr_mesh_func.F90
  Sets the values of a cell-based array over an unstr_mesh object using
  user-specified scalar_func objects.  Think of this as a time-independent
  mesh function.  This does not derive from bndry_func.


LINEAR ALGEBRA MODULES

csr_matrix_type.F90
  Provides a type for representing a sparse matrix stored in compressed
  sparse row format (CSR).

hypre_c_binding.F90
  Defines bindings to a subset of the HYPRE C library.

hypre_ext.c
  Some higher-level "extensions" to the HYPRE C library interface.  These
  bundle together some bits of boiler plate calls, which were more easily
  done in C than in Fortran.

fhypre.F90
  A Fortran-centric interface to HYPRE.  This uses hypre_c_binding and
  hypre_ext, but tweaks the interface to feel more comfortable to Fortran.

hypre_pcg_type.F90
  Provides a preconditioned CG solver class for csr_matrix objects and is built
  on the ParCSR PCG solver from Hypre that uses BoomerAMG as the preconditioner.

hypre_fgmres_type.F90
  Provides an FGMRES solver class for csr_matrix objects and is built on the
  ParCSR FlexGMRES solver from Hypre that uses BoomerAMG as the preconditioner.

hypre_hybrid_type.F90
  Provides a preconditioned CG solver class for csr_matrix objects and is built
  on the ParCSR Hybrid solver from Hypre that initially uses diagonal scaling
  as the preconditioner, switching to BoomerAMG only when the convergence
  rate is too slow.

csr_precon_class.F90
  Defines an abstract base class that provides an interface to a generic
  preconditioner for a csr_matrix object.  Application code will form
  equation-based preconditioners that make use of preconditioners for CSR
  matrix pieces of the full matrix.  This abstraction makes the application
  code independent of the specific preconditioner used.

csr_precon_ssor_type.F90
  A concrete implementation of csr_precon that uses SSOR sweeps.

csr_precon_boomer_type.F90
  A concrete implementation of csr_precon that uses BoomerAMG from Hypre.

upper_packed_matrix.F90
  Provides Lapack-like routines for symmetric positive definite matrices
  stored in an upper packed format.  This is intended for the small matrices
  (e.g., 6x6) arising from local cell-based equations.  No special data
  structure is defined to hold the matrix -- the upper triangle is stored in
  a normal rank-1 array.


MIMETIC DISCRETIZATION MODULES

mfd_disc_type.F90
  Defines an object that provides computational procedures for some basic
  mimetic finite difference operations.  This wraps an unstr_mesh object.

mfd_diff_matrix_type.F90
  Defines an object that encapsulates the MFD diffusion matrix over an
  unstr_mesh object.

mfd_diff_precon_type.F90
  Defines a preconditioner object for a MFD diffusion matrix.  This is built
  in part on a csr_precon class object.


ODE INTEGRATOR MODULES

idaesol_type.F90
  Provides the BDF2 integrator object for implicit systems of DAE (index-1).

nka_type.F90
  Provides the nonlinear Krylov accelerator object that the BDF2 integrator
  uses in its nonlinear solver.


OUTPUT MODULES

gmvwrite.c
  GMV-provided C library for writing GMV-format visualization files
  
gmvwrite_fortran_binding.F90
  Fortran bindings to gmvwrite.c.

unstr_mesh_gmv.F90
  Application-specific layer over gmvwrite_fortran_bindings.  Provides
  procedures for writing an unstr_mesh object and cell and node based
  fields over that mesh.
