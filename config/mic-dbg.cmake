set(CMAKE_C_COMPILER icc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER ifort CACHE STRING "Fortran Compiler")

# The C/C++ compiler flags default to acceptable CMake-defined values based
# on the compiler ID and build type.  Not so for Fortran Flags -- we must
# specify them (or get none).
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")
set(CMAKE_C_FLAGS "-mmic" CACHE STRING "C compile flags")
set(CMAKE_Fortran_FLAGS "-g -u -mmic -assume realloc_lhs" CACHE STRING "Fortran compile flags")
