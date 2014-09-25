set(CMAKE_C_COMPILER gcc CACHE STRING "C Compiler")
set(CMAKE_Fortran_COMPILER nagfor CACHE STRING "Fortran Compiler")

# The C/C++ compiler flags default to acceptable CMake-defined values that
# depend on the compiler ID and build type.  Not so for Fortran Flags -- we
# must specify them (or get none).
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")
set(CMAKE_Fortran_FLAGS "-u -C -C=dangling -g90 -gline -nan" CACHE STRING "Fortran compile flags")
