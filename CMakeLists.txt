cmake_minimum_required(VERSION 2.8.11)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
                      "${CMAKE_SOURCE_DIR}/cmake/Modules"
                      "${CMAKE_SOURCE_DIR}/cmake/Utils")

project(Pececillo Fortran C)

# Override the default installation prefix (/usr/local).
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${Pececillo_SOURCE_DIR}/install
      CACHE PATH "Installation prefix" FORCE)
endif()

# External project (netCDF, HDF5, HYPRE,etc.) install prefix
if (NOT EXTERNAL_INSTALL_PREFIX)
  set(EXTERNAL_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/external)
endif()

# ---------------------------------------------------------------------------- #
# Build Options
# ---------------------------------------------------------------------------- #

option(ENABLE_TESTS "Build test programs" ON)

# Build the tests
if (ENABLE_TESTS)
  enable_testing()
endif()

#-------------------------------------------------------------------------------
# Compiler Definitions
#-------------------------------------------------------------------------------

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

message(STATUS "CMAKE_Fortran_COMPILER_ID=${CMAKE_Fortran_COMPILER_ID}")
if (CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
 set (CMAKE_Fortran_FLAGS_RELEASE "$-u -O3 -w=uda -DNDEBUG")
 set (CMAKE_Fortran_FLAGS_DEBUG   "-u -g90")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -assume realloc_lhs")
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback -assume realloc_lhs -fp-stack-check -ftrapuv -check bounds,uninit,noarg-temp-created")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback -assume realloc_lhs -fp-stack-check -check bounds,noarg-temp-created")
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-O3 -DNDEBUG -assume realloc_lhs")
 message(STATUS "${CMAKE_Fortran_COMPILER_VERSION}")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_C_FLAGS}")
endif (OPENMP_FOUND)

#-------------------------------------------------------------------------------
# TPLs
#-------------------------------------------------------------------------------

# Find and/or build the required TPLs.
add_subdirectory(external)

# We need to incorporate netcdf and exodus into the TPL build process, but for
# now we require them to be prebuilt and installed someplace.  EXODUS_LIB_DIR
# must be set on the command line to the directory containing the libraries.
# Source tarballs for the packages are in external/tarballs.
if(NOT EXODUS_LIB_DIR)
  message(FATAL_ERROR "ERROR! EXODUS_LIB_DIR must be set to the directory "
          "containing the netcdf and exoIIv2c libraries.")
endif()
find_library(EXODUS_LIBRARY exoIIv2c ${EXODUS_LIB_DIR} NO_DEFAULT_PATH)
find_library(NETCDF_LIBRARY netcdf   ${EXODUS_LIB_DIR} NO_DEFAULT_PATH)
if(EXODUS_LIBRARY STREQUAL "EXODUS_LIBRARY-NOTFOUND" OR
   NETCDF_LIBRARY STREQUAL "NETCDF_LIBRARY-NOTFOUND")
  message(FATAL_ERROR "ERROR! Could not find the netcdf and/or exoIIv2c "
          "libraries in EXODUS_LIB_DIR=${EXODUS_LIB_DIR}")
endif()
add_library(exodus STATIC IMPORTED)
set_target_properties(exodus PROPERTIES
    IMPORTED_LOCATION ${EXODUS_LIBRARY}
    IMPORTED_LINK_INTERFACE_LIBRARIES ${NETCDF_LIBRARY})
message(STATUS "NETCDF_LIBRARY=${NETCDF_LIBRARY}")

#-------------------------------------------------------------------------------
# Installation Definitions
#-------------------------------------------------------------------------------

add_subdirectory(petaca)
set(Petaca_MODULE_DIR ${Petaca_BINARY_DIR}/mod_files/)

set(Pececillo_MODULE_DIR ${Pececillo_BINARY_DIR}/mod_files/)
set(CMAKE_Fortran_MODULE_DIRECTORY ${Pececillo_MODULE_DIR})
add_subdirectory(src)

