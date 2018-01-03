cmake_minimum_required(VERSION 3.5)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
                      "${CMAKE_SOURCE_DIR}/cmake/Modules"
                      "${CMAKE_SOURCE_DIR}/cmake/Utils")

include(NoInSourceBuild)

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
  # note: need to disable openmp below for traceback support
  set (CMAKE_Fortran_FLAGS_RELEASE "-u -O3 -w=uda -DNDEBUG")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g90 -g -gline -O0 -w=uda -C=all")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -g -DNDEBUG -assume realloc_lhs -traceback")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-g -assume realloc_lhs -traceback -fp-stack-check -ftrapuv -check bounds,uninit,noarg-temp-created,pointer")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_C_FLAGS}")
  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_C_FLAGS}")
endif (OPENMP_FOUND)

#-------------------------------------------------------------------------------
# External dependencies
#-------------------------------------------------------------------------------

set(CMAKE_PREFIX_PATH ${PECECILLO_TPL_DIR})

# HYPRE library
find_package(HYPRE 2.9.0)
if(HYPRE_FOUND)
  if(HYPRE_IS_PARALLEL)
    set(HYPRE_FOUND False)
    message(STATUS "Require serial HYPRE library but found unsuitable parallel library")
  endif()
endif()
if(NOT HYPRE_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Hypre installation")
endif()

# HDF5 library
find_package(HDF5 "1.8.8" COMPONENTS C HL)
if(HDF5_FOUND)
  if(HDF5_IS_PARALLEL)
    set(HDF5_FOUND False)
    message(STATUS "Require serial HDF5 library but found unsuitable parallel library")
  endif()
endif()
if(HDF5_FOUND)
  add_library(hdf5 INTERFACE)
  target_link_libraries(hdf5 INTERFACE ${HDF5_LIBRARIES})
  target_include_directories(hdf5 INTERFACE ${HDF5_INCLUDE_DIRS})
  target_compile_definitions(hdf5 INTERFACE ${HDF5_DEFINITIONS})
else()
  message(FATAL_ERROR "Could NOT find a suitable HDF5 installation")
endif()

# NetCDF library
find_package(NetCDF "4.1.3")
if(NETCDF_FOUND)
  if(NOT NETCDF_HAS_NC4)
    message(STATUS "Found unsuitable NetCDF without the required netcdf-4 feature")
    set(NETCDF_FOUND False)
  endif()
endif()
if(NOT NETCDF_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable NetCDF installation")
endif()

# ExodusII library
find_package(Exodus "602" EXACT)
if(NOT EXODUS_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable Exodus installation")
endif()

# YAJL library
find_package(YAJL)
if(NOT YAJL_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable YAJL library")
endif()

# LAPACK library
find_package(LAPACK)
if(NOT LAPACK_FOUND)
  message(FATAL_ERROR "Could NOT find a suitable LAPACK library")
endif()

#-------------------------------------------------------------------------------
# Installation Definitions
#-------------------------------------------------------------------------------

# RPATH handling for installed shared libraries and dynamically-linked
# executables.  See http://www.cmake.org/Wiki/CMake_RPATH_handling
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_subdirectory(src)
