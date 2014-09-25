project(PececilloTPL)

message(STATUS "Configuring the external packages")

# Add
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

# Add ExternalProject_Add to scope
include(ExternalProject)

# If not specified, the prefix for the external package installation path
# defaults to the "external" subdirectory of CMAKE_INSTALL_PREFIX.  This
# is also one location cmake will look for pre-installed packages.
if(NOT TPL_INSTALL_PREFIX)
  set(TPL_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/external)
endif()

# --- Useful directory variables (GET RID OF THESE)
set(TPL_TEMPLATES_DIR ${CMAKE_CURRENT_LIST_DIR}/templates)
set(TPL_ARCHIVE_DIR   ${CMAKE_CURRENT_LIST_DIR}/tarfiles)
set(TPL_PATCHES_DIR   ${CMAKE_CURRENT_LIST_DIR}/patches)

# GNU AUTOCONF configure static/shared
set(TPL_SHARED_SWITCH --disable-shared)
if(BUILD_SHARED_LIBS)
  set(TPL_SHARED_SWITCH --enable-shared)
endif()

# Log all activity -- log files for each step in the timestamp directory
set(TPL_LOG_OPTS
    LOG_DOWNLOAD  1
    LOG_UPDATE    1
    LOG_CONFIGURE 1
    LOG_BUILD     1
    LOG_TEST      1
    LOG_INSTALL   1
)

#include(Verify_NetCDF)
#if(NOT NetCDF_VERIFIED)
#  include(ExternalProject_NetCDF)
#endif()

include(Verify_HYPRE)
if(NOT HYPRE_VERIFIED)
  include(ExternalProject_HYPRE)
endif()

# Should be some way of doing an imported library thing, no?

#set(HYPRE_FOUND        ${HYPRE_FOUND} PARENT_SCOPE)
#set(HYPRE_BUILD_TARGET ${HYPRE_BUILD_TARGET} PARENT_SCOPE)
set(HYPRE_VERSION      ${HYPRE_VERSION} PARENT_SCOPE)
set(HYPRE_IS_PARALLEL  ${HYPRE_IS_PARALLEL} PARENT_SCOPE)
#set(HYPRE_LIBRARY      ${HYPRE_LIBRARY} PARENT_SCOPE)
set(HYPRE_LIBRARIES    ${HYPRE_LIBRARIES} PARENT_SCOPE)
#set(HYPRE_INCLUDE_DIR  ${HYPRE_INCLUDE_DIR} PARENT_SCOPE)
set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIRS} PARENT_SCOPE)
