# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/pmsb_group4/apps/lagan_svn/CMakeLists.txt
#
# CMakeLists.txt file for lagan_svn.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_pmsb_group4_apps_lagan_svn)
message (STATUS "Configuring sandbox/pmsb_group4/apps/lagan_svn")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
#set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your application.
add_executable (lagan_svn lagan_svn.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (lagan_svn ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (lagan_svn)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install lagan_svn in ${PREFIX}/bin directory
install (TARGETS lagan_svn
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/lagan_svn for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (lagan_svn)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:lagan_svn")
  set (CPACK_PACKAGE_NAME "lagan_svn")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "lagan_svn")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (lagan_svn "lagan_svn")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:lagan_svn")

