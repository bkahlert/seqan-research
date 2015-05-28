# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox\\my_sandbox\\apps\\query_fai/CMakeLists.txt
#
# CMakeLists.txt file for query_fai.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_my_sandbox_apps_query_fai)
message (STATUS "Configuring sandbox\\my_sandbox\\apps\\query_fai")

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
add_executable (query_fai query_fai.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (query_fai ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (query_fai)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install query_fai in ${PREFIX}/bin directory
install (TARGETS query_fai
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/query_fai for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (query_fai)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:query_fai")
  set (CPACK_PACKAGE_NAME "query_fai")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "query_fai")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (query_fai "query_fai")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:query_fai")

