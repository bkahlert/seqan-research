# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/h4nn3s/apps/scarlet/CMakeLists.txt
#
# CMakeLists.txt file for scarlet.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_h4nn3s_apps_scarlet)
message (STATUS "Configuring sandbox/h4nn3s/apps/scarlet")

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
add_executable (scarlet scarlet.cpp
			scarlet.hpp
			misc.hpp
			finder.hpp
			options.hpp)
add_executable (scarlet_indexer scarlet_indexer.cpp
			options.hpp
			misc.hpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (scarlet ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (scarlet)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install scarlet in ${PREFIX}/bin directory
install (TARGETS scarlet scarlet_indexer
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/scarlet for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (scarlet)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:scarlet")
  set (CPACK_PACKAGE_NAME "scarlet")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "scarlet")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (scarlet "scarlet")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:scarlet")

