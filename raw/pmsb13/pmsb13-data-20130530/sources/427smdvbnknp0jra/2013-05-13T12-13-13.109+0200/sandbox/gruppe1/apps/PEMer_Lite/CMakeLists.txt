# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox\\gruppe1\\apps\\PEMer_Lite/CMakeLists.txt
#
# CMakeLists.txt file for PEMer_Lite.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_gruppe1_apps_PEMer_Lite)
message (STATUS "Configuring sandbox\\gruppe1\\apps\\PEMer_Lite")

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
add_executable (PEMer_Lite PEMer_Lite.cpp)
add_executable (PEMer_Lite PEMer_Lite.h)
add_executable (PEMer_Lite_test PEMer_Lite_test.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (PEMer_Lite ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (PEMer_Lite)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install PEMer_Lite in ${PREFIX}/bin directory
install (TARGETS PEMer_Lite
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/PEMer_Lite for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (PEMer_Lite)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:PEMer_Lite")
  set (CPACK_PACKAGE_NAME "PEMer_Lite")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "PEMer_Lite")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (PEMer_Lite "PEMer_Lite")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:PEMer_Lite")

