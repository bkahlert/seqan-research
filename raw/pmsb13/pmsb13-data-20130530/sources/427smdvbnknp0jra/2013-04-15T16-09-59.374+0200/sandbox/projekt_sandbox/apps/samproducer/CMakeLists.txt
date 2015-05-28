# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox\\projekt_sandbox\\apps\\samproducer/CMakeLists.txt
#
# CMakeLists.txt file for samproducer.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_projekt_sandbox_apps_samproducer)
message (STATUS "Configuring sandbox\\projekt_sandbox\\apps\\samproducer")

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
add_executable (samproducer samproducer.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (samproducer ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (samproducer)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install samproducer in ${PREFIX}/bin directory
install (TARGETS samproducer
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/samproducer for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (samproducer)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:samproducer")
  set (CPACK_PACKAGE_NAME "samproducer")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "samproducer")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (samproducer "samproducer")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:samproducer")

