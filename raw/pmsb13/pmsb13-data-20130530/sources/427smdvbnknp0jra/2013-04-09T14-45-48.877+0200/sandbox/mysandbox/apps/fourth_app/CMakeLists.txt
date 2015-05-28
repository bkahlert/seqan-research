# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox\\mysandbox\\apps\\fourth_app/CMakeLists.txt
#
# CMakeLists.txt file for fourth_app.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_mysandbox_apps_fourth_app)
message (STATUS "Configuring sandbox\\mysandbox\\apps\\fourth_app")

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
add_executable (fourth_app fourth_app.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (fourth_app ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (fourth_app)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install fourth_app in ${PREFIX}/bin directory
install (TARGETS fourth_app
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/fourth_app for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (fourth_app)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:fourth_app")
  set (CPACK_PACKAGE_NAME "fourth_app")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "fourth_app")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Your Name <your.email@example.net>")
  set (CPACK_PACKAGE_VENDOR "Your Name <your.email@example.net>")

  seqan_configure_cpack_app (fourth_app "fourth_app")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:fourth_app")

