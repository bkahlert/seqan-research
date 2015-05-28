# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/my_sandbox/CMakeLists.txt
#
# CMakeLists.txt file for sandbox/my_sandbox.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_my_sandbox)

# Add the paths sandbox/my_sandbox/include and sandbox/my_sandbox/include to the paths that CMake
# searches for libraries.
set (CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/../../extras/include
                                              ${CMAKE_CURRENT_SOURCE_DIR}/../../core/include)

# Setup the library modules from core and extras.
seqan_setup_library (sandbox_my_sandbox core extras)

# Add subdirectory for apps.
add_subdirectory (apps)

add_subdirectory (flexbar-trunk)

# Demos are required when doing a Whole SeqAn Release (copy demos) or when
# developing (build demos).
if (("${SEQAN_BUILD_SYSTEM}" STREQUAL "SEQAN_RELEASE") OR
    ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP"))
    add_subdirectory (demos)
endif ()

# Tests are only built when building in DEVLOP mode.
if ("${SEQAN_BUILD_SYSTEM}" STREQUAL "DEVELOP")
    add_subdirectory (tests)
endif ()
