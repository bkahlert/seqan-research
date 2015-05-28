# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox\\gruppe1\\tests\\PEMer_Lite_header/CMakeLists.txt
#
# CMakeLists.txt file for the PEMer_Lite_header module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_gruppe1_tests_PEMer_Lite_header)
message (STATUS "Configuring sandbox\\gruppe1\\tests\\PEMer_Lite_header")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES NONE)
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Update the list of file names below if you add source files to your test.
add_executable (test_PEMer_Lite_header
                test_PEMer_Lite_header.cpp
                test_PEMer_Lite_header.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_PEMer_Lite_header ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_PEMer_Lite_header COMMAND $<TARGET_FILE:test_PEMer_Lite_header>)
