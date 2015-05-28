# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/PMSB_group6/tests/GSearch/CMakeLists.txt
#
# CMakeLists.txt file for the GSearch module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_PMSB_group6_tests_GSearch)
message (STATUS "Configuring sandbox/PMSB_group6/tests/GSearch")

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
add_executable (test_GSearch
                test_GSearch.cpp
                test_GSearch_systemTest.h
                test_GCheckClusterMaster.h
                test_GComputeMatch.h
                test_GParseMasterId.h)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_GSearch ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_GSearch COMMAND $<TARGET_FILE:test_GSearch>)
