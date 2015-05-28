# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/my_sandbox/tests/test_lagan/CMakeLists.txt
#
# CMakeLists.txt file for the test_lagan module tests.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_my_sandbox_tests_test_lagan)
message (STATUS "Configuring sandbox/my_sandbox/tests/test_lagan")

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
add_executable (test_test_lagan
                test_test_lagan.cpp
                test_test_lagan.h
                ../../apps/lagan_omp/lagan_functions.h
                ../../apps/lagan_omp/lagan_functions.cpp)

# Add dependencies found by find_package (SeqAn).
target_link_libraries (test_test_lagan ${SEQAN_LIBRARIES})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Register with CTest
# ----------------------------------------------------------------------------

add_test (NAME test_test_test_lagan COMMAND $<TARGET_FILE:test_test_lagan>)
