cmake_minimum_required (VERSION 2.6)
project (sandbox_meyerclp_tests_ChIPSeq)

# Update the list of file names below if you add source files to your test.
seqan_add_test_executable(test_ChIPSeq test_ChIPSeq.cpp test_ChIPSeq.h)
