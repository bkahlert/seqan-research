cmake_minimum_required (VERSION 2.6)
project (sandbox_meyerclp_apps_dreme)

# Update the list of file names below if you add source files to your application.
seqan_add_executable(dreme dreme.cpp dreme.h)

# Tests for DREME code.
seqan_add_test_executable(dreme_test dreme_test.cpp)
