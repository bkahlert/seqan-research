# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox\\mysandbox\\demos/CMakeLists.txt
#
# CMakeLists.txt file for the sandbox\\mysandbox\\demos demos.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_mysandbox_demos)
message (STATUS "Configuring sandbox\\mysandbox\\demos")

# Recursively register all demos in core and prefix their target/binary names
# with "demo_".  The demos will be built in DEVELOP mode, copied in
# SEQAN_RELEASE mode and left alone otherwise.
seqan_register_demos (demo_)
