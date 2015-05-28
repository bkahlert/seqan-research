# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox\\mysandbox\\apps/CMakeLists.txt
#
# CMakeLists.txt file for the sandbox\\mysandbox\\apps apps.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_mysandbox_apps)

# Register all apps in core.  The apps will be built in SEQAN_RELEASE mode or
# the apps called ${app} in APP:${app} mode.
seqan_register_apps ()
