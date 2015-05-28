# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /sandbox/group2/apps/CMakeLists.txt
#
# CMakeLists.txt file for the sandbox/group2/apps apps.
# ===========================================================================

cmake_minimum_required (VERSION 2.8.2)
project (seqan_sandbox_group2_apps)

# Register all apps in core.  The apps will be built in SEQAN_RELEASE mode or
# the apps called ${app} in APP:${app} mode.
seqan_register_apps ()
