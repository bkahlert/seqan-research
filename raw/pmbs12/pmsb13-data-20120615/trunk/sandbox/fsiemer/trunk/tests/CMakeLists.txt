cmake_minimum_required (VERSION 2.6)
project (seqan_sandbox_fsiemer_tests)

seqan_setup_tests(sandbox_fsiemer_tests)  # Setup flags etc, add meta target.

seqan_add_all_subdirectories()  # include_directory(*)
