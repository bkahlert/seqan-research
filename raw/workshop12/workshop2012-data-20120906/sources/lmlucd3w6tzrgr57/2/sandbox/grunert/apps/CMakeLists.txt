cmake_minimum_required (VERSION 2.6)
project (seqan_sandbox_grunert_apps)

seqan_setup_apps(sandbox_grunert_apps)  # Add meta target and make executables depend on it.

seqan_add_all_subdirectories()  # include_directory(*)