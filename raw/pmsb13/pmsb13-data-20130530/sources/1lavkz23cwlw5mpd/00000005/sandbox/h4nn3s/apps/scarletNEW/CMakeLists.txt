cmake_minimum_required (VERSION 2.6)
project (sandbox_h4nn3s_apps_scarlet)
message (STATUS "Configuring sandbox/h4nn3s/apps/scarlet")

# Update the list of file names below if you add source files to your application.
seqan_add_executable(scarlet_indexer scarlet_indexer.cpp)

seqan_add_executable(scarlet scarlet.cpp)
