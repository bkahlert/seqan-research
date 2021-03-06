# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/alphons/SeqAn/seqan-trunk

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/alphons/SeqAn/seqan-trunk

# Include any dependencies generated for this target.
include sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/depend.make

# Include the progress variables for this target.
include sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/progress.make

# Include the compile flags for this target's objects.
include sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/flags.make

sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o: sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/flags.make
sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o: sandbox/my_sandbox/apps/lagan/lagan.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/alphons/SeqAn/seqan-trunk/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o"
	cd /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/lagan.dir/lagan.cpp.o -c /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan/lagan.cpp

sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lagan.dir/lagan.cpp.i"
	cd /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan/lagan.cpp > CMakeFiles/lagan.dir/lagan.cpp.i

sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lagan.dir/lagan.cpp.s"
	cd /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan/lagan.cpp -o CMakeFiles/lagan.dir/lagan.cpp.s

sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o.requires:
.PHONY : sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o.requires

sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o.provides: sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o.requires
	$(MAKE) -f sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/build.make sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o.provides.build
.PHONY : sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o.provides

sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o.provides.build: sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o

# Object files for target lagan
lagan_OBJECTS = \
"CMakeFiles/lagan.dir/lagan.cpp.o"

# External object files for target lagan
lagan_EXTERNAL_OBJECTS =

bin/lagan: sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o
bin/lagan: /usr/lib/i386-linux-gnu/libz.so
bin/lagan: sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/build.make
bin/lagan: sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable ../../../../bin/lagan"
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Pre Build Instrumentation..."
	cd /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan && /usr/bin/python /home/alphons/SeqAn/seqan-trunk/misc/seqan_instrumentation/bin/seqan_instrumentation.py pre_build /home/alphons/SeqAn/seqan-trunk /home/alphons/SeqAn/seqan-trunk lagan
	cd /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lagan.dir/link.txt --verbose=$(VERBOSE)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Post Build Instrumentation..."
	cd /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan && /usr/bin/python /home/alphons/SeqAn/seqan-trunk/misc/seqan_instrumentation/bin/seqan_instrumentation.py post_build /home/alphons/SeqAn/seqan-trunk /home/alphons/SeqAn/seqan-trunk lagan

# Rule to build all files generated by this target.
sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/build: bin/lagan
.PHONY : sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/build

sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/requires: sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/lagan.cpp.o.requires
.PHONY : sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/requires

sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/clean:
	cd /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan && $(CMAKE_COMMAND) -P CMakeFiles/lagan.dir/cmake_clean.cmake
.PHONY : sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/clean

sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/depend:
	cd /home/alphons/SeqAn/seqan-trunk && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/alphons/SeqAn/seqan-trunk /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan /home/alphons/SeqAn/seqan-trunk /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan /home/alphons/SeqAn/seqan-trunk/sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : sandbox/my_sandbox/apps/lagan/CMakeFiles/lagan.dir/depend

