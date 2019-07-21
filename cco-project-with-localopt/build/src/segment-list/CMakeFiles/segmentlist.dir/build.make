# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.11

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build

# Include any dependencies generated for this target.
include src/segment-list/CMakeFiles/segmentlist.dir/depend.make

# Include the progress variables for this target.
include src/segment-list/CMakeFiles/segmentlist.dir/progress.make

# Include the compile flags for this target's objects.
include src/segment-list/CMakeFiles/segmentlist.dir/flags.make

src/segment-list/CMakeFiles/segmentlist.dir/segment-list.cpp.o: src/segment-list/CMakeFiles/segmentlist.dir/flags.make
src/segment-list/CMakeFiles/segmentlist.dir/segment-list.cpp.o: ../src/segment-list/segment-list.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/segment-list/CMakeFiles/segmentlist.dir/segment-list.cpp.o"
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/segment-list && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/segmentlist.dir/segment-list.cpp.o -c /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/src/segment-list/segment-list.cpp

src/segment-list/CMakeFiles/segmentlist.dir/segment-list.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/segmentlist.dir/segment-list.cpp.i"
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/segment-list && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/src/segment-list/segment-list.cpp > CMakeFiles/segmentlist.dir/segment-list.cpp.i

src/segment-list/CMakeFiles/segmentlist.dir/segment-list.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/segmentlist.dir/segment-list.cpp.s"
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/segment-list && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/src/segment-list/segment-list.cpp -o CMakeFiles/segmentlist.dir/segment-list.cpp.s

# Object files for target segmentlist
segmentlist_OBJECTS = \
"CMakeFiles/segmentlist.dir/segment-list.cpp.o"

# External object files for target segmentlist
segmentlist_EXTERNAL_OBJECTS =

src/segment-list/libsegmentlist.a: src/segment-list/CMakeFiles/segmentlist.dir/segment-list.cpp.o
src/segment-list/libsegmentlist.a: src/segment-list/CMakeFiles/segmentlist.dir/build.make
src/segment-list/libsegmentlist.a: src/segment-list/CMakeFiles/segmentlist.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libsegmentlist.a"
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/segment-list && $(CMAKE_COMMAND) -P CMakeFiles/segmentlist.dir/cmake_clean_target.cmake
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/segment-list && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/segmentlist.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/segment-list/CMakeFiles/segmentlist.dir/build: src/segment-list/libsegmentlist.a

.PHONY : src/segment-list/CMakeFiles/segmentlist.dir/build

src/segment-list/CMakeFiles/segmentlist.dir/clean:
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/segment-list && $(CMAKE_COMMAND) -P CMakeFiles/segmentlist.dir/cmake_clean.cmake
.PHONY : src/segment-list/CMakeFiles/segmentlist.dir/clean

src/segment-list/CMakeFiles/segmentlist.dir/depend:
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/src/segment-list /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/segment-list /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/segment-list/CMakeFiles/segmentlist.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/segment-list/CMakeFiles/segmentlist.dir/depend

