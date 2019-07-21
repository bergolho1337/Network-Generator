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
include src/local_optimization_library/CMakeFiles/default_local_optimization.dir/depend.make

# Include the progress variables for this target.
include src/local_optimization_library/CMakeFiles/default_local_optimization.dir/progress.make

# Include the compile flags for this target's objects.
include src/local_optimization_library/CMakeFiles/default_local_optimization.dir/flags.make

src/local_optimization_library/CMakeFiles/default_local_optimization.dir/local_optimization.cpp.o: src/local_optimization_library/CMakeFiles/default_local_optimization.dir/flags.make
src/local_optimization_library/CMakeFiles/default_local_optimization.dir/local_optimization.cpp.o: ../src/local_optimization_library/local_optimization.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/local_optimization_library/CMakeFiles/default_local_optimization.dir/local_optimization.cpp.o"
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/local_optimization_library && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/default_local_optimization.dir/local_optimization.cpp.o -c /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/src/local_optimization_library/local_optimization.cpp

src/local_optimization_library/CMakeFiles/default_local_optimization.dir/local_optimization.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/default_local_optimization.dir/local_optimization.cpp.i"
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/local_optimization_library && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/src/local_optimization_library/local_optimization.cpp > CMakeFiles/default_local_optimization.dir/local_optimization.cpp.i

src/local_optimization_library/CMakeFiles/default_local_optimization.dir/local_optimization.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/default_local_optimization.dir/local_optimization.cpp.s"
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/local_optimization_library && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/src/local_optimization_library/local_optimization.cpp -o CMakeFiles/default_local_optimization.dir/local_optimization.cpp.s

# Object files for target default_local_optimization
default_local_optimization_OBJECTS = \
"CMakeFiles/default_local_optimization.dir/local_optimization.cpp.o"

# External object files for target default_local_optimization
default_local_optimization_EXTERNAL_OBJECTS =

../shared-libs/libdefault_local_optimization.so: src/local_optimization_library/CMakeFiles/default_local_optimization.dir/local_optimization.cpp.o
../shared-libs/libdefault_local_optimization.so: src/local_optimization_library/CMakeFiles/default_local_optimization.dir/build.make
../shared-libs/libdefault_local_optimization.so: src/local_optimization_library/CMakeFiles/default_local_optimization.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module ../../../shared-libs/libdefault_local_optimization.so"
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/local_optimization_library && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/default_local_optimization.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/local_optimization_library/CMakeFiles/default_local_optimization.dir/build: ../shared-libs/libdefault_local_optimization.so

.PHONY : src/local_optimization_library/CMakeFiles/default_local_optimization.dir/build

src/local_optimization_library/CMakeFiles/default_local_optimization.dir/clean:
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/local_optimization_library && $(CMAKE_COMMAND) -P CMakeFiles/default_local_optimization.dir/cmake_clean.cmake
.PHONY : src/local_optimization_library/CMakeFiles/default_local_optimization.dir/clean

src/local_optimization_library/CMakeFiles/default_local_optimization.dir/depend:
	cd /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/src/local_optimization_library /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/local_optimization_library /home/bergolho/Documentos/Github/Network-Generator/cco-project-with-localopt/build/src/local_optimization_library/CMakeFiles/default_local_optimization.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/local_optimization_library/CMakeFiles/default_local_optimization.dir/depend

