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
CMAKE_SOURCE_DIR = /home/berg/Github/Network-Generator/cco-project/backup

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/berg/Github/Network-Generator/cco-project/backup/build

# Include any dependencies generated for this target.
include CMakeFiles/Cco.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/Cco.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/Cco.dir/flags.make

CMakeFiles/Cco.dir/src/cco.cpp.o: CMakeFiles/Cco.dir/flags.make
CMakeFiles/Cco.dir/src/cco.cpp.o: ../src/cco.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/berg/Github/Network-Generator/cco-project/backup/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/Cco.dir/src/cco.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Cco.dir/src/cco.cpp.o -c /home/berg/Github/Network-Generator/cco-project/backup/src/cco.cpp

CMakeFiles/Cco.dir/src/cco.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cco.dir/src/cco.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/berg/Github/Network-Generator/cco-project/backup/src/cco.cpp > CMakeFiles/Cco.dir/src/cco.cpp.i

CMakeFiles/Cco.dir/src/cco.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cco.dir/src/cco.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/berg/Github/Network-Generator/cco-project/backup/src/cco.cpp -o CMakeFiles/Cco.dir/src/cco.cpp.s

CMakeFiles/Cco.dir/src/graph.cpp.o: CMakeFiles/Cco.dir/flags.make
CMakeFiles/Cco.dir/src/graph.cpp.o: ../src/graph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/berg/Github/Network-Generator/cco-project/backup/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/Cco.dir/src/graph.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Cco.dir/src/graph.cpp.o -c /home/berg/Github/Network-Generator/cco-project/backup/src/graph.cpp

CMakeFiles/Cco.dir/src/graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cco.dir/src/graph.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/berg/Github/Network-Generator/cco-project/backup/src/graph.cpp > CMakeFiles/Cco.dir/src/graph.cpp.i

CMakeFiles/Cco.dir/src/graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cco.dir/src/graph.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/berg/Github/Network-Generator/cco-project/backup/src/graph.cpp -o CMakeFiles/Cco.dir/src/graph.cpp.s

CMakeFiles/Cco.dir/src/main.cpp.o: CMakeFiles/Cco.dir/flags.make
CMakeFiles/Cco.dir/src/main.cpp.o: ../src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/berg/Github/Network-Generator/cco-project/backup/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/Cco.dir/src/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Cco.dir/src/main.cpp.o -c /home/berg/Github/Network-Generator/cco-project/backup/src/main.cpp

CMakeFiles/Cco.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cco.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/berg/Github/Network-Generator/cco-project/backup/src/main.cpp > CMakeFiles/Cco.dir/src/main.cpp.i

CMakeFiles/Cco.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cco.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/berg/Github/Network-Generator/cco-project/backup/src/main.cpp -o CMakeFiles/Cco.dir/src/main.cpp.s

CMakeFiles/Cco.dir/src/options.cpp.o: CMakeFiles/Cco.dir/flags.make
CMakeFiles/Cco.dir/src/options.cpp.o: ../src/options.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/berg/Github/Network-Generator/cco-project/backup/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/Cco.dir/src/options.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Cco.dir/src/options.cpp.o -c /home/berg/Github/Network-Generator/cco-project/backup/src/options.cpp

CMakeFiles/Cco.dir/src/options.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Cco.dir/src/options.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/berg/Github/Network-Generator/cco-project/backup/src/options.cpp > CMakeFiles/Cco.dir/src/options.cpp.i

CMakeFiles/Cco.dir/src/options.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Cco.dir/src/options.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/berg/Github/Network-Generator/cco-project/backup/src/options.cpp -o CMakeFiles/Cco.dir/src/options.cpp.s

# Object files for target Cco
Cco_OBJECTS = \
"CMakeFiles/Cco.dir/src/cco.cpp.o" \
"CMakeFiles/Cco.dir/src/graph.cpp.o" \
"CMakeFiles/Cco.dir/src/main.cpp.o" \
"CMakeFiles/Cco.dir/src/options.cpp.o"

# External object files for target Cco
Cco_EXTERNAL_OBJECTS =

../bin/Cco: CMakeFiles/Cco.dir/src/cco.cpp.o
../bin/Cco: CMakeFiles/Cco.dir/src/graph.cpp.o
../bin/Cco: CMakeFiles/Cco.dir/src/main.cpp.o
../bin/Cco: CMakeFiles/Cco.dir/src/options.cpp.o
../bin/Cco: CMakeFiles/Cco.dir/build.make
../bin/Cco: /usr/local/lib/libvtkDomainsChemistryOpenGL2-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersFlowPaths-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersGeneric-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersHyperTree-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersParallelImaging-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersPoints-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersProgrammable-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersSMP-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersSelection-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersTexture-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersVerdict-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkverdict-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkGeovisCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkproj4-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOAMR-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOEnSight-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOExodus-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOExport-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingGL2PSOpenGL2-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkgl2ps-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOImport-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOInfovis-7.1.so.1
../bin/Cco: /usr/local/lib/libvtklibxml2-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOLSDyna-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOMINC-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOMovie-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkoggtheora-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOPLY-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOParallel-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkjsoncpp-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOParallelXML-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOSQL-7.1.so.1
../bin/Cco: /usr/local/lib/libvtksqlite-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOTecplotTable-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOVideo-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingMorphological-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingStatistics-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingStencil-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkInteractionImage-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingContextOpenGL2-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingImage-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingLOD-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingVolumeOpenGL2-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkViewsContext2D-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkViewsInfovis-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkDomainsChemistry-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersAMR-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersParallel-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkexoIIc-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOGeometry-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIONetCDF-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkNetCDF_cxx-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkNetCDF-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkhdf5_hl-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkhdf5-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkParallelCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOLegacy-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingOpenGL2-7.1.so.1
../bin/Cco: /usr/lib64/libSM.so
../bin/Cco: /usr/lib64/libICE.so
../bin/Cco: /usr/lib64/libX11.so
../bin/Cco: /usr/lib64/libXext.so
../bin/Cco: /usr/lib64/libXt.so
../bin/Cco: /usr/local/lib/libvtkglew-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingMath-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkChartsCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingContext2D-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersImaging-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkInfovisLayout-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkInfovisCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkViewsCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkInteractionWidgets-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersHybrid-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingGeneral-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingSources-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersModeling-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingHybrid-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOImage-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkDICOMParser-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkmetaio-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkpng-7.1.so.1
../bin/Cco: /usr/local/lib/libvtktiff-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkjpeg-7.1.so.1
../bin/Cco: /usr/lib64/libm.so
../bin/Cco: /usr/local/lib/libvtkInteractionStyle-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersExtraction-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersStatistics-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingFourier-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkalglib-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingAnnotation-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingColor-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingVolume-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkImagingCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOXML-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOXMLParser-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkIOCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkexpat-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingLabel-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingFreeType-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkRenderingCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkCommonColor-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersGeometry-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersSources-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersGeneral-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkCommonComputationalGeometry-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkFiltersCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkCommonExecutionModel-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkCommonDataModel-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkCommonTransforms-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkCommonMisc-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkCommonMath-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkCommonSystem-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkCommonCore-7.1.so.1
../bin/Cco: /usr/local/lib/libvtksys-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkfreetype-7.1.so.1
../bin/Cco: /usr/local/lib/libvtkzlib-7.1.so.1
../bin/Cco: CMakeFiles/Cco.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/berg/Github/Network-Generator/cco-project/backup/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable ../bin/Cco"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Cco.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/Cco.dir/build: ../bin/Cco

.PHONY : CMakeFiles/Cco.dir/build

CMakeFiles/Cco.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/Cco.dir/cmake_clean.cmake
.PHONY : CMakeFiles/Cco.dir/clean

CMakeFiles/Cco.dir/depend:
	cd /home/berg/Github/Network-Generator/cco-project/backup/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/berg/Github/Network-Generator/cco-project/backup /home/berg/Github/Network-Generator/cco-project/backup /home/berg/Github/Network-Generator/cco-project/backup/build /home/berg/Github/Network-Generator/cco-project/backup/build /home/berg/Github/Network-Generator/cco-project/backup/build/CMakeFiles/Cco.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/Cco.dir/depend
