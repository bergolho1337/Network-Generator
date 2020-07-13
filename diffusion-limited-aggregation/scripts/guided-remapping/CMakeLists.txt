cmake_minimum_required(VERSION 3.0.1)
project(Guided_Remapping)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/stl-reader)
ADD_SUBDIRECTORY(src/graph)
ADD_SUBDIRECTORY(src/utils)

ADD_EXECUTABLE(Guided_Remapping src/main.cpp)

ADD_DEPENDENCIES(Guided_Remapping reader graph utils)

TARGET_LINK_LIBRARIES(Guided_Remapping ${VTK_LIBRARIES} reader graph utils dl m)

