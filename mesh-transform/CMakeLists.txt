cmake_minimum_required(VERSION 3.0.1)
project(MeshTransform)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/vtk_utils)

ADD_EXECUTABLE(MeshTransform src/main.cpp)

ADD_DEPENDENCIES(MeshTransform vtk_utils)

TARGET_LINK_LIBRARIES(MeshTransform ${VTK_LIBRARIES} vtk_utils dl m)
