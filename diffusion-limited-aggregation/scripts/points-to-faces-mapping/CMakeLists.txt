cmake_minimum_required(VERSION 3.0.1)
project(Points_To_Faces_Mapping)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/stl-reader)
ADD_SUBDIRECTORY(src/utils)

ADD_EXECUTABLE(Points_To_Faces_Mapping src/main.cpp)

ADD_DEPENDENCIES(Points_To_Faces_Mapping reader utils)

TARGET_LINK_LIBRARIES(Points_To_Faces_Mapping ${VTK_LIBRARIES} reader utils dl m)

