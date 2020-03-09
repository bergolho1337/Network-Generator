cmake_minimum_required(VERSION 3.0.1)
project(Point-Remapper-v2)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/stl-reader)
ADD_SUBDIRECTORY(src/graph)
#ADD_SUBDIRECTORY(src/ini_parser)
#ADD_SUBDIRECTORY(src/utils)

ADD_EXECUTABLE(Point-Remapper-v2 src/main.cpp)

ADD_DEPENDENCIES(Point-Remapper-v2 reader graph)

TARGET_LINK_LIBRARIES(Point-Remapper-v2 ${VTK_LIBRARIES} reader graph dl m)

