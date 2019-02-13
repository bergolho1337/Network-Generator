cmake_minimum_required(VERSION 3.0.1)
project(Cco)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

ADD_SUBDIRECTORY(src/point-list)
ADD_SUBDIRECTORY(src/segment-list)
ADD_SUBDIRECTORY(src/options)
ADD_SUBDIRECTORY(src/cco)

ADD_EXECUTABLE(Cco src/main.cpp)

ADD_DEPENDENCIES(Cco pointlist segmentlist options cco)

TARGET_LINK_LIBRARIES(Cco pointlist segmentlist options cco dl m)

