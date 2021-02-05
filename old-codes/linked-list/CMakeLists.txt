cmake_minimum_required(VERSION 3.0.1)
project(LinkedList)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

ADD_SUBDIRECTORY(src/linked-list)

ADD_EXECUTABLE(LinkedList src/main.cpp)

ADD_DEPENDENCIES(LinkedList linkedlist)

TARGET_LINK_LIBRARIES(LinkedList linkedlist dl m)

