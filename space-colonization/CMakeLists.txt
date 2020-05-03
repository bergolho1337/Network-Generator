cmake_minimum_required(VERSION 3.0.1)
project(SpaceColonization)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/options)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/network)

ADD_EXECUTABLE(SpaceColonization src/main.cpp)

ADD_DEPENDENCIES(SpaceColonization options ini_parser utils network)

TARGET_LINK_LIBRARIES(SpaceColonization ${VTK_LIBRARIES} options ini_parser utils network dl m)

