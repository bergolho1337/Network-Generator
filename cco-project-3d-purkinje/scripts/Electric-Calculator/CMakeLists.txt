cmake_minimum_required(VERSION 3.0.1)
project(PurkinjeError)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/reader)
ADD_SUBDIRECTORY(src/graph)
ADD_SUBDIRECTORY(src/utils)

ADD_EXECUTABLE(PurkinjeError src/main.cpp)

ADD_DEPENDENCIES(PurkinjeError reader graph utils)

TARGET_LINK_LIBRARIES(PurkinjeError ${VTK_LIBRARIES} reader graph utils dl m)
