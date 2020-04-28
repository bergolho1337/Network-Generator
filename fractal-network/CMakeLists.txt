cmake_minimum_required(VERSION 3.0.1)
project(FractalTree)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/options)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/fractal)

ADD_EXECUTABLE(FractalTree src/main.cpp)

ADD_DEPENDENCIES(FractalTree options ini_parser utils fractal)

TARGET_LINK_LIBRARIES(FractalTree ${VTK_LIBRARIES} options ini_parser utils fractal dl m)
