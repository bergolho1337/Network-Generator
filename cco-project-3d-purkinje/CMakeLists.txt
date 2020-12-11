cmake_minimum_required(VERSION 3.0.1)
project(Cco_3D)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/point)
ADD_SUBDIRECTORY(src/segment)
#ADD_SUBDIRECTORY(src/face-list)
ADD_SUBDIRECTORY(src/options)
ADD_SUBDIRECTORY(src/ini_parser)
ADD_SUBDIRECTORY(src/utils)
ADD_SUBDIRECTORY(src/test)
ADD_SUBDIRECTORY(src/cco)

ADD_SUBDIRECTORY(src/cost_function_library)
ADD_SUBDIRECTORY(src/local_optimization_library)
#ADD_SUBDIRECTORY(src/pruning_library)

ADD_EXECUTABLE(Cco_3D src/main.cpp)

ADD_DEPENDENCIES(Cco_3D costfunction localopt cco point segment options ini_parser test utils)

TARGET_LINK_LIBRARIES(Cco_3D ${VTK_LIBRARIES} costfunction localopt cco point segment options ini_parser test utils dl m)
