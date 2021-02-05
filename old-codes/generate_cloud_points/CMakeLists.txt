cmake_minimum_required(VERSION 3.0.1)
project(CloudGenerator)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
    endif ()
else ()
    set (CMAKE_CXX_STANDARD 11)
endif ()

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/shared-libs)

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fno-strict-aliasing -Wall -Wno-unused-function")

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_SUBDIRECTORY(src/generator)
ADD_SUBDIRECTORY(src/config)
ADD_SUBDIRECTORY(src/utils)

ADD_SUBDIRECTORY(src/cloud_library)

ADD_EXECUTABLE(CloudGenerator src/main.cpp)

ADD_DEPENDENCIES(CloudGenerator generator config utils)

TARGET_LINK_LIBRARIES(CloudGenerator ${VTK_LIBRARIES} generator config utils dl m)