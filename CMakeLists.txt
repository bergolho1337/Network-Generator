cmake_minimum_required(VERSION 3.0.1)
project(NetworkGenerator)

if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
    endif ()
else ()
    set (CMAKE_CXX_STANDARD 11)
endif ()

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin )

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -fno-strict-aliasing -Wall -Wno-unused-function")

FIND_PACKAGE(OpenMP)
IF (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    #message("${CMAKE_CXX_FLAGS}")
ENDIF(OPENMP_FOUND)

INCLUDE_DIRECTORIES(include)

FILE(GLOB SOURCES "src/*.cpp")

ADD_EXECUTABLE(NetworkGenerator ${SOURCES}) 

#ADD_DEPENDENCIES(FDMMonodomain1D utils solvers hashes ini_parser string )

TARGET_LINK_LIBRARIES(NetworkGenerator m)

#ADD_SUBDIRECTORY(src/models_library)

#SET(COMMON_SRC ../libraries_common/config_helpers.c ../libraries_common/config_helpers.h)

#ADD_SUBDIRECTORY(src/stimuli_library)
###


