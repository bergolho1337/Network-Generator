set(SOURCE_FILES_POINT point-list.h point-list.cpp)

ADD_LIBRARY(pointlist STATIC ${SOURCE_FILES_POINT})

if(NOT MSVC)
    set_property(TARGET pointlist PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()