set(SOURCE_FILES_SEGMENT segment-list.h segment-list.cpp)

ADD_LIBRARY(segmentlist STATIC ${SOURCE_FILES_SEGMENT})

if(NOT MSVC)
    set_property(TARGET segmentlist PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()