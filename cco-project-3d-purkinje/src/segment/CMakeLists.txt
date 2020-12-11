set(SOURCE_FILES_SEGMENT segment.h segment.cpp)

ADD_LIBRARY(segment STATIC ${SOURCE_FILES_SEGMENT})

TARGET_LINK_LIBRARIES(segment point utils)

if(NOT MSVC)
    set_property(TARGET segment PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()