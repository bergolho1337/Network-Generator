set(SOURCE_FILES_CCO cco.h cco.cpp)

ADD_LIBRARY(cco STATIC ${SOURCE_FILES_CCO})

TARGET_LINK_LIBRARIES(cco pointlist segmentlist)

if(NOT MSVC)
    set_property(TARGET cco PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()