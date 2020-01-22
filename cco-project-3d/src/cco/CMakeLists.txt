set(SOURCE_FILES_CCO cco.h cco.cpp
                     cco_helper.h cco_helper.cpp)

ADD_LIBRARY(cco STATIC ${SOURCE_FILES_CCO})

TARGET_LINK_LIBRARIES(cco pointlist segmentlist utils)

if(NOT MSVC)
    set_property(TARGET cco PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()