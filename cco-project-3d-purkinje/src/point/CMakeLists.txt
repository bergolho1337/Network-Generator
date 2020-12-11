set(SOURCE_FILES_POINT point.h point.cpp)

ADD_LIBRARY(point STATIC ${SOURCE_FILES_POINT})

#TARGET_LINK_LIBRARIES(cco costfunction localopt utils)

if(NOT MSVC)
    set_property(TARGET point PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()