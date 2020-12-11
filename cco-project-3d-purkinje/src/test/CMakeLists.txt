set(SOURCE_FILES_TEST test.h test.cpp)

ADD_LIBRARY(test STATIC ${SOURCE_FILES_TEST})

TARGET_LINK_LIBRARIES(test cco point segment utils)

if(NOT MSVC)
    set_property(TARGET test PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()