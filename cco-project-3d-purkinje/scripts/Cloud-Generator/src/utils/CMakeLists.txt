set(SOURCE_FILES_OPTIONS utils.h utils.cpp)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_OPTIONS})

if(NOT MSVC)
    set_property(TARGET utils PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()