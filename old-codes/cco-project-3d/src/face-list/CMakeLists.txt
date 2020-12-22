set(SOURCE_FILES_FACE face-list.h face-list.cpp)

ADD_LIBRARY(facelist STATIC ${SOURCE_FILES_FACE})

if(NOT MSVC)
    set_property(TARGET facelist PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()