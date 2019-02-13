set(SOURCE_FILES_LINKEDLIST linked-list.h linked-list.cpp)

ADD_LIBRARY(linkedlist STATIC ${SOURCE_FILES_LINKEDLIST})

if(NOT MSVC)
    set_property(TARGET linkedlist PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()