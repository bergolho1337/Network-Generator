set(SOURCE_FILES_CCO cco.h cco.cpp
                     pmj.h pmj.cpp
                     cloud.h cloud.cpp
                     helper.h helper.cpp
                     constants.h)

ADD_LIBRARY(cco STATIC ${SOURCE_FILES_CCO})

TARGET_LINK_LIBRARIES(cco costfunction localopt point segment utils)

if(NOT MSVC)
    set_property(TARGET cco PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()