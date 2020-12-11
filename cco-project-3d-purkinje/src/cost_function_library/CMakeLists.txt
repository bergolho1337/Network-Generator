set(SOURCE_FILES_COST_FUNCTION cost_function.h cost_function.cpp)

ADD_LIBRARY(costfunction STATIC ${SOURCE_FILES_COST_FUNCTION})

TARGET_LINK_LIBRARIES(costfunction cco)

if(NOT MSVC)
    set_property(TARGET costfunction PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()