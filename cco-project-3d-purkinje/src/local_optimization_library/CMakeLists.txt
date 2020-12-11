set(SOURCE_FILES_COST_FUNCTION local_optimization.h local_optimization.cpp)

ADD_LIBRARY(localopt STATIC ${SOURCE_FILES_COST_FUNCTION})

TARGET_LINK_LIBRARIES(localopt cco)

if(NOT MSVC)
    set_property(TARGET localopt PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()