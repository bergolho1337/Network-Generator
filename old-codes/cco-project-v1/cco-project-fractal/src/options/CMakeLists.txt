set(SOURCE_FILES_OPTIONS user_options.h user_options.cpp
                         cost_function_config.h cost_function_config.cpp)

ADD_LIBRARY(options STATIC ${SOURCE_FILES_OPTIONS})

if(NOT MSVC)
    set_property(TARGET options PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()