set(SOURCE_FILES_CONFIG config.cpp config.h)

ADD_LIBRARY(config STATIC ${SOURCE_FILES_CONFIG})

set_property(TARGET config PROPERTY POSITION_INDEPENDENT_CODE ON)