set(SOURCE_FILES_GENERATOR generator.cpp generator.h)

ADD_LIBRARY(generator STATIC ${SOURCE_FILES_GENERATOR})

set_property(TARGET generator PROPERTY POSITION_INDEPENDENT_CODE ON)