set(SOURCE_FILES_READER reader.h reader.cpp
                graph.h graph.cpp
                cloud.h cloud.cpp)

ADD_LIBRARY(reader STATIC ${SOURCE_FILES_READER})

set_property(TARGET reader PROPERTY POSITION_INDEPENDENT_CODE ON)
