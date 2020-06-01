set(SOURCE_FILES_NETWORK tree.h tree.cpp
                        leaf.h leaf.cpp
                        branch.h branch.cpp
                        point.h point.cpp)

ADD_LIBRARY(network STATIC ${SOURCE_FILES_NETWORK})

if(NOT MSVC)
    set_property(TARGET network PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()