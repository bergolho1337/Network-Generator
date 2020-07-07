set(SOURCE_FILES_TREE tree.h tree.cpp
                    walker.h walker.cpp
                    segment.h segment.cpp
                    network.h network.cpp)

ADD_LIBRARY(tree STATIC ${SOURCE_FILES_TREE})

if(NOT MSVC)
    set_property(TARGET tree PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()