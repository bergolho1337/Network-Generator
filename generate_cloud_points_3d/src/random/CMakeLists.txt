set(SOURCE_FILES_RANDOM random.h random.cpp
    dSFMT.cpp dSFMT.h dSFMT-params.h dSFMT-params19937.h)

ADD_LIBRARY(random STATIC ${SOURCE_FILES_RANDOM})

if(NOT MSVC)
    set_property(TARGET random PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()
