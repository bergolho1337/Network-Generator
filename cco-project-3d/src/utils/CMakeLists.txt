set(SOURCE_FILES_UTILS utils.h utils.cpp
    dSFMT.cpp dSFMT.h dSFMT-params.h dSFMT-params19937.h)

ADD_LIBRARY(utils STATIC ${SOURCE_FILES_UTILS})

if(NOT MSVC)
    set_property(TARGET utils PROPERTY POSITION_INDEPENDENT_CODE ON)
endif()