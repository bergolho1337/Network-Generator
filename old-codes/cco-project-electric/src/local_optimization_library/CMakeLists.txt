SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT NEWTON-COTES LIBRARY ######################################
ADD_LIBRARY(default_local_optimization ${LIB_TYPE} local_optimization.cpp)
TARGET_LINK_LIBRARIES(default_local_optimization)
##################################################################################