SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT NEWTON-COTES LIBRARY ######################################
ADD_LIBRARY(default_cost_function ${LIB_TYPE} cost_function.cpp)
TARGET_LINK_LIBRARIES(default_cost_function)
##################################################################################