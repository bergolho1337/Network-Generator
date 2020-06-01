SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT MINIMIZE VOLUME LIBRARY ######################################
ADD_LIBRARY(minimize_volume ${LIB_TYPE} minimize_volume.cpp)
TARGET_LINK_LIBRARIES(minimize_volume)
#####################################################################################

############## DEFAULT MINIMIZE ACTIVATION TIME LIBRARY #############################
ADD_LIBRARY(minimize_activation_time ${LIB_TYPE} minimize_activation_time.cpp)
TARGET_LINK_LIBRARIES(minimize_activation_time)
#####################################################################################

############## DEFAULT MINIMIZE CUSTOM FUNCTION LIBRARY #############################
ADD_LIBRARY(minimize_custom_function ${LIB_TYPE} minimize_custom_function.cpp)
TARGET_LINK_LIBRARIES(minimize_custom_function)
#####################################################################################