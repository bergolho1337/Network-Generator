SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT MINIMIZE VOLUME LIBRARY ######################################
ADD_LIBRARY(minimize_volume ${LIB_TYPE} minimize_volume.cpp)
TARGET_LINK_LIBRARIES(minimize_volume)
#####################################################################################
