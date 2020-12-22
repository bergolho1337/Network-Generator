SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT PRUNING LIBRARY ######################################
ADD_LIBRARY(default_pruning ${LIB_TYPE} pruning.cpp)
TARGET_LINK_LIBRARIES(default_pruning)
##################################################################################