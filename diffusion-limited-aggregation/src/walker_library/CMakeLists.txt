SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT SQUARE WALKER LIBRARY ########################################
ADD_LIBRARY(square_walker ${LIB_TYPE} square.cpp)
TARGET_LINK_LIBRARIES(square_walker)
#####################################################################################
