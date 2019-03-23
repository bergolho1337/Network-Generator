SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT CLOUD CIRCLE LIBRARY ##############################
ADD_LIBRARY(circle_cloud_generator ${LIB_TYPE} circle.cpp)
TARGET_LINK_LIBRARIES(circle_cloud_generator)
##################################################################################