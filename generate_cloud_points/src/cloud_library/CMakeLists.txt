SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT CLOUD CIRCLE LIBRARY ##############################
ADD_LIBRARY(circle_cloud_generator ${LIB_TYPE} circle.cpp)
TARGET_LINK_LIBRARIES(circle_cloud_generator)
##################################################################################

############## DEFAULT CLOUD SQUARE LIBRARY ##############################
ADD_LIBRARY(square_cloud_generator ${LIB_TYPE} square.cpp)
TARGET_LINK_LIBRARIES(square_cloud_generator)
##################################################################################

############## DEFAULT CLOUD TRIANGLE LIBRARY ##############################
ADD_LIBRARY(triangle_cloud_generator ${LIB_TYPE} triangle.cpp)
TARGET_LINK_LIBRARIES(triangle_cloud_generator)
##################################################################################
