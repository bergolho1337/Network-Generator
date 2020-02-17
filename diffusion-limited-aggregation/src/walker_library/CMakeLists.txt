SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT SQUARE WALKER LIBRARY ########################################
ADD_LIBRARY(square_walker ${LIB_TYPE} square.cpp)
TARGET_LINK_LIBRARIES(square_walker)
#####################################################################################

############## DEFAULT CIRCLE WALKER LIBRARY ########################################
ADD_LIBRARY(circle_walker ${LIB_TYPE} circle.cpp)
TARGET_LINK_LIBRARIES(circle_walker)
#####################################################################################

############## DEFAULT BOX WALKER LIBRARY ########################################
ADD_LIBRARY(box_walker ${LIB_TYPE} box.cpp)
TARGET_LINK_LIBRARIES(box_walker)
#####################################################################################

############## DEFAULT SPHERE WALKER LIBRARY ########################################
ADD_LIBRARY(sphere_walker ${LIB_TYPE} sphere.cpp)
TARGET_LINK_LIBRARIES(sphere_walker)
#####################################################################################

############## DEFAULT CUSTOM LIBRARY ###############################################
ADD_LIBRARY(custom_walker ${LIB_TYPE} custom.cpp)
TARGET_LINK_LIBRARIES(custom_walker)
#####################################################################################
