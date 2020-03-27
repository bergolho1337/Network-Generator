SET(LIB_TYPE MODULE) #set the lib type

############## DEFAULT CLOUD SPHERE LIBRARY ##############################
ADD_LIBRARY(sphere_cloud_generator ${LIB_TYPE} sphere.cpp)
TARGET_LINK_LIBRARIES(sphere_cloud_generator)
##################################################################################

############## DEFAULT CLOUD HEXAEDRON LIBRARY ##############################
ADD_LIBRARY(hexaedron_cloud_generator ${LIB_TYPE} hexaedron.cpp)
TARGET_LINK_LIBRARIES(hexaedron_cloud_generator)
##################################################################################

############## DEFAULT CLOUD PARABOLOID LIBRARY ##############################
ADD_LIBRARY(paraboloid_cloud_generator ${LIB_TYPE} paraboloid.cpp)
TARGET_LINK_LIBRARIES(paraboloid_cloud_generator)
##################################################################################

############## DEFAULT CLOUD SQUARE LIBRARY ##############################
ADD_LIBRARY(square_cloud_generator ${LIB_TYPE} square.cpp)
TARGET_LINK_LIBRARIES(square_cloud_generator)
##################################################################################