
#########################################################################################
# CDF
#########################################################################################

#set paths where CDF can be found
set(CDF_INC_HINTS "/opt/local/include" "/usr/include" "/usr/local/include")
set(CDF_LIB_HINTS "/opt/local/lib" "/usr/lib" "/usr/local/lib")

# Try to find include direcory
find_path(CDF_INCLUDES cdf.h
        HINTS ${CDF_INC_HINTS})

# Choose to find the static or the shared version
if(CDF_USE_STATIC_LIBS)
    find_library(CDF_LIBRARY NAMES ${CMAKE_STATIC_LIBRARY_PREFIX}cdf${CMAKE_STATIC_LIBRARY_SUFFIX}
                 HINTS ${CDF_LIB_HINTS})
else()
    find_library(CDF_LIBRARY NAMES ${CMAKE_SHARED_LIBRARY_PREFIX}cdf${CMAKE_SHARED_LIBRARY_SUFFIX}
                 HINTS ${CDF_LIB_HINTS})
endif(CDF_USE_STATIC_LIBS)

# Handle the arguments, if CDF_LIBRARY and CDF_INCLUDES is set then CDF_FOUND is set
find_package_handle_standard_args(CDF  DEFAULT_MSG
                                  CDF_LIBRARY CDF_INCLUDES)

if(CDF_FOUND) 
	set(CDF_INCLUDE_DIRS ${CDF_INCLUDES})
	set(CDF_LIBRARIES ${CDF_LIBRARY})
endif(CDF_FOUND)