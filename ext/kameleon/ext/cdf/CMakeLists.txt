# core kameleon library (kameleon-plus/src/ccmc)

#########################################################################################
# General settings
#########################################################################################
if (NOT CDF_ROOT_DIR)
    set(CDF_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cdf35_0-dist)
endif ()

#########################################################################################
# Files and group settings
#########################################################################################
file(GLOB CDF_SOURCE
    ${CDF_ROOT_DIR}/src/lib/*.cpp ${CDF_ROOT_DIR}/src/lib/*.c
)
file(GLOB CDF_HEADER
    ${CDF_ROOT_DIR}/src/include/*.h ${CDF_ROOT_DIR}/src/include/*.inc
)

if(NOT WIN32)
	list(REMOVE_ITEM CDF_SOURCE ${CDF_ROOT_DIR}/src/lib/libmain.c)
endif(NOT WIN32)

source_group(Cdf FILES ${CDF_SOURCE})  
source_group(Cdf FILES ${CDF_HEADER})  

set(ALL_SOURCE ${CDF_SOURCE})
set(ALL_HEADER ${CDF_HEADER})

option(CDF_USE_ZLIB "Use ZLib" ON)
if(CDF_BUILD_ZLIB)
	file(GLOB ZLIB_SOURCE
    	${CDF_ROOT_DIR}/src/lib/zlib/*.cpp ${CDF_ROOT_DIR}/src/lib/zlib/*.c
	)
	file(GLOB ZLIB_HEADER
    	${CDF_ROOT_DIR}/src/lib/zlib/*.h 
	)

	source_group(Zlib FILES ${ZLIB_HEADER}) 
	source_group(Zlib FILES ${ZLIB_SOURCE}) 

	set(ALL_SOURCE ${ALL_SOURCE} ${ZLIB_SOURCE})
endif(CDF_BUILD_ZLIB)

add_definitions(-Di386 -D__MACH__ -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -D_LARGEFILE_SOURCE -DHAVE_STDARG_H -DHAVE_UNISTD_H -DZ_PREFIX)
if(__APPLE__)
	add_definitions(-Dmac )
endif(__APPLE__)

#########################################################################################
# Create CDF Library
#########################################################################################
if(BUILD_SHARED_LIBS)
	add_library(cdf SHARED ${ALL_SOURCE} ${ALL_HEADER})
else()
	add_library(cdf STATIC ${ALL_SOURCE} ${ALL_HEADER})
endif(BUILD_SHARED_LIBS)

# CDF
set(CDF_INCLUDE_DIRS ${CDF_ROOT_DIR}/src/include)
set(CDF_LIBRARIES cdf )

include_directories(${CDF_INCLUDE_DIRS})

# hack to make variables accessible from top-level CMakeLists
set(CDF_INCLUDE_DIRS ${CDF_INCLUDE_DIRS} PARENT_SCOPE)
set(CDF_LIBRARIES ${CDF_LIBRARIES} PARENT_SCOPE)




