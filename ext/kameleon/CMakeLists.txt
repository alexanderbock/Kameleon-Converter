#########################################################################################
# General settings
#########################################################################################

# top level CMakeLists.txt (kameleon-plus/)
cmake_minimum_required (VERSION 2.8) 
project (Kameleon)

#Set version number
set(Kameleon_VERSION_MAJOR "6")
set(Kameleon_VERSION_MINOR "0")
set(Kameleon_VERSION_REVISION "0")
set(KAMELEON_PLUS_VERSION ${Kameleon_VERSION_MAJOR}.${Kameleon_VERSION_MINOR}.${Kameleon_VERSION_REVISION})
message(STATUS "Kameleon_plus_version " ${KAMELEON_PLUS_VERSION})

# Set Kameleon paths
set(KAMELEON_BASE_DIR "${PROJECT_SOURCE_DIR}")
set(KAMELEON_SRC_DIR "${KAMELEON_BASE_DIR}/src")
set(KAMELEON_EXT_DIR "${KAMELEON_BASE_DIR}/ext")

# Link Kameleon as library (do not build wrappers and examples)
option(KAMELEON_LIBRARY_ONLY "Build with Kameleon as library only" OFF)

# add path to external libraries
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${KAMELEON_EXT_DIR})

#########################################################################################
# External Third-party software options
#########################################################################################

#boost settings
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME ON) 

# If libraries needs to be built, make them shared
option(BUILD_SHARED_LIBS "Build Shared Libraries" ON)

# HDF5 options
option(HDF5_USE_STATIC_LIBRARIES "Use HDF5 static libraries" OFF)
option(KAMELEON_USE_HDF5 "Build Kameleon with HDF5 library" ON)

# CDF options
option(CDF_USE_STATIC_LIBS "Use CDF static libraries" OFF)
option(CDF_BUILD_ZLIB "Build CDF with zlib" ON)

#########################################################################################
# External Third-party software
#########################################################################################

# Boost
find_package(Boost REQUIRED)
set(CCMC_DEPENDENT_LIBS ${CCMC_DEPENDENT_LIBS} ${Boost_LIBRARIES})

# HDF5
find_package(HDF5 COMPONENTS C CXX)
if(KAMELEON_USE_HDF5)
	if(HDF5_FOUND)
		message(STATUS "HDF5 Found: Building Kameleon with HDF 5")
		add_definitions(-DHAVE_HDF5)
		set(CCMC_DEPENDENT_LIBS ${CCMC_DEPENDENT_LIBS} ${HDF5_LIBRARIES})
	else(HDF5_FOUND)
		message(STATUS "Buildning HDF5 locally")
		set(HDF5_ROOT ${KAMELEON_EXT_DIR}/hdf5/hdf5-1.8.12)
		set(HDF5_EXTERNALLY_CONFIGURED 1)
		set(HDF5_BUILD_CPP_LIB ON CACHE BOOL "Build HDF5 C++ Library" FORCE)
		add_subdirectory(${HDF5_ROOT})
		set(HDF5_INCLUDE_DIR ${HDF5_ROOT} ${HDF5_ROOT}/src ${HDF5_ROOT}/c++/src)
		SET(HDF5_LIBRARIES hdf5)
		SET(HDF5_LIBRARIES hdf5_cpp)
		set(CCMC_DEPENDENT_LIBS ${CCMC_DEPENDENT_LIBS} ${HDF5_LIBRARIES})
	endif(HDF5_FOUND)
else()
	message(STATUS "Not Building HDF5 - CDF files only supported")
endif(KAMELEON_USE_HDF5)

# CDF
find_package(CDF)
if(NOT CDF_FOUND)
	message(STATUS "Buildning CDF locally")
	add_subdirectory(${KAMELEON_EXT_DIR}/cdf)
endif(NOT CDF_FOUND)
set(CCMC_DEPENDENT_LIBS ${CCMC_DEPENDENT_LIBS} ${CDF_LIBRARIES})

#########################################################################################
# Configuration
#########################################################################################
configure_file (
	${PROJECT_SOURCE_DIR}/kameleon-plus-Config.h.in
	${PROJECT_SOURCE_DIR}/src/ccmc/Kameleon-plus-Config.h
)
include_directories("${PROJECT_SOURCE_DIR}")

#########################################################################################
# Executable
#########################################################################################
# CCMC
set(CCMC_ROOT_DIR "${KAMELEON_SRC_DIR}/ccmc")
add_subdirectory(${CCMC_ROOT_DIR})

#########################################################################################
# C++ Examples
#########################################################################################
set(KAMELEON_EXAMPLES_DIR "${PROJECT_SOURCE_DIR}/examples")
set(EXAMPLES_ROOT_DIR "${KAMELEON_BASE_DIR}/src/example/c++")
if(NOT KAMELEON_LIBRARY_ONLY)
	add_subdirectory(${EXAMPLES_ROOT_DIR})
endif(NOT KAMELEON_LIBRARY_ONLY)

