##################################################################################################
#                                                                                                #
# Kameleon Converter                                                                             #
#                                                                                                #
# Copyright (c) 2016, Alexander Bock                                                             #
# All rights reserved.                                                                           #
#                                                                                                #
# Redistribution and use in source and binary forms, with or without modification, are permitted #
# provided that the following conditions are met:                                                #
#                                                                                                #
# 1. Redistributions of source code must retain the above copyright notice, this list of         #
# conditions and the following disclaimer.                                                       #
#                                                                                                #
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of      #
# conditions and the following disclaimer in the documentation and/or other materials provided   #
# with the distribution.                                                                         #
#                                                                                                #
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to   #
# endorse or promote products derived from this software without specific prior written          #
# permission.                                                                                    #
#                                                                                                #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR #
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY   #
# AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR #
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR            #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR       #
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY   #
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR   #
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE          #
# POSSIBILITY OF SUCH DAMAGE.                                                                    #
##################################################################################################

####################
# General settings #
####################
project(KameleonConverter)
cmake_minimum_required(VERSION 3.0)

if (MSVC)
    if (MSVC_VERSION LESS 1900)
        message(FATAL_ERROR "Visual Studio 2015 required")
    endif ()
endif ()

set_property(GLOBAL PROPERTY USE_FOLDERS On)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER CMake)

add_executable(KameleonConverter main.cpp)

set_property(TARGET KameleonConverter PROPERTY CXX_STANDARD 14)
set_property(TARGET KameleonConverter PROPERTY CXX_STANDARD_REQUIRED On)

add_subdirectory(ext/ghoul)
target_link_libraries(KameleonConverter Ghoul)
get_property(GHOUL_INCLUDE_DIR TARGET Ghoul PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
target_include_directories(KameleonConverter PUBLIC ${GHOUL_INCLUDE_DIR})
get_property(GHOUL_DEFINITIONS TARGET Ghoul PROPERTY INTERFACE_COMPILE_DEFINITIONS)
target_compile_definitions(KameleonConverter PUBLIC ${GHOUL_DEFINITIONS})
set_property(TARGET Lua PROPERTY FOLDER "External")
set_property(TARGET lz4 PROPERTY FOLDER "External")
set_property(TARGET tinyobjloader PROPERTY FOLDER "External")


option(KAMELEON_LIBRARY_ONLY "Build with Kameleon as library only" ON)
if (WIN32)
    option(BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
else ()
    option(BUILD_SHARED_LIBS "Build Shared Libraries" ON)
endif ()
mark_as_advanced(BUILD_SHARED_LIBS)
option(KAMELEON_USE_HDF5 "Kameleon use HDF5" OFF)
set(KAMELEON_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/kameleon)
set(KAMELEON_INCLUDES ${KAMELEON_ROOT_DIR}/src)
add_subdirectory(${KAMELEON_ROOT_DIR})
target_include_directories(KameleonConverter SYSTEM PUBLIC ${KAMELEON_INCLUDES})
target_link_libraries(KameleonConverter ccmc)
set_property(TARGET ccmc PROPERTY FOLDER "External")
if (TARGET cdf)
    set_property(TARGET cdf PROPERTY FOLDER "External")
endif ()
