#############################################################################
#
# This file is part of the libPR software.
# Copyright (C) 2017 by MIS lab (UPJV). All rights reserved.
#
# See http://mis.u-picardie.fr/~g-caron/fr/index.php?page=7 for more information.
#
# This software was developed at:
# MIS - UPJV
# 33 rue Saint-Leu
# 80039 AMIENS CEDEX
# France
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# Description:
# CMake package config file for libPR. Freely inspired from the one of the 
# ViSP library. 
#
# Authors:
# Guillaume Caron
#
#############################################################################
#
# Description:
# CMake package config file for PER.
#
# ** File generated automatically, do not modify **
#
# This file will define the following CMake variables:
#   - PER_INCLUDE_DIRS   : PER and third-party include directories
#   - PER_LIBRARIES      : PER library to link against. Third-party libraries are
#                           linked automatically thanks to cmake export file PERTargets.cmake
#   - PER_VERSION_STRING : Full PER version that is build. Example: "2.10.0"
#   - PER_VERSION_MAJOR  : Major version part of PER_VERSION. Example: "2"
#   - PER_VERSION_MINOR  : Minor version part of PER_VERSION. Example: "10"
#   - PER_VERSION_PATCH  : Patch version part of PER_VERSION. Example: "0"
#
# Advanced variables:
#   - PER_SHARED        : Use PER as shared library
#   - PER_CONFIG_PATH   : Path to this PERConfig.cmake
#   - PER_FIND_QUIETLY  : If set to TRUE turn off messages during configuration
#   - PER_USE_FILE      : File to include to use PER without specific cmake code
#
# Windows specific variables:
#   - PER_STATIC        : If set to TRUE uses PER static library (.lib) rather then dynamic (.dll) 
#
# Typical usage in user project:
#
#   find_package(PER)
#   include_directories(${PER_INCLUDE_DIRS})
#   target_link_libraries(MY_TARGET_NAME ${PER_LIBRARIES})
#
# It is also possible to build your project using PER_USE_FILE.
#
#   find_package(PER)
#   if(PER_FOUND)
#     include(${PER_USE_FILE})
#   endif()
#
# Authors:
# Fabien Spindler
#
#############################################################################

# similar code exist in PERDetectPlatform.cmake
if(MSVC)
  if(CMAKE_CL_64)
    set(PER_ARCH x64)
  else()
    set(PER_ARCH x86)
  endif()
  if(MSVC_VERSION EQUAL 1400)
    set(PER_RUNTIME vc8)
  elseif(MSVC_VERSION EQUAL 1500)
    set(PER_RUNTIME vc9)
  elseif(MSVC_VERSION EQUAL 1600)
    set(PER_RUNTIME vc10)
  elseif(MSVC_VERSION EQUAL 1700)
    set(PER_RUNTIME vc11)
  elseif(MSVC_VERSION EQUAL 1800)
    set(PER_RUNTIME vc12)
  endif()
elseif(MINGW)
  set(PER_RUNTIME mingw)

  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpmachine
                  OUTPUT_VARIABLE PER_GCC_TARGET_MACHINE
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(PER_GCC_TARGET_MACHINE MATCHES "64")
    set(MINGW64 1)
    set(PER_ARCH x64)
  else()
    set(PER_ARCH x86)
  endif()
endif()

if(CMAKE_VERSION VERSION_GREATER 2.6.2)
  unset(PER_CONFIG_PATH CACHE)
endif()

if(NOT PER_FIND_QUIETLY)
  message(STATUS "PER ARCH: ${PER_ARCH}")
  message(STATUS "PER RUNTIME: ${PER_RUNTIME}")
endif()

get_filename_component(PER_CONFIG_PATH "${CMAKE_CURRENT_LIST_FILE}" PATH CACHE)
if(PER_RUNTIME AND PER_ARCH)
  if(NOT DEFINED PER_STATIC AND EXISTS "${PER_CONFIG_PATH}/${PER_ARCH}/${PER_RUNTIME}/lib/PERConfig.cmake")
    set(PER_LIB_PATH "${PER_CONFIG_PATH}/${PER_ARCH}/${PER_RUNTIME}/lib")
  elseif(NOT DEFINED PER_STATIC AND EXISTS "${PER_CONFIG_PATH}/${PER_ARCH}/${PER_RUNTIME}/staticlib/PERConfig.cmake")
    set(PER_LIB_PATH "${PER_CONFIG_PATH}/${PER_ARCH}/${PER_RUNTIME}/staticlib")
  elseif(PER_STATIC AND EXISTS "${PER_CONFIG_PATH}/${PER_ARCH}/${PER_RUNTIME}/staticlib/PERConfig.cmake")
    set(PER_LIB_PATH "${PER_CONFIG_PATH}/${PER_ARCH}/${PER_RUNTIME}/staticlib")
  elseif(PER_STATIC EXISTS "${PER_CONFIG_PATH}/${PER_ARCH}/${PER_RUNTIME}/lib/PERConfig.cmake")
    set(PER_LIB_PATH "${PER_CONFIG_PATH}/${PER_ARCH}/${PER_RUNTIME}/lib")
  endif()
endif()

if(PER_LIB_PATH AND EXISTS "${PER_LIB_PATH}/PERConfig.cmake")
  include("${PER_LIB_PATH}/PERConfig.cmake")

  set(PER_FOUND TRUE CACHE BOOL "" FORCE)

  if(NOT PER_FIND_QUIETLY)
    message(STATUS "Found PER ${PER_VERSION} in ${PER_LIB_PATH}")
    if(NOT PER_LIB_PATH MATCHES "/staticlib")
      get_filename_component(_PER_LIB_PATH "${PER_LIB_PATH}/../bin" ABSOLUTE)
      file(TO_NATIVE_PATH "${_PER_LIB_PATH}" _PER_LIB_PATH)
      message(STATUS "You might need to add ${_PER_LIB_PATH} to your PATH to be able to run your applications.")
    endif()
  endif()
else()
  if(NOT PER_FIND_QUIETLY)
    message(WARNING
"Found PER for Windows but it has no binaries compatible with your configuration.
You should manually point CMake variable PER_DIR to your build of PER library."
    )
  endif()
  set(PER_FOUND FALSE CACHE BOOL "" FORCE)
endif()
