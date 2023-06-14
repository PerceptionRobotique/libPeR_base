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
# libPR overall configuration file. Freely inspired from the CMakeList tree of the 
# ViSP library. Detect third party libraries (X11, GTK, ...)
#
# Authors:
# Guillaume Caron
#
#############################################################################

# ----------------------------------------------------------------------------
# Similar code exist in PERConfig.cmake
# ----------------------------------------------------------------------------

if(NOT DEFINED PER_STATIC)
  # look for global setting
  if(NOT DEFINED BUILD_SHARED_LIBS OR BUILD_SHARED_LIBS)
    add_definitions(-DPER_SHARED_EXPORT)
    set(PER_STATIC OFF)
  else()
    set(PER_STATIC ON)
  endif()
endif()

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
  elseif(MSVC_VERSION EQUAL 1900)
    set(PER_RUNTIME vc14)
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
