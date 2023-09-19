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
# ViSP library. 
#
# Authors:
# Guillaume Caron
#
#############################################################################

project(PER-examples)

cmake_minimum_required(VERSION 2.6)

find_package(PER REQUIRED)

if(MSVC)
  if(NOT PER_SHARED)
    foreach(flag_var
            CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
            CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
            CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
            CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      if(${flag_var} MATCHES "/MD")
        string(REGEX REPLACE "/MD" "/MT" ${flag_var} "${${flag_var}}")
      endif()
      if(${flag_var} MATCHES "/MDd")
        string(REGEX REPLACE "/MDd" "/MTd" ${flag_var} "${${flag_var}}")
      endif()
    endforeach(flag_var)

    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /NODEFAULTLIB:atlthunk.lib /NODEFAULTLIB:msvcrt.lib /NODEFAULTLIB:msvcrtd.lib")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} /NODEFAULTLIB:libcmt.lib")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /NODEFAULTLIB:libcmtd.lib")
  endif()
endif()

per_add_subdirectory(pointPerspectiveProjection REQUIRED_DEPS per_core)
per_add_subdirectory(saveAndLoadStereocameraModelXML REQUIRED_DEPS per_core per_io)
per_add_subdirectory(saveAndLoadColoredPointCloud REQUIRED_DEPS per_core per_core_extended per_io per_io_extended)
per_add_subdirectory(FisheyeEquidistant2OmniCameraModelConversion REQUIRED_DEPS per_core per_io)
per_add_subdirectory(Omni2PolyCartCameraModelConversion REQUIRED_DEPS per_core per_io)
per_add_subdirectory(RationalPolynomial2PlumbBoCameraModelConversion REQUIRED_DEPS per_core per_io)
