#############################################################################
#
# This file is part of the libPR software.
# Copyright (C) 2022 by MIS lab (UPJV) and JRL (CNRS-AIST). All rights reserved.
#
# See http://mis.u-picardie.fr/~g-caron/fr/index.php?page=7 for more information.
#
# This software was developed at:
# CNRS-AIST JRL
# 1-1-1 Umezono
# Tsukuba
# Japan
#
# This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
# WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#
# Description:
# configuration file for libPeR example. Freely inspired from the CMakeList tree of the 
# ViSP library. 
#
# Authors:
# Guillaume Caron
#
#############################################################################

project(example-FisheyeEquidistant2OmniCameraModelConversion)

cmake_minimum_required(VERSION 2.6)

find_package(PER REQUIRED per_core per_io)

set(example_cpp
  FisheyeEquidistant2OmniCameraModelConversion.cpp
)

foreach(cpp ${example_cpp})
  per_add_target(${cpp})
  if(COMMAND pr_add_dependency)
    pr_add_dependency(${cpp} "examples")
  endif()
endforeach()

