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
# libPR configuration file. Freely inspired from the 
# ViSP library. 
#
# Authors:
# Guillaume Caron
#
#############################################################################

if(NOT PER_MODULES_PATH)
  set(PER_MODULES_PATH "${CMAKE_CURRENT_SOURCE_DIR}")
endif()

# extra modules are outside PER. Useful for contrib or specific additional modules
pr_glob_modules(${PER_MODULES_PATH} EXTRA ${PER_CONTRIB_MODULES_PATH})
