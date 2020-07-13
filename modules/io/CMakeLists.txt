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

# Add optional 3rd parties
set(opt_incs "")
set(opt_libs "")

if(USE_XML2)
  list(APPEND opt_incs ${XML2_INCLUDE_DIRS})
  list(APPEND opt_libs ${XML2_LIBRARIES})
endif()

pr_add_module(io per_core)
pr_glob_module_sources()
pr_module_include_directories(${opt_incs})
pr_create_module(${opt_libs})
#pr_add_tests(DEPENDS_ON pr_features)
