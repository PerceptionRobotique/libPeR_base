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
# Freely inspired from the CMakeList tree of the 
# ViSP library.
#
# Authors:
# Guillaume Caron
#
#############################################################################

# platform-specific config file
#configure_file("${PER_SOURCE_DIR}/cmake/templates/vpConfig.h.in" "${PER_INCLUDE_DIR}/per/core/vpConfig.h")
#install(FILES "${PER_INCLUDE_DIR}/per/core/vpConfig.h"
#  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/per/core
#  COMPONENT dev
#)

#----------------------------------------------------------------------
# information file
#----------------------------------------------------------------------
configure_file(${PER_SOURCE_DIR}/cmake/templates/PER-third-party.txt.in "${PER_BINARY_DIR}/PER-third-party.txt")

# ----------------------------------------------------------------------------
#  per_modules.h based on actual modules list
# ----------------------------------------------------------------------------
set(PER_MODULE_DEFINITIONS_CONFIGMAKE "#ifndef __per_modules_h__\n#define __per_modules_h__\n\n")

set(PER_MOD_LIST ${PER_MODULES_PUBLIC})
pr_list_sort(PER_MOD_LIST)
foreach(m ${PER_MOD_LIST})
  if(m MATCHES "^per_")
    string(REGEX REPLACE "^per_" "" m "${m}")
  endif()
  string(TOUPPER "${m}" m)
  set(PER_MODULE_DEFINITIONS_CONFIGMAKE "${PER_MODULE_DEFINITIONS_CONFIGMAKE}#define PER_HAVE_MODULE_${m}\n")
endforeach()

set(PER_MODULE_DEFINITIONS_CONFIGMAKE "${PER_MODULE_DEFINITIONS_CONFIGMAKE}\n#endif\n")

configure_file("${PER_SOURCE_DIR}/cmake/templates/per_modules.h.in" "${PER_INCLUDE_DIR}/per/per_modules.h")
install(FILES "${PER_INCLUDE_DIR}/per/per_modules.h"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/per
  COMPONENT dev
)

# ----------------------------------------------------------------------------
#  install old headers
# ----------------------------------------------------------------------------
file(GLOB old_hdrs "${PER_INCLUDE_DIR}/per/*.h")
install(FILES ${old_hdrs}
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/per
  COMPONENT dev
)
