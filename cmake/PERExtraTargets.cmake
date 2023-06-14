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
#   Uninstall target, for "make uninstall"
# ----------------------------------------------------------------------------
configure_file(
  cmake/templates/cmake_uninstall.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

add_custom_target(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake")

if(ENABLE_SOLUTION_FOLDERS)
  set_target_properties(uninstall PROPERTIES FOLDER "CMakeTargets")
endif()

# ----------------------------------------------------------------------------
#   Doxygen documentation target, for "make per_doc" and "make html-doc" (to keep compat with previous versions)
# ----------------------------------------------------------------------------
if(DOXYGEN_FOUND)
  add_custom_target(html-doc ${DOXYGEN_EXECUTABLE} ${PER_DOC_DIR}/config-doxygen) # for compat with previous versions
  add_custom_target(per_doc ${DOXYGEN_EXECUTABLE} ${PER_DOC_DIR}/config-doxygen)
  if(ENABLE_SOLUTION_FOLDERS)
    set_target_properties(per_doc PROPERTIES FOLDER "extra")
    set_target_properties(html-doc PROPERTIES FOLDER "extra")
  endif()
endif()

# ----------------------------------------------------------------------------
#   Tests target, for make per_tests
# ----------------------------------------------------------------------------
if(BUILD_TESTS)
  add_custom_target(per_tests)
  if(ENABLE_SOLUTION_FOLDERS)
    set_target_properties(per_tests PROPERTIES FOLDER "extra")
  endif()
endif()

# ----------------------------------------------------------------------------
#   Tests target, for make per_examples
# ----------------------------------------------------------------------------
if(BUILD_EXAMPLES)
  add_custom_target(per_examples)
  if(ENABLE_SOLUTION_FOLDERS)
    set_target_properties(per_examples PROPERTIES FOLDER "extra")
  endif()
endif()

# ----------------------------------------------------------------------------
#   Target building all libPR modules
# ----------------------------------------------------------------------------
add_custom_target(per_modules)
if(ENABLE_SOLUTION_FOLDERS)
  set_target_properties(per_modules PROPERTIES FOLDER "extra")
endif()

