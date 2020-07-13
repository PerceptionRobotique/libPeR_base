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
# This file generates the PERConfig.cmake file: 
#  Part 1/3: ${BIN_DIR}/PERConfig.cmake              -> For use *without* "make install"
#  Part 2/3: ${BIN_DIR}/unix-install/PERConfig.cmake -> For use with "make install"
#  Part 3/3: ${BIN_DIR}/win-install/PERConfig.cmake  -> For use within binary installers/packages
# Freely inspired from the ViSP library. 
#
# Authors:
# Guillaume Caron
#
#############################################################################

# Macro that returns the relative path to go from a child folder to the parent folder
# input: path_to_child
# output: path_to_parent, the relative path to go from path_to_child to parent
# example: if input =lib/x86_64-linux-gnu, then output=../..
macro(get_path_to_parent path_to_child path_to_parent)
  set(${path_to_parent} "")
  set(input_ "${path_to_child}")
  while(input_)
    if(input_)
      set(${path_to_parent} "${${path_to_parent}}../")
    endif()
    get_filename_component(input_ "${input_}" PATH)
  endwhile(input_)
endmacro()

# Here we determine the relative path from ./${CMAKE_INSTALL_LIBDIR} to its parent folder
# if CMAKE_INSTALL_LIBDIR=lib, then PER_INSTALL_LIBDIR_TO_PARENT=../
# if CMAKE_INSTALL_LIBDIR=lib/x86_64-linux-gnu, then PER_INSTALL_LIBDIR_TO_PARENT=../..
get_path_to_parent(${CMAKE_INSTALL_LIBDIR} PER_INSTALL_LIBDIR_TO_PARENT)

#build list of modules available for the PER user
set(PER_LIB_COMPONENTS "")
foreach(m ${PER_MODULES_PUBLIC})
  list(INSERT PER_LIB_COMPONENTS 0 ${${m}_MODULE_DEPS_OPT} ${m})
endforeach()
pr_list_unique(PER_LIB_COMPONENTS)
#message("PER_LIB_COMPONENTS: ${PER_LIB_COMPONENTS}")
set(PER_MODULES_CONFIGCMAKE ${PER_LIB_COMPONENTS})
pr_list_filterout(PER_LIB_COMPONENTS "^per_")
if(PER_LIB_COMPONENTS)
  list(REMOVE_ITEM PER_MODULES_CONFIGCMAKE ${PER_LIB_COMPONENTS})
endif()


# -------------------------------------------------------------------------------------------
#  Part 1/3: ${BIN_DIR}/PERConfig.cmake              -> For use *without* "make install"
# -------------------------------------------------------------------------------------------

# Export the library
export(TARGETS ${PERModules_TARGETS} FILE "${PROJECT_BINARY_DIR}/PERModules.cmake")

## Update include dirs
set(PER_INCLUDE_DIRS_CONFIGCMAKE "${PER_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}")
foreach(m ${PER_MODULES_BUILD})
  if(EXISTS "${PER_MODULE_${m}_LOCATION}/include")
    list(APPEND PER_INCLUDE_DIRS_CONFIGCMAKE "${PER_MODULE_${m}_LOCATION}/include")
  endif()
  list(APPEND PER_INCLUDE_DIRS_CONFIGCMAKE ${PER_MODULE_${m}_INC_DEPS})
endforeach()
pr_list_unique(PER_INCLUDE_DIRS_CONFIGCMAKE)

configure_file(
  cmake/templates/PERConfig.cmake.in
  ${PER_BINARY_DIR}/PERConfig.cmake
  IMMEDIATE @ONLY
)

configure_file(
  cmake/templates/PERConfigVersion.cmake.in
  ${PER_BINARY_DIR}/PERConfigVersion.cmake
  IMMEDIATE @ONLY
)

configure_file(
  cmake/PERUse.cmake.in
  ${PER_BINARY_DIR}/PERUse.cmake
  IMMEDIATE @ONLY
)

# --------------------------------------------------------------------------------------------
#  Part 2/3: ${BIN_DIR}/unix-install/PERConfig.cmake -> For use *with* "make install"
# -------------------------------------------------------------------------------------------

if(UNIX)
  set(PER_INCLUDE_DIRS_CONFIGCMAKE "\${PER_INSTALL_PATH}/${CMAKE_INSTALL_INCLUDEDIR}")
  foreach(m ${PER_MODULES_BUILD})
    list(APPEND PER_INCLUDE_DIRS_CONFIGCMAKE ${PER_MODULE_${m}_INC_DEPS})
  endforeach()
  pr_list_unique(PER_INCLUDE_DIRS_CONFIGCMAKE)

  configure_file(
    cmake/templates/PERConfig.cmake.in
    ${PER_BINARY_DIR}/unix-install/PERConfig.cmake
    IMMEDIATE @ONLY
  )

  configure_file(
    cmake/templates/PERConfigVersion.cmake.in
    ${PER_BINARY_DIR}/unix-install/PERConfigVersion.cmake
    IMMEDIATE @ONLY
  )

  configure_file(
    cmake/PERUse.cmake.in
    ${PER_BINARY_DIR}/unix-install/PERUse.cmake
    IMMEDIATE @ONLY
  )

  install(FILES
    ${PER_BINARY_DIR}/unix-install/PERConfig.cmake
    ${PER_BINARY_DIR}/unix-install/PERConfigVersion.cmake
    ${PER_BINARY_DIR}/unix-install/PERUse.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/per"
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_WRITE
    COMPONENT dev
  )

  # Install the export set for use with the install-tree
  install(EXPORT PERModules
    FILE PERModules.cmake
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/per"
    COMPONENT dev
  )
endif()

# --------------------------------------------------------------------------------------------
#  Part 3/3: ${BIN_DIR}/win-install/PERConfig.cmake  -> For use within binary installers/packages
# --------------------------------------------------------------------------------------------
if(WIN32)
  set(PER_INCLUDE_DIRS_CONFIGCMAKE "\${PER_CONFIG_PATH}/${CMAKE_INSTALL_INCLUDEDIR}")
  foreach(m ${PER_MODULES_BUILD})
    list(APPEND PER_INCLUDE_DIRS_CONFIGCMAKE ${PER_MODULE_${m}_INC_DEPS})
  endforeach()
  pr_list_unique(PER_INCLUDE_DIRS_CONFIGCMAKE)

  configure_file(
    cmake/templates/PERConfig.cmake.in
    ${PER_BINARY_DIR}/win-install/PERConfig.cmake
    IMMEDIATE @ONLY
  )

  configure_file(
    cmake/templates/PERConfigVersion.cmake.in
    ${PER_BINARY_DIR}/win-install/PERConfigVersion.cmake
    IMMEDIATE @ONLY
  )

  configure_file(
    cmake/PERUse.cmake.in
    ${PER_BINARY_DIR}/win-install/PERUse.cmake
    IMMEDIATE @ONLY
  )

  if(BUILD_SHARED_LIBS)
    install(FILES
      "${CMAKE_BINARY_DIR}/win-install/PERConfig.cmake"
      "${CMAKE_BINARY_DIR}/win-install/PERUse.cmake"
      DESTINATION "${PER_INSTALL_BINARIES_PREFIX}${CMAKE_INSTALL_LIBDIR}"
      COMPONENT dev)
    install(EXPORT PERModules 
      DESTINATION "${PER_INSTALL_BINARIES_PREFIX}${CMAKE_INSTALL_LIBDIR}"
      FILE PERModules.cmake 
      COMPONENT dev)
  else()
    install(FILES
      "${CMAKE_BINARY_DIR}/win-install/PERConfig.cmake"
      "${CMAKE_BINARY_DIR}/win-install/PERUse.cmake"
      DESTINATION "${PER_INSTALL_BINARIES_PREFIX}static${CMAKE_INSTALL_LIBDIR}"
      COMPONENT dev)
    install(EXPORT PERModules 
      DESTINATION "${PER_INSTALL_BINARIES_PREFIX}static${CMAKE_INSTALL_LIBDIR}"
      FILE PERModules.cmake 
      COMPONENT dev)
  endif()

  install(FILES
    "cmake/PERConfig.cmake"
    "${PER_BINARY_DIR}/win-install/PERConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_PREFIX}"
    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_WRITE
    COMPONENT dev
  )
endif()
