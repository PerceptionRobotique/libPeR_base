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

# adds include directories in such way that directories from the PER source tree go first
function(pr_include_directories)
  pr_debug_message("pr_include_directories( ${ARGN} )")
  set(__add_before "")
  foreach(dir ${ARGN})
    get_filename_component(__abs_dir "${dir}" ABSOLUTE)
    string(REPLACE "+" "\\+" __PER_BINARY_DIR_filtered ${PER_BINARY_DIR})
#    if("${__abs_dir}" MATCHES "^${PER_SOURCE_DIR}" OR "${__abs_dir}" MATCHES "^${PER_BINARY_DIR}")
    if("${__abs_dir}" MATCHES "^${PER_SOURCE_DIR}" OR "${__abs_dir}" MATCHES "^${__PER_BINARY_DIR_filtered}")
      list(APPEND __add_before "${dir}")
    else()
      include_directories(AFTER SYSTEM "${dir}")
    endif()
  endforeach()
  include_directories(BEFORE ${__add_before})
endfunction()

# adds include directories in such way that directories from the PER source tree go first
function(pr_target_include_directories target)
  set(__params "")
  foreach(dir ${ARGN})
    get_filename_component(__abs_dir "${dir}" ABSOLUTE)
    string(REPLACE "+" "\\+" __PER_BINARY_DIR_filtered ${PER_BINARY_DIR})
#    if("${__abs_dir}" MATCHES "^${PER_SOURCE_DIR}" OR "${__abs_dir}" MATCHES "^${PER_BINARY_DIR}")
    if("${__abs_dir}" MATCHES "^${PER_SOURCE_DIR}" OR "${__abs_dir}" MATCHES "^${__PER_BINARY_DIR_filtered}")
      list(APPEND __params "${__abs_dir}")
    else()
      list(APPEND __params "${dir}")
    endif()
  endforeach()
  if(__params)
    if(CMAKE_VERSION VERSION_LESS 2.8.11)
      include_directories(${__params})
    else()
      if(TARGET ${target})
        target_include_directories(${target} PRIVATE ${__params})
      else()
        set(__new_inc "${PR_TARGET_INCLUDE_DIRS_${target}};${__params}")
        set(PR_TARGET_INCLUDE_DIRS_${target} "${__new_inc}" CACHE INTERNAL "")
      endif()
    endif()
  endif()
endfunction()

# clears all passed variables
macro(pr_clear_vars)
  foreach(_var ${ARGN})
    unset(${_var} CACHE)
  endforeach()
endmacro()

# print message
macro(pr_debug_message)
  #string(REPLACE ";" " " __msg "${ARGN}")
  #message(STATUS "${__msg}")
endmacro()

# remove all matching elements from the list
macro(pr_list_filterout lst regex)
  foreach(item ${${lst}})
    if(item MATCHES "${regex}")
      list(REMOVE_ITEM ${lst} "${item}")
    endif()
  endforeach()
endmacro()

# stable & safe duplicates removal macro
macro(pr_list_unique __lst)
  if(${__lst})
    list(REMOVE_DUPLICATES ${__lst})
  endif()
endmacro()

# safe list sorting macro
macro(pr_list_sort __lst)
  if(${__lst})
    list(SORT ${__lst})
  endif()
endmacro()

# gets and removes the first element from list
macro(pr_list_pop_front LST VAR)
  if(${LST})
    list(GET ${LST} 0 ${VAR})
    list(REMOVE_AT ${LST} 0)
  else()
    set(${VAR} "")
  endif()
endmacro()

# remove cmake ; list separator
macro(pr_list_remove_separator __lst)
  if(${__lst})
    set(__lst_reformated "")
    foreach(element ${${__lst}})
      set(__lst_reformated "${__lst_reformated} ${element}")
    endforeach()
    set(${__lst} ${__lst_reformated})
  endif()
endmacro()

# Provides an option that the user can optionally select.
# Can accept condition to control when option is available for user.
# Usage:
#   pr_option(<option_variable>
#             <package to find to activate the option>
#             <QUIET or "">
#             "help string describing the option"
#             "advanced list of vars separated by ;"
#             <initial value or boolean expression> [IF <condition>])
macro(PR_OPTION variable package quiet description advanced value)
  set(__option TRUE)
  set(__value ${value})
  set(__condition "")
  set(__varname "__value")
  set(__p ${package})

  # get the first package considered as the main package from a list: ie "Zlib;MyZlib"
  set(__first_package "")
  foreach(p ${package})
    if(${p} MATCHES "^My")
      string(REGEX REPLACE "^My" "" p "${p}")
    endif()
    set(__first_package ${p})
    break()
  endforeach()

  if(NOT ${__first_package} STREQUAL "")
    # get the first package name from the list
    #list(GET ${package} 0 ${__package})
    string(TOLOWER "${__first_package}" __package_lower)
    string(TOUPPER "${__first_package}" __package_upper) # useful for Qt -> QT_FOUND
  endif()

  # set PER_HAVE_<package>_FOUND="no"
  set(__alias_have_found_str PER_HAVE_${__package_upper}_FOUND)
  set(${__alias_have_found_str} "no")

  foreach(arg ${ARGN})
    if(arg STREQUAL "IF" OR arg STREQUAL "if")
      set(__varname "__condition")
    else()
      list(APPEND ${__varname} ${arg})
    endif()
  endforeach()
  unset(__varname)
  if(__condition STREQUAL "")
    set(__condition 2 GREATER 1)
  endif()
  if(${__condition})

    if(NOT ${__first_package} STREQUAL "")
      foreach(p ${package})
        if("${quiet}" STREQUAL "")
          find_package(${p})
        else()
          find_package(${p} ${quiet})
        endif()
        if(${__package_upper}_FOUND OR ${__first_package}_FOUND)
          set(__option TRUE)
          break()
        else()
          set(__option FALSE)
        endif()
      endforeach()
    endif()
    if(${__option})
      if(__value MATCHES ";")
        if(${__value})
          option(${variable} "${description}" ON)
        else()
          option(${variable} "${description}" OFF)
        endif()
      elseif(DEFINED ${__value})
        if(${__value})
          option(${variable} "${description}" ON)
        else()
          option(${variable} "${description}" OFF)
        endif()
      else()
        option(${variable} "${description}" ${__value})
      endif()
    else()
      unset(${variable} CACHE)
    endif()
    unset(__condition)
    unset(__value)

  else()
    set(${variable} OFF)
  endif()
  foreach(a ${advanced})
    mark_as_advanced(${a})
  endforeach()
  if(${variable} AND NOT ${__first_package} STREQUAL "")
    # set PER_HAVE_<package>=TRUE and PER_HAVE_<package>_FOUND="yes"
    message(STATUS "${__package_lower} found")
    set(${__alias_have_found_str} "yes") # for PER-third-party.txt
  endif()
  unset(__option)
  unset(__alias_have)
  unset(__alias_have_found_str)
endmacro()

# Provides a macro to set a var.
# Can accept condition to set var.
# Usage:
#   pr_set(<option_variable>
#          <initial value or boolean expression> [IF <condition>])
macro(PR_SET variable value)
  set(__value ${value})
  set(__condition "")
  set(__varname "__value")

  foreach(arg ${ARGN})
    if(arg STREQUAL "IF" OR arg STREQUAL "if")
      set(__varname "__condition")
    else()
      list(APPEND ${__varname} ${arg})
    endif()
  endforeach()
  unset(__varname)
  if(__condition STREQUAL "")
    set(__condition 2 GREATER 1)
  endif()
  if(${__condition})
    if(__value MATCHES ";")
      if(${__value})
        set(${variable} TRUE)
      else()
        set(${variable} FALSE)
      endif()
    elseif(DEFINED ${__value})
      if(${__value})
        set(${variable} TRUE)
      else()
        set(${variable} FALSE)
      endif()
    else()
      set(${variable} ${__value})
    endif()
  endif()
  unset(__condition)
  unset(__value)
endmacro()

# short command to setup source group
function(pr_source_group group)
  cmake_parse_arguments(SG "" "DIRBASE" "GLOB;GLOB_RECURSE;FILES" ${ARGN})
  set(files "")
  if(SG_FILES)
    list(APPEND files ${SG_FILES})
  endif()
  if(SG_GLOB)
    file(GLOB srcs ${SG_GLOB})
    list(APPEND files ${srcs})
  endif()
  if(SG_GLOB_RECURSE)
    file(GLOB_RECURSE srcs ${SG_GLOB_RECURSE})
    list(APPEND files ${srcs})
  endif()
  if(SG_DIRBASE)
    foreach(f ${files})
      file(RELATIVE_PATH fpart "${SG_DIRBASE}" "${f}")
      if(fpart MATCHES "^\\.\\.")
        message(AUTHOR_WARNING "Can't detect subpath for source_group command: Group=${group} FILE=${f} DIRBASE=${SG_DIRBASE}")
        set(fpart "")
      else()
        get_filename_component(fpart "${fpart}" PATH)
        if(fpart)
          set(fpart "/${fpart}") # add '/'
          string(REPLACE "/" "\\" fpart "${fpart}")
        endif()
      endif()
      source_group("${group}${fpart}" FILES ${f})
    endforeach()
  else()
    source_group(${group} FILES ${files})
  endif()
endfunction()

# convert list of paths to full paths
macro(pr_convert_to_full_paths VAR)
  if(${VAR})
    set(__tmp "")
    foreach(path ${${VAR}})
      get_filename_component(${VAR} "${path}" ABSOLUTE)
      list(APPEND __tmp "${${VAR}}")
    endforeach()
    set(${VAR} ${__tmp})
    unset(__tmp)
  endif()
endmacro()

# add install command
function(pr_install_target)
  install(TARGETS ${ARGN})

  set(isPackage 0)
  unset(__package)
  unset(__target)
  foreach(e ${ARGN})
    if(NOT DEFINED __target)
      set(__target "${e}")
    endif()
    if(isPackage EQUAL 1)
      set(__package "${e}")
      break()
    endif()
    if(e STREQUAL "EXPORT")
      set(isPackage 1)
    endif()
  endforeach()

  if(DEFINED __package)
    list(APPEND ${__package}_TARGETS ${__target})
    set(${__package}_TARGETS "${${__package}_TARGETS}" CACHE INTERNAL "List of ${__package} targets")
  endif()

  if(INSTALL_CREATE_DISTRIB)
    if(MSVC AND NOT BUILD_SHARED_LIBS)
      set(__target "${ARGV0}")

      set(isArchive 0)
      set(isDst 0)
      unset(__dst)
      foreach(e ${ARGN})
        if(isDst EQUAL 1)
          set(__dst "${e}")
          break()
        endif()
        if(isArchive EQUAL 1 AND e STREQUAL "DESTINATION")
          set(isDst 1)
        endif()
        if(e STREQUAL "ARCHIVE")
          set(isArchive 1)
        else()
          set(isArchive 0)
        endif()
      endforeach()

#      message(STATUS "Process ${__target} dst=${__dst}...")
      if(DEFINED __dst)
        if(CMAKE_VERSION VERSION_LESS 2.8.12)
          get_target_property(fname ${__target} LOCATION_DEBUG)
          if(fname MATCHES "\\.lib$")
            string(REGEX REPLACE "\\.lib$" ".pdb" fname "${fname}")
            install(FILES ${fname} DESTINATION ${__dst} CONFIGURATIONS Debug)
          endif()

          get_target_property(fname ${__target} LOCATION_RELEASE)
          if(fname MATCHES "\\.lib$")
            string(REGEX REPLACE "\\.lib$" ".pdb" fname "${fname}")
            install(FILES ${fname} DESTINATION ${__dst} CONFIGURATIONS Release)
          endif()
        else()
          # CMake 2.8.12 brokes PDB support in STATIC libraries for MSVS
        endif()
      endif()
    endif()
  endif()
endfunction()

function(pr_target_link_libraries target)
  set(LINK_DEPS ${ARGN})
  target_link_libraries(${target} ${LINK_DEPS})
endfunction()

function(_pr_append_target_includes target)
  if(DEFINED PR_TARGET_INCLUDE_DIRS_${target})
    target_include_directories(${target} PRIVATE ${PR_TARGET_INCLUDE_DIRS_${target}})
    unset(PR_TARGET_INCLUDE_DIRS_${target} CACHE)
  endif()
endfunction()

function(pr_add_executable target)
  add_executable(${target} ${ARGN})
  _pr_append_target_includes(${target})
endfunction()

function(pr_add_library target)
  add_library(${target} ${ARGN})

  _pr_append_target_includes(${target})
endfunction()

# Macros that checks if package have been installed.
# After it set vars:
#   <package>_FOUND
#   <package>_INCLUDE_DIRS
#   <package>_LIBRARIES
#   PER_HAVE_<package>
macro(PR_CHECK_PACKAGE package)
  set(ALIAS                 ${package})
  string(TOUPPER "${ALIAS}" ALIAS_UPPER) # useful for OpenGL
  set(ALIAS_FOUND           ${ALIAS}_FOUND)
  set(ALIAS_UPPER_FOUND     ${ALIAS_UPPER}_FOUND)
  set(ALIAS_INCLUDE_DIRS    ${ALIAS}_INCLUDE_DIRS)
  set(ALIAS_LIBRARIES       ${ALIAS}_LIBRARIES)
  set(ALIAS_PER_HAVE       PER_HAVE_${ALIAS})
  set(ALIAS_UPPER_PER_HAVE PER_HAVE_${ALIAS_UPPER})

  find_package(${ALIAS})

  if(${ALIAS_FOUND} OR ${ALIAS_UPPER_FOUND})
    set(${ALIAS_PER_HAVE} 1)
    set(${ALIAS_UPPER_PER_HAVE} 1)
  endif()
endmacro()

# Macro the get the list of subdirs from the path
# var: returned variable name
# path: path from witch relative subdirs are
macro(pr_get_relative_subdirs var path)
  set(ALIAS                 ${var})
  file(GLOB_RECURSE rel_path_lst_ RELATIVE ${path} ${path}/*)
    set(${ALIAS} "")
    foreach(f ${rel_path_lst_})
      get_filename_component(d ${f} PATH)
      list(APPEND ${ALIAS} ${d})
    endforeach()
    list(REMOVE_DUPLICATES ${ALIAS})
endmacro()
