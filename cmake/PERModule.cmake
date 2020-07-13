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
#
# Authors:
# Guillaume Caron
#
#############################################################################

# Local variables (set for each module):
#
# name       - short name in lower case i.e. core
# the_module - full name in lower case i.e. per_core

# Global variables:
#
# PER_MODULE_${the_module}_LOCATION
# PER_MODULE_${the_module}_BINARY_DIR
# PER_MODULE_${the_module}_DESCRIPTION
# PER_MODULE_${the_module}_CLASS - PUBLIC|INTERNAL|BINDINGS
# PER_MODULE_${the_module}_HEADERS
# PER_MODULE_${the_module}_SOURCES
# PER_MODULE_${the_module}_DEPS - final flattened set of module dependencies
# PER_MODULE_${the_module}_DEPS_TO_LINK - differs from above for world build only
# PER_MODULE_${the_module}_DEPS_EXT - non-module dependencies
# PER_MODULE_${the_module}_REQ_DEPS
# PER_MODULE_${the_module}_OPT_DEPS
# PER_MODULE_${the_module}_PRIVATE_REQ_DEPS
# PER_MODULE_${the_module}_PRIVATE_OPT_DEPS
# PER_MODULE_${the_module}_CHILDREN - list of submodules for compound modules (cmake >= 2.8.8)
# HAVE_${the_module} - for fast check of module availability

# To control the setup of the module you could also set:
# the_description - text to be used as current module description
# PER_MODULE_TYPE - STATIC|SHARED - set to force override global settings for current module
# BUILD_${the_module}_INIT - ON|OFF (default ON) - initial value for BUILD_${the_module}
# PER_MODULE_CHILDREN - list of submodules

# The verbose template for libPR module:
#
#   pr_add_module(modname <dependencies>)
#   pr_glob_module_sources((<extra sources&headers>)
#                          or glob them manually and pr_set_module_sources(...)
#   pr_module_include_directories(<extra include directories>)
#   pr_create_module()
#
# If module have no "extra" then you can define it in one line:
#
#   pr_define_module(modname <dependencies>)

# clean flags for modules enabled on previous cmake run
# this is necessary to correctly handle modules removal
foreach(mod ${PER_MODULES_BUILD} ${PER_MODULES_DISABLED_USER} ${PER_MODULES_DISABLED_AUTO} ${PER_MODULES_DISABLED_FORCE})
  if(HAVE_${mod})
    unset(HAVE_${mod} CACHE)
  endif()
  unset(PER_MODULE_${mod}_REQ_DEPS CACHE)
  unset(PER_MODULE_${mod}_OPT_DEPS CACHE)
  unset(PER_MODULE_${mod}_PRIVATE_REQ_DEPS CACHE)
  unset(PER_MODULE_${mod}_PRIVATE_OPT_DEPS CACHE)
  unset(PER_MODULE_${mod}_LINK_DEPS CACHE)
  unset(PER_MODULE_${mod}_INC_DEPS CACHE)
endforeach()

# clean modules info which needs to be recalculated
set(PER_MODULES_PUBLIC         "" CACHE INTERNAL "List of PER modules marked for export")
set(PER_MODULES_BUILD          "" CACHE INTERNAL "List of PER modules included into the build")
set(PER_MODULES_DISABLED_USER  "" CACHE INTERNAL "List of PER modules explicitly disabled by user")
set(PER_MODULES_DISABLED_AUTO  "" CACHE INTERNAL "List of PER modules implicitly disabled due to dependencies")
set(PER_MODULES_DISABLED_FORCE "" CACHE INTERNAL "List of PER modules which can not be build in current configuration")

# adds dependencies to PER module
# Usage:
#   add_dependencies(per_<name> [REQUIRED] [<list of dependencies>] [OPTIONAL <list of modules>])
# Notes:
# * <list of dependencies> - can include full names of modules or full pathes to shared/static libraries or cmake targets
macro(pr_add_dependencies full_modname)
  pr_debug_message("pr_add_dependencies(" ${full_modname} ${ARGN} ")")
  #we don't clean the dependencies here to allow this macro several times for every module
  foreach(d "REQUIRED" ${ARGN})
    if(d STREQUAL "REQUIRED")
      set(__depsvar PER_MODULE_${full_modname}_REQ_DEPS)
    elseif(d STREQUAL "OPTIONAL")
      set(__depsvar PER_MODULE_${full_modname}_OPT_DEPS)
    elseif(d STREQUAL "PRIVATE_REQUIRED")
      set(__depsvar PER_MODULE_${full_modname}_PRIVATE_REQ_DEPS)
    elseif(d STREQUAL "PRIVATE_OPTIONAL")
      set(__depsvar PER_MODULE_${full_modname}_PRIVATE_OPT_DEPS)
    else()
      list(APPEND ${__depsvar} "${d}")
    endif()
  endforeach()
  unset(__depsvar)

  pr_list_unique(PER_MODULE_${full_modname}_REQ_DEPS)
  pr_list_unique(PER_MODULE_${full_modname}_OPT_DEPS)
  pr_list_unique(PER_MODULE_${full_modname}_PRIVATE_REQ_DEPS)
  pr_list_unique(PER_MODULE_${full_modname}_PRIVATE_OPT_DEPS)

  set(PER_MODULE_${full_modname}_REQ_DEPS ${PER_MODULE_${full_modname}_REQ_DEPS}
    CACHE INTERNAL "Required dependencies of ${full_modname} module")
  set(PER_MODULE_${full_modname}_OPT_DEPS ${PER_MODULE_${full_modname}_OPT_DEPS}
    CACHE INTERNAL "Optional dependencies of ${full_modname} module")
  set(PER_MODULE_${full_modname}_PRIVATE_REQ_DEPS ${PER_MODULE_${full_modname}_PRIVATE_REQ_DEPS}
    CACHE INTERNAL "Required private dependencies of ${full_modname} module")
  set(PER_MODULE_${full_modname}_PRIVATE_OPT_DEPS ${PER_MODULE_${full_modname}_PRIVATE_OPT_DEPS}
    CACHE INTERNAL "Optional private dependencies of ${full_modname} module")
endmacro()

# declare new PER module in current folder
# Usage:
#   pr_add_module(<name> [INTERNAL|BINDINGS] [REQUIRED] [<list of dependencies>] [OPTIONAL <list of optional dependencies>])
# Example:
#   pr_add_module(mymodule INTERNAL per_core OPTIONAL per_ar)
macro(pr_add_module _name)
  pr_debug_message("pr_add_module(" ${_name} ${ARGN} ")")
  string(TOLOWER "${_name}" name)
  set(the_module per_${name})
  #message("Found module: ${the_module}")
  
  # the first pass - collect modules info, the second pass - create targets
  if(PER_INITIAL_PASS)
    #guard agains redefinition
    if(";${PER_MODULES_BUILD};${PER_MODULES_DISABLED_USER};" MATCHES ";${the_module};")
      message(FATAL_ERROR "Redefinition of the ${the_module} module.
  at:                    ${CMAKE_CURRENT_SOURCE_DIR}
  previously defined at: ${PER_MODULE_${the_module}_LOCATION}
")
    endif()

    if(NOT DEFINED the_description)
      set(the_description "The PER ${name} module")
    endif()

    if(NOT DEFINED BUILD_${the_module}_INIT)
      set(BUILD_${the_module}_INIT ON)
    endif()

    # create option to enable/disable this module
    option(BUILD_MODULE_${the_module} "Include ${the_module} module into PER build" ${BUILD_${the_module}_INIT})
    
    # remember the module details
    set(PER_MODULE_${the_module}_DESCRIPTION "${the_description}" CACHE INTERNAL "Brief description of ${the_module} module")
    set(PER_MODULE_${the_module}_LOCATION    "${CMAKE_CURRENT_SOURCE_DIR}" CACHE INTERNAL "Location of ${the_module} module sources")

    set(PER_MODULE_${the_module}_LINK_DEPS "" CACHE INTERNAL "")
    set(PER_MODULE_${the_module}_INC_DEPS "" CACHE INTERNAL "")

    # parse list of dependencies
    if("${ARGV1}" STREQUAL "INTERNAL" OR "${ARGV1}" STREQUAL "BINDINGS")
      set(PER_MODULE_${the_module}_CLASS "${ARGV1}" CACHE INTERNAL "The category of the module")
      set(__pr_argn__ ${ARGN})
      list(REMOVE_AT __pr_argn__ 0)
      pr_add_dependencies(${the_module} ${__pr_argn__})
      unset(__pr_argn__)
    else()
      set(PER_MODULE_${the_module}_CLASS "PUBLIC" CACHE INTERNAL "The category of the module")
      pr_add_dependencies(${the_module} ${ARGN})
      if(BUILD_MODULE_${the_module})
        set(PER_MODULES_PUBLIC ${PER_MODULES_PUBLIC} "${the_module}" CACHE INTERNAL "List of PER modules marked for export")
      endif()
    endif()

    if(BUILD_MODULE_${the_module})
      set(PER_MODULES_BUILD ${PER_MODULES_BUILD} "${the_module}" CACHE INTERNAL "List of PER modules included into the build")
    else()
      set(PER_MODULES_DISABLED_USER ${PER_MODULES_DISABLED_USER} "${the_module}" CACHE INTERNAL "List of PER modules explicitly disabled by user")
    endif()

    # add submodules if any
    set(PER_MODULE_${the_module}_CHILDREN "${PER_MODULE_CHILDREN}" CACHE INTERNAL "List of ${the_module} submodules")

    # stop processing of current file
    return()
  else()
    set(PER_MODULE_${the_module}_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE INTERNAL "")
    if(NOT BUILD_MODULE_${the_module})
      return() # extra protection from redefinition
    endif()
  endif()
endmacro()

# remove per_ prefix from name
macro(pr_short_module_name name)
  if(${name} MATCHES "^per_")
    string(REGEX REPLACE "^per_" "" ${name} "${${name}}")
  endif()
endmacro()

# collect modules from specified directories
# NB: must be called only once!
macro(pr_glob_modules)
  if(DEFINED PER_INITIAL_PASS)
    message(FATAL_ERROR "PER has already loaded its modules. Calling pr_glob_modules second time is not allowed.")
  endif()
  set(__directories_observed "")

  # collect modules
  set(PER_INITIAL_PASS ON)
  set(PER_PROCESSING_EXTRA_MODULES 0)
  foreach(__path ${ARGN})
    if("${__path}" STREQUAL "EXTRA")
      set(PER_PROCESSING_EXTRA_MODULES 1)
    endif()
    get_filename_component(__path "${__path}" ABSOLUTE)

    list(FIND __directories_observed "${__path}" __pathIdx)
    if(__pathIdx GREATER -1)
      message(FATAL_ERROR "The directory ${__path} is observed for PER modules second time.")
    endif()
    list(APPEND __directories_observed "${__path}")

    file(GLOB __prmodules RELATIVE "${__path}" "${__path}/*")

    if(__prmodules)
      list(SORT __prmodules)
      foreach(mod ${__prmodules})
        get_filename_component(__modpath "${__path}/${mod}" ABSOLUTE)
        if(EXISTS "${__modpath}/CMakeLists.txt")

          list(FIND __directories_observed "${__modpath}" __pathIdx)
          if(__pathIdx GREATER -1)
            message(FATAL_ERROR "The module from ${__modpath} is already loaded.")
          endif()
          list(APPEND __directories_observed "${__modpath}")

          add_subdirectory("${__modpath}" "${CMAKE_CURRENT_BINARY_DIR}/${mod}/.${mod}")
        else()
          # modules in tracker
          get_filename_component(__subpath "${__path}/${mod}" ABSOLUTE)
          file(GLOB __prsubmodules RELATIVE "${__subpath}" "${__subpath}/*")
          if(__prsubmodules)
            list(SORT __prsubmodules)
            foreach(submod ${__prsubmodules})
              get_filename_component(__submodpath "${__subpath}/${submod}" ABSOLUTE)
              if(EXISTS "${__submodpath}/CMakeLists.txt")
                list(FIND __directories_observed "${__submodpath}" __pathIdx)
                if(__pathIdx GREATER -1)
                  message(FATAL_ERROR "The module from ${__submodpath} is already loaded.")
                endif()
                list(APPEND __directories_observed "${__submodpath}")
                add_subdirectory("${__submodpath}" "${CMAKE_CURRENT_BINARY_DIR}/${submod}/.${submod}")
              endif()
            endforeach()
          endif()
        endif()
      endforeach()
    endif()
  endforeach()
  pr_clear_vars(__prmodules __directories_observed __path __modpath __pathIdx __prsubmodules __subpath __submodpath)

  # resolve dependencies
  __pr_resolve_dependencies()

  # create modules
  set(PER_INITIAL_PASS OFF PARENT_SCOPE)
  set(PER_INITIAL_PASS OFF)

    foreach(m ${PER_MODULES_BUILD})
      if(m MATCHES "^per_")
        string(REGEX REPLACE "^per_" "" __shortname "${m}")
        add_subdirectory("${PER_MODULE_${m}_LOCATION}" "${CMAKE_CURRENT_BINARY_DIR}/${__shortname}")
      else()
        message(WARNING "Check module name: ${m}")
        add_subdirectory("${PER_MODULE_${m}_LOCATION}" "${CMAKE_CURRENT_BINARY_DIR}/${m}")
      endif()
    endforeach()

  unset(__shortname)
endmacro()

# disables PER module with missing dependencies
function(__pr_module_turn_off the_module)
  list(REMOVE_ITEM PER_MODULES_DISABLED_AUTO "${the_module}")
  list(APPEND PER_MODULES_DISABLED_AUTO "${the_module}")
  list(REMOVE_ITEM PER_MODULES_BUILD "${the_module}")
  list(REMOVE_ITEM PER_MODULES_PUBLIC "${the_module}")
  set(HAVE_${the_module} OFF CACHE INTERNAL "Module ${the_module} can not be built in current configuration")

  set(PER_MODULES_DISABLED_AUTO "${PER_MODULES_DISABLED_AUTO}" CACHE INTERNAL "")
  set(PER_MODULES_BUILD "${PER_MODULES_BUILD}" CACHE INTERNAL "")
  set(PER_MODULES_PUBLIC "${PER_MODULES_PUBLIC}" CACHE INTERNAL "")
endfunction()

# sort modules by dependencies
function(__pr_sort_modules_by_deps __lst)
  pr_list_sort(${__lst})
  set(input ${${__lst}})
  set(result "")
  while(input)
    list(LENGTH input length_before)
    foreach (m ${input})
      # check if module is in the result already
      if (NOT ";${result};" MATCHES ";${m};")
        # scan through module dependencies...
        set(unresolved_deps_found FALSE)
        foreach (d ${PER_MODULE_${m}_CHILDREN} ${PER_MODULE_${m}_DEPS})
          # ... which are not already in the result and are enabled
          if ((NOT ";${result};" MATCHES ";${d};") AND HAVE_${d})
            set(unresolved_deps_found TRUE)
            break()
          endif()
        endforeach()
        # chek if all dependencies for this module has been resolved
        if (NOT unresolved_deps_found)
          list(APPEND result ${m})
          list(REMOVE_ITEM input ${m})
        endif()
      endif()
    endforeach()
    list(LENGTH input length_after)
    # check for infinite loop or unresolved dependencies
    if (NOT length_after LESS length_before)
      message(WARNING "Unresolved dependencies or loop in dependency graph (${length_after})\n"
        "Processed ${__lst}: ${${__lst}}\n"
        "Good modules: ${result}\n"
        "Bad modules: ${input}"
      )
      list(APPEND result ${input})
      break()
    endif()
  endwhile()
  set(${__lst} "${result}" PARENT_SCOPE)
endfunction()

# resolve dependensies
function(__pr_resolve_dependencies)
  foreach(m ${PER_MODULES_DISABLED_USER})
    set(HAVE_${m} OFF CACHE INTERNAL "Module ${m} will not be built in current configuration")
  endforeach()
  foreach(m ${PER_MODULES_BUILD})
    set(HAVE_${m} ON CACHE INTERNAL "Module ${m} will be built in current configuration")
  endforeach()

  # disable MODULES with unresolved dependencies
  set(has_changes ON)
  while(has_changes)
    set(has_changes OFF)
    foreach(m ${PER_MODULES_BUILD})
      set(__deps ${PER_MODULE_${m}_REQ_DEPS} ${PER_MODULE_${m}_PRIVATE_REQ_DEPS})
      while(__deps)
        pr_list_pop_front(__deps d)
        string(TOLOWER "${d}" upper_d)
        if(NOT (HAVE_${d} OR HAVE_${upper_d} OR TARGET ${d} OR EXISTS ${d}))
          if(d MATCHES "^per_") # TODO Remove this condition in the future and use HAVE_ variables only
            message(STATUS "Module ${m} disabled because ${d} dependency can't be resolved!")
            __pr_module_turn_off(${m})
            set(has_changes ON)
            break()
          else()
            message(STATUS "Assume that non-module dependency is available: ${d} (for module ${m})")
          endif()
        endif()
      endwhile()
    endforeach()
  endwhile()

#  message(STATUS "List of active modules: ${PER_MODULES_BUILD}")

  foreach(m ${PER_MODULES_BUILD})
    set(deps_${m} ${PER_MODULE_${m}_REQ_DEPS})
    foreach(d ${PER_MODULE_${m}_OPT_DEPS})
      if(NOT (";${deps_${m}};" MATCHES ";${d};"))
        if(HAVE_${d} OR TARGET ${d})
          list(APPEND deps_${m} ${d})
        endif()
      endif()
    endforeach()
#    message(STATUS "Initial deps of ${m} (w/o private deps): ${deps_${m}}")
  endforeach()

  # propagate dependencies
  set(has_changes ON)
  while(has_changes)
    set(has_changes OFF)
    foreach(m2 ${PER_MODULES_BUILD}) # transfer deps of m2 to m
      foreach(m ${PER_MODULES_BUILD})
        if((NOT m STREQUAL m2) AND ";${deps_${m}};" MATCHES ";${m2};")
          foreach(d ${deps_${m2}})
            if(NOT (";${deps_${m}};" MATCHES ";${d};"))
#              message(STATUS "  Transfer dependency ${d} from ${m2} to ${m}")
              list(APPEND deps_${m} ${d})
              set(has_changes ON)
            endif()
          endforeach()
        endif()
      endforeach()
    endforeach()
  endwhile()

  # process private deps
  foreach(m ${PER_MODULES_BUILD})
    foreach(d ${PER_MODULE_${m}_PRIVATE_REQ_DEPS})
      if(NOT (";${deps_${m}};" MATCHES ";${d};"))
        list(APPEND deps_${m} ${d})
      endif()
    endforeach()
    foreach(d ${PER_MODULE_${m}_PRIVATE_OPT_DEPS})
      if(NOT (";${deps_${m}};" MATCHES ";${d};"))
        if(HAVE_${d} OR TARGET ${d})
          list(APPEND deps_${m} ${d})
        endif()
      endif()
    endforeach()
  endforeach()

  pr_list_sort(PER_MODULES_BUILD)

  foreach(m ${PER_MODULES_BUILD})
    #message(STATUS "FULL deps of ${m}: ${deps_${m}}")
    set(PER_MODULE_${m}_DEPS ${deps_${m}})
    set(PER_MODULE_${m}_DEPS_EXT ${deps_${m}})
    pr_list_filterout(PER_MODULE_${m}_DEPS_EXT "^per_[^ ]+$")
    if(PER_MODULE_${m}_DEPS_EXT AND PER_MODULE_${m}_DEPS)
      list(REMOVE_ITEM PER_MODULE_${m}_DEPS ${PER_MODULE_${m}_DEPS_EXT})
    endif()
  endforeach()

  # reorder dependencies
  foreach(m ${PER_MODULES_BUILD})
    __pr_sort_modules_by_deps(PER_MODULE_${m}_DEPS)
    pr_list_sort(PER_MODULE_${m}_DEPS_EXT)

    set(LINK_DEPS ${PER_MODULE_${m}_DEPS})

    set(PER_MODULE_${m}_DEPS ${PER_MODULE_${m}_DEPS} CACHE INTERNAL "Flattened dependencies of ${m} module")
    set(PER_MODULE_${m}_DEPS_EXT ${PER_MODULE_${m}_DEPS_EXT} CACHE INTERNAL "Extra dependencies of ${m} module")
    set(PER_MODULE_${m}_DEPS_TO_LINK ${LINK_DEPS} CACHE INTERNAL "Flattened dependencies of ${m} module (for linker)")

#    message(STATUS "  module deps of ${m}: ${PER_MODULE_${m}_DEPS}")
#    message(STATUS "  module link deps of ${m}: ${PER_MODULE_${m}_DEPS_TO_LINK}")
#    message(STATUS "  extra deps of ${m}: ${PER_MODULE_${m}_DEPS_EXT}")
#    message(STATUS "")
  endforeach()

  __pr_sort_modules_by_deps(PER_MODULES_BUILD)

  set(PER_MODULES_PUBLIC        ${PER_MODULES_PUBLIC}        CACHE INTERNAL "List of PER modules marked for export")
  set(PER_MODULES_BUILD         ${PER_MODULES_BUILD}         CACHE INTERNAL "List of PER modules included into the build")
  set(PER_MODULES_DISABLED_AUTO ${PER_MODULES_DISABLED_AUTO} CACHE INTERNAL "List of PER modules implicitly disabled due to dependencies")
endfunction()

# setup include paths for the list of passed modules
macro(pr_target_include_modules target)
  foreach(d ${ARGN})
    if(d MATCHES "^per_" AND HAVE_${d})
      if (EXISTS "${PER_MODULE_${d}_LOCATION}/include")
        pr_target_include_directories(${target} "${PER_MODULE_${d}_LOCATION}/include")
        # Work arround to be able to build the modules without INTERFACE_INCLUDE_DIRECTORIES
        # that was only introduces since CMake 2.8.12
        if (CMAKE_VERSION VERSION_LESS 2.8.12)
          pr_target_include_directories(${target} "${PER_MODULE_${d}_INC_DEPS}")
        endif()
      endif()
    elseif(EXISTS "${d}")
      # FS keep external deps inc
      set(PER_MODULE_${the_module}_INC_DEPS "${PER_MODULE_${the_module}_INC_DEPS};${d}" CACHE INTERNAL "")
      pr_target_include_directories(${target} "${d}")
    endif()
  endforeach()
  pr_list_unique(PER_MODULE_${the_module}_INC_DEPS)
endmacro()

# setup include path for PER headers for specified module
# pr_module_include_directories(<extra include directories/extra include modules>)
macro(pr_module_include_directories)
  pr_target_include_directories(${the_module}
      "${PER_MODULE_${the_module}_LOCATION}/include"
      "${PER_MODULE_${the_module}_LOCATION}/src"
      )
  pr_target_include_modules(${the_module} ${PER_MODULE_${the_module}_DEPS} ${ARGN})
endmacro()

# sets header and source files for the current module
# NB: all files specified as headers will be installed
# Usage:
# ocv_set_module_sources([HEADERS] <list of files> [SOURCES] <list of files>)
macro(pr_set_module_sources)
  pr_debug_message("pr_set_module_sources(" ${ARGN} ")")

  set(PER_MODULE_${the_module}_HEADERS "")
  set(PER_MODULE_${the_module}_SOURCES "")

  foreach(f "HEADERS" ${ARGN})
    if(f STREQUAL "HEADERS" OR f STREQUAL "SOURCES")
      set(__filesvar "PER_MODULE_${the_module}_${f}")
    else()
      list(APPEND ${__filesvar} "${f}")
    endif()
  endforeach()

  # use full paths for module to be independent from the module location
  pr_convert_to_full_paths(PER_MODULE_${the_module}_HEADERS)

  if(${the_module} MATCHES per_core)
#    list(APPEND PER_MODULE_${the_module}_HEADERS "${PER_INCLUDE_DIR}/per/core/prConfig.h")
    list(APPEND PER_MODULE_${the_module}_HEADERS "${PER_INCLUDE_DIR}/per/per_modules.h")
  endif()


  set(PER_MODULE_${the_module}_HEADERS ${PER_MODULE_${the_module}_HEADERS} CACHE INTERNAL "List of header files for ${the_module}")
  set(PER_MODULE_${the_module}_SOURCES ${PER_MODULE_${the_module}_SOURCES} CACHE INTERNAL "List of source files for ${the_module}")
endmacro()

# finds and sets headers and sources for the standard PER module
# Usage:
# pr_glob_module_sources(<extra sources&headers in the same format as used in pr_set_module_sources>)
macro(pr_glob_module_sources)
  pr_debug_message("pr_glob_module_sources(" ${ARGN} ")")
  set(_argn ${ARGN})

  file(GLOB_RECURSE lib_srcs
       "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp"
  )
  file(GLOB_RECURSE lib_int_hdrs
       "${CMAKE_CURRENT_LIST_DIR}/src/*.hpp"
       "${CMAKE_CURRENT_LIST_DIR}/src/*.h"
  )
  file(GLOB lib_hdrs
       "${CMAKE_CURRENT_LIST_DIR}/include/per/*.h"
       "${CMAKE_CURRENT_LIST_DIR}/include/per/${name}/*.h"
  )

  pr_source_group("Src" DIRBASE "${CMAKE_CURRENT_LIST_DIR}/src" FILES ${lib_srcs} ${lib_int_hdrs})
  pr_source_group("Include" DIRBASE "${CMAKE_CURRENT_LIST_DIR}/include" FILES ${lib_hdrs})

  pr_set_module_sources(${_argn} HEADERS ${lib_hdrs}
                        SOURCES ${lib_srcs} ${lib_int_hdrs})
endmacro()

# finds and copy data from a source to a destination
# Usage:
# pr_glob_module_data(<source> <destination>)
macro(pr_glob_module_copy_data src dst)
  set(__data "")
  file(GLOB_RECURSE __data
       "${CMAKE_CURRENT_LIST_DIR}/${src}"
  )

  foreach(__d ${__data})
    file(COPY ${__d}
       DESTINATION "${PER_BINARY_DIR}/${dst}"
       FILE_PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
       OWNER_WRITE
    )

    # install
    if(UNIX)
      set(__install_dst "${CMAKE_INSTALL_DATAROOTDIR}/per-${PER_VERSION}/${dst}")
    else()
      set(__install_dst "${dst}")
    endif()

    install(FILES ${__d}
       DESTINATION "${__install_dst}"
       PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
       OWNER_WRITE
    )
  endforeach()
endmacro()

# creates old headers for compat with previous releases in include/per
# Usage:
#   pr_create_compat_headers(<list of new headers>)
macro(pr_create_compat_headers)
  pr_debug_message("pr_create_compat_headers(" ${ARGN} ")")

  set(PER_HEADER_CONTENT_CONFIGMAKE "")

  foreach(h ${ARGN})
    get_filename_component(__h_name_we ${h} NAME_WE)
    get_filename_component(__h_name ${h} NAME)
    set(PER_HEADER_CONTENT_CONFIGMAKE "#ifndef __${__h_name_we}_h_\n#define __${__h_name_we}_h_\n\n#include <per/${name}/${__h_name}>\n\n#endif\n")
    set(__compat_header_dst "${PER_INCLUDE_DIR}/per/${__h_name_we}.h")
    configure_file("${PER_SOURCE_DIR}/cmake/templates/prHeader.h.in" ${__compat_header_dst})
  endforeach()

  unset(__h_name_we)
  unset(__h_name)
  unset(__compat_header_dst)
endmacro()

# creates headers for modules include/per/<module>/pr<module>.h
# Usage:
#   pr_create_global_module_header(<module>)
macro(pr_create_global_module_header module)
  pr_debug_message("pr_create_global_module_header(" ${module} ")")

  set(__name ${module})
  pr_short_module_name(__name)
  set(__module_header_dst "${PER_INCLUDE_DIR}/per/${__name}/${module}.h")
  set(PER_HEADER_CONTENT_CONFIGMAKE "#ifndef __${module}_h_\n#define __${module}_h_\n")

  # when core, include also prConfig.h
#  if(__name MATCHES "core")
#    set(PER_HEADER_CONTENT_CONFIGMAKE "${PER_HEADER_CONTENT_CONFIGMAKE}\n#include <per/${__name}/prConfig.h>")
#  endif()

  # include the modules we depend on
  if(PER_MODULE_${module}_REQ_DEPS)
    foreach(dep ${PER_MODULE_${module}_REQ_DEPS})
      pr_short_module_name(dep)
    set(PER_HEADER_CONTENT_CONFIGMAKE "${PER_HEADER_CONTENT_CONFIGMAKE}\n#include <per/${dep}/pr${dep}.h>")
    endforeach()
  endif()

  foreach(h ${PER_MODULE_${module}_HEADERS})
    get_filename_component(__h_name_we ${h} NAME_WE)
    get_filename_component(__h_name ${h} NAME)
    set(PER_HEADER_CONTENT_CONFIGMAKE "${PER_HEADER_CONTENT_CONFIGMAKE}\n#include <per/${__name}/${__h_name}>")
  endforeach()

  set(PER_HEADER_CONTENT_CONFIGMAKE "${PER_HEADER_CONTENT_CONFIGMAKE}\n\n#endif\n")
  configure_file("${PER_SOURCE_DIR}/cmake/templates/prHeader.h.in" ${__module_header_dst})

  unset(__h_name_we)
  unset(__h_name)
  unset(__module_header_dst)
endmacro()

# creates PER module in current folder
# creates new target, configures standard dependencies, compilers flags, install rules
# Usage:
#   pr_create_module(<extra link dependencies>)
#   pr_create_module()
macro(pr_create_module)
  pr_debug_message("pr_create_module(" ${ARGN} ")")
  set(PER_MODULE_${the_module}_LINK_DEPS "${PER_MODULE_${the_module}_LINK_DEPS};${ARGN}" CACHE INTERNAL "")
  _pr_create_module(${ARGN})
  set(the_module_target ${the_module})
endmacro()

macro(_pr_create_module)
  pr_create_compat_headers(${PER_MODULE_${the_module}_HEADERS})
  pr_create_global_module_header(${the_module})

  pr_add_library(${the_module} ${PER_MODULE_TYPE} ${PER_MODULE_${the_module}_HEADERS} ${PER_MODULE_${the_module}_SOURCES})

  pr_target_link_libraries(${the_module} ${PER_MODULE_${the_module}_DEPS_TO_LINK})
  #pr_target_link_libraries(${the_module} LINK_INTERFACE_LIBRARIES ${PER_MODULE_${the_module}_DEPS_TO_LINK})
  pr_target_link_libraries(${the_module} ${PER_MODULE_${the_module}_DEPS_EXT} ${PER_LINKER_LIBS} ${ARGN})

  add_dependencies(per_modules ${the_module})

  if(ENABLE_SOLUTION_FOLDERS)
    set_target_properties(${the_module} PROPERTIES FOLDER "modules")
  endif()

  set_target_properties(${the_module} PROPERTIES
    OUTPUT_NAME "${the_module}${PER_DLLVERSION}"
    DEBUG_POSTFIX "${PER_DEBUG_POSTFIX}"
    ARCHIVE_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH}
    LIBRARY_OUTPUT_DIRECTORY ${LIBRARY_OUTPUT_PATH}
    RUNTIME_OUTPUT_DIRECTORY ${BINARY_OUTPUT_PATH}
  )

  set_property(TARGET ${the_module} APPEND PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${PER_MODULE_${the_module}_INC_DEPS}
  )

  # For dynamic link numbering convenions
  if(NOT ANDROID)
    # Android SDK build scripts can include only .so files into final .apk
    # As result we should not set version properties for Android
    set_target_properties(${the_module} PROPERTIES
      VERSION ${PER_VERSION}
      SOVERSION ${PER_VERSION_MAJOR}.${PER_VERSION_MINOR}
    )
  endif()

  if((NOT DEFINED PER_MODULE_TYPE AND BUILD_SHARED_LIBS)
      OR (DEFINED PER_MODULE_TYPE AND PER_MODULE_TYPE STREQUAL SHARED))
    set_target_properties(${the_module} PROPERTIES COMPILE_DEFINITIONS per_EXPORTS)
    set_target_properties(${the_module} PROPERTIES DEFINE_SYMBOL per_EXPORTS)
  endif()

  if(MSVC)
    if(CMAKE_CROSSCOMPILING)
      set_target_properties(${the_module} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:secchk")
    endif()
    set_target_properties(${the_module} PROPERTIES LINK_FLAGS "/NODEFAULTLIB:libc /DEBUG")
  endif()

  pr_install_target(${the_module} EXPORT PERModules OPTIONAL
    RUNTIME DESTINATION ${PER_BIN_INSTALL_PATH} COMPONENT libs
    LIBRARY DESTINATION ${PER_LIB_INSTALL_PATH} COMPONENT libs
    ARCHIVE DESTINATION ${PER_LIB_INSTALL_PATH} COMPONENT dev
    )

  foreach(m ${PER_MODULE_${the_module}_CHILDREN} ${the_module})
    # only "public" headers need to be installed
    if(PER_MODULE_${m}_HEADERS AND ";${PER_MODULES_PUBLIC};" MATCHES ";${m};")
      foreach(hdr ${PER_MODULE_${m}_HEADERS})
        string(REGEX REPLACE "^.*per/" "per/" hdr2 "${hdr}")
        if(NOT hdr2 MATCHES "per/${m}/private.*" AND hdr2 MATCHES "^(per/?.*)/[^/]+.h(..)?$" )
          install(FILES ${hdr} OPTIONAL DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_MATCH_1}" COMPONENT dev)
        endif()
      endforeach()
    endif()
  endforeach()

endmacro()


# short command for adding simple PER module
# see pr_add_module for argument details
# Usage:
# pr_define_module(module_name  [INTERNAL] [REQUIRED] [<list of dependencies>] [OPTIONAL <list of optional dependencies>])
macro(pr_define_module module_name)
  pr_debug_message("pr_define_module(" ${module_name} ${ARGN} ")")
  set(_argn ${ARGN})

  pr_add_module(${module_name} ${_argn})
  pr_glob_module_sources()
  pr_module_include_directories()
  pr_create_module()
endmacro()

# ensures that all passed modules are available
# sets PR_DEPENDENCIES_FOUND variable to TRUE/FALSE
macro(pr_check_dependencies)
  set(PR_DEPENDENCIES_FOUND TRUE)
  foreach(d ${ARGN})
    if(d MATCHES "^per_[^ ]+$" AND NOT HAVE_${d})
      set(PR_DEPENDENCIES_FOUND FALSE)
      break()
    endif()
  endforeach()
endmacro()

# auxiliary macro to parse arguments of pr_add_tests commands
macro(__pr_parse_test_sources tests_type)
  set(PER_${tests_type}_${the_module}_SOURCES "")
  set(PER_${tests_type}_${the_module}_SOURCES_EXCLUDE "")
  set(PER_${tests_type}_${the_module}_DEPS "")
  set(PER_${tests_type}_${the_module}_CTEST_EXCLUDE_FOLDER "")
  set(__file_group_name "")
  set(__file_group_sources "")
  foreach(arg "DEPENDS_ON" ${ARGN} "FILES")
    if(arg STREQUAL "FILES")
      set(__currentvar "__file_group_sources")
      if(__file_group_name AND __file_group_sources)
        source_group("${__file_group_name}" FILES ${__file_group_sources})
        list(APPEND PER_${tests_type}_${the_module}_SOURCES ${__file_group_sources})
      endif()
      set(__file_group_name "")
      set(__file_group_sources "")
    elseif(arg STREQUAL "DEPENDS_ON")
      set(__currentvar "PER_${tests_type}_${the_module}_DEPS")
    elseif(" ${__currentvar}" STREQUAL " __file_group_sources" AND NOT __file_group_name) # spaces to avoid CMP0054
      set(__file_group_name "${arg}")
    elseif(arg STREQUAL "CTEST_EXCLUDE_PATH")
      set(__currentvar "PER_${tests_type}_${the_module}_CTEST_EXCLUDE_FOLDER")
    elseif(arg STREQUAL "SOURCES_EXCLUDE")
      set(__currentvar "PER_${tests_type}_${the_module}_SOURCES_EXCLUDE")
    else()
      list(APPEND ${__currentvar} "${arg}")
    endif()
  endforeach()
  unset(__file_group_name)
  unset(__file_group_sources)
  unset(__currentvar)
endmacro()

# this is a command for adding PER tests to the module
# pr_add_tests([FILES <source group name> <list of sources>]
#              [FILES_EXCLUDE <list of sources>]
#              [DEPENDS_ON] <list of extra dependencies>
#              [CTEST_EXCLUDE_FOLDER] <list of folder to exclude from ctest>)
macro(pr_add_tests)
  pr_debug_message("pr_add_tests(" ${ARGN} ")")

  set(test_path "${CMAKE_CURRENT_LIST_DIR}/test")
  if(BUILD_TESTS AND EXISTS "${test_path}")
    __pr_parse_test_sources(TEST ${ARGN})

    set(__exclude_ctest "")
    foreach(__folder ${PER_TEST_${the_module}_CTEST_EXCLUDE_FOLDER} )
      file(GLOB_RECURSE __files "${CMAKE_CURRENT_LIST_DIR}/test/${__folder}/*.cpp")
      list(APPEND __exclude_ctest ${__files})
    endforeach()
    set(__exclude_sources "")
    foreach(__source ${PER_TEST_${the_module}_SOURCES_EXCLUDE} )
      file(GLOB __files "${CMAKE_CURRENT_LIST_DIR}/test/${__source}")
      list(APPEND __exclude_sources ${__files})
    endforeach()

    set(test_deps ${the_module} ${PER_MODULE_${the_module}_DEPS})

    foreach(d ${PER_TEST_${the_module}_DEPS})
      list(APPEND test_deps ${d})
      list(APPEND test_deps ${PER_MODULE_${d}_DEPS})
      # Work arround to be able to build the modules without INTERFACE_INCLUDE_DIRECTORIES
      # that was only introduces since CMake 2.8.12
      if(CMAKE_VERSION VERSION_LESS 2.8.12)
        list(APPEND test_deps "${PER_MODULE_${__m}_INC_DEPS}")
      endif()
    endforeach()

    pr_check_dependencies(${test_deps})
    if(PR_DEPENDENCIES_FOUND)
      if(NOT PER_TEST_${the_module}_SOURCES)
        file(GLOB_RECURSE test_srcs "${test_path}/*.cpp")
        pr_source_group("Src" DIRBASE "${test_path}" FILES ${test_srcs})
        set(PER_TEST_${the_module}_SOURCES ${test_srcs})
      endif()

      foreach(t ${PER_TEST_${the_module}_SOURCES})
        # check if source is not in exclude list
        list(FIND __exclude_sources ${t} __to_exclude_from_sources)
        if(${__to_exclude_from_sources} EQUAL -1)
          # Compute the name of the binary to create
          get_filename_component(the_target ${t} NAME_WE)
          # From source compile the binary and add link rules
          pr_add_executable(${the_target} ${t})
          pr_target_include_modules(${the_target} ${test_deps})
          pr_target_link_libraries(${the_target} ${test_deps} ${PER_MODULE_${the_module}_DEPS} ${PER_LINKER_LIBS})

          # ctest only if not in the exclude list
          list(FIND __exclude_ctest ${t} __to_exclude_from_ctest)
          if(${__to_exclude_from_ctest} EQUAL -1)
            add_test(${the_target} ${the_target} -c ${OPTION_TO_DESACTIVE_DISPLAY})
          endif()
          # TODO FS add per_test_${name} target to group all the tests
          add_dependencies(per_tests ${the_target})
          if(ENABLE_SOLUTION_FOLDERS)
            set_target_properties(${the_target} PROPERTIES FOLDER "tests")
          endif()
        endif()
      endforeach()

    else(PR_DEPENDENCIES_FOUND)
      # TODO: warn about unsatisfied dependencies
    endif(PR_DEPENDENCIES_FOUND)

  endif()
endmacro()

# setup include paths for the list of passed modules
macro(pr_include_modules)
  foreach(d ${ARGN})
    if(d MATCHES "^per_" AND HAVE_${d})
      if (EXISTS "${PER_MODULE_${d}_LOCATION}/include")
        pr_include_directories("${PER_MODULE_${d}_LOCATION}/include")
      endif()
    elseif(EXISTS "${d}")
      pr_include_directories("${d}")
    endif()
  endforeach()
endmacro()

# same as previous but with dependencies
macro(pr_include_modules_recurse)
  pr_include_modules(${ARGN})
  foreach(d ${ARGN})
    if(d MATCHES "^per_" AND HAVE_${d} AND DEFINED PER_MODULE_${d}_DEPS)
      foreach (sub ${PER_MODULE_${d}_DEPS})
        pr_include_modules(${sub})
      endforeach()
    endif()
  endforeach()
endmacro()

# This is a command to configure files as include headers of the corresponding module.
# pr_add_config_file(<list of header config files>)
#
# If the input config filename is suffixed by .in or .cmake the suffix is removed
# in the configured file.
#
# Example:
#   pr_add_config_file(cmake/template/prConfigModule.h.in)
#   creates include/per/module_name/prConfigModule.h
macro(pr_add_config_file)
  foreach(d ${ARGN})
    # Removes first "/" if it exists
    string(FIND ${d} "/" FIRST_SEPARATOR_POS)
    if(${FIRST_SEPARATOR_POS} EQUAL 0)
      string(SUBSTRING ${d} 1 -1 d)
    endif()

    # Find start of file name
    string(FIND ${d} "/" LAST_SEPARATOR_POS REVERSE)
    if(${LAST_SEPARATOR_POS} EQUAL -1)
      set(START 0)
    else()
      math(EXPR START "${LAST_SEPARATOR_POS}+1")
    endif()

    # Save entire path
    set(FILENAME_CONFIG ${d})

    # Find file name
    string(FIND ${d} "." EXTENSION_POS REVERSE)

    if(${EXTENSION_POS} EQUAL -1)
      string(SUBSTRING ${d} ${START} -1 FILENAME_CONFIG_SHORT)
    else()
      string(SUBSTRING ${d} ${EXTENSION_POS} -1 EXT_CONFIG_FILE)
      if(EXT_CONFIG_FILE MATCHES ".cmake" OR EXT_CONFIG_FILE MATCHES ".in")
        math(EXPR LENGTH "${EXTENSION_POS} - ${START}")
        string(SUBSTRING ${d} ${START} ${LENGTH} FILENAME_CONFIG_SHORT)
      else()
        string(SUBSTRING ${d} ${START} -1 FILENAME_CONFIG_SHORT)
      endif()
    endif()

    set(MODULE_NAME ${the_module})
    if(MODULE_NAME MATCHES "^per_")
      string(REGEX REPLACE "^per_" "" MODULE_NAME "${MODULE_NAME}")
    endif()

    configure_file("${PER_MODULE_${the_module}_LOCATION}/${FILENAME_CONFIG}" "${PER_INCLUDE_DIR}/per/${MODULE_NAME}/${FILENAME_CONFIG_SHORT}")

    pr_create_compat_headers("${PER_INCLUDE_DIR}/per/${MODULE_NAME}/${FILENAME_CONFIG_SHORT}")

    install(FILES "${PER_INCLUDE_DIR}/per/${MODULE_NAME}/${FILENAME_CONFIG_SHORT}"
      DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/per/${MODULE_NAME}
      COMPONENT dev
    )

  endforeach()
endmacro()

# This is a command to add a list of paths associated to the corresponding module
# to the CMAKE_MODULE_PATH global var to find specific cmake material
# pr_add_cmake_module_path(<list of cmake module paths>)
# Example:
#   pr_add_cmake_module_path(cmake)
#   Appends the cmake full path to CMAKE_MODULE_PATH var.
macro(pr_add_cmake_module_path)
  foreach(d ${ARGN})
    # Removes first "/" if it exists
    string(FIND ${d} "/" FIRST_SEPARATOR_POS)
    if(${FIRST_SEPARATOR_POS} EQUAL 0)
      string(SUBSTRING ${d} 1 -1 d)
    endif()
    if(EXISTS "${PER_MODULE_${the_module}_LOCATION}/${d}")
      list(APPEND CMAKE_MODULE_PATH "${PER_MODULE_${the_module}_LOCATION}/${d}")
    endif()
  endforeach()
endmacro()
