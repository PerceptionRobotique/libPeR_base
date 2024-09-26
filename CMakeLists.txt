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
# ViSP library. Detect third party libraries (ViSP, ...)
#
# Authors:
# Guillaume Caron
#
#############################################################################

# Need to be befor project(PER) to work
if(WIN32)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Installation Directory")
endif()

project(PER C CXX)

cmake_minimum_required(VERSION 2.8.3)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

include(cmake/PERUtils.cmake)

pr_clear_vars(PERModules_TARGETS)

#-----------------------------------------------------------------------------
# PER version number.  An even minor number corresponds to releases.
set(PER_VERSION_MAJOR "0")
set(PER_VERSION_MINOR "4")
set(PER_VERSION_PATCH "0")
set(PER_VERSION "${PER_VERSION_MAJOR}.${PER_VERSION_MINOR}.${PER_VERSION_PATCH}")
set(PROJECT_FULL_VERSION ${PER_VERSION})
# Package revision number
set(PER_REVISION "1")

find_file(GNU_INSTALL_DIRS_FROM_CMAKE NAMES GNUInstallDirs.cmake PATHS ${CMAKE_ROOT}/Modules)
mark_as_advanced(GNU_INSTALL_DIRS_FROM_CMAKE)
if(GNU_INSTALL_DIRS_FROM_CMAKE)
  include(${CMAKE_ROOT}/Modules/GNUInstallDirs.cmake)
else()
  include(cmake/GNUInstallDirs.cmake)
endif()

# the include directory we depend on
set(PER_INCLUDE_DIR ${PER_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR})
set(PER_DOC_DIR "${PER_BINARY_DIR}/doc")

# The location in which to install PER libraries.
set(LIBRARY_OUTPUT_PATH ${PER_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
# The location in which to install some PER binaries.
set(BINARY_OUTPUT_PATH ${PER_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

include_directories(${PER_BINARY_DIR}/include)

if(WIN32)
  # Postfix of .lib and .dll
  set(PER_DEBUG_POSTFIX "d")
  set(PER_DLLVERSION "${PER_VERSION_MAJOR}${PER_VERSION_MINOR}${PER_VERSION_PATCH}")
else()
  set(PER_DEBUG_POSTFIX "")
  set(PER_DLLVERSION "")
endif()

# ----------------------------------------------------------------------------
# Handle always full RPATH
# http://www.cmake.org/Wiki/CMake_RPATH_handling
# ----------------------------------------------------------------------------

# avoid CMP0042 warning
if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif()

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
endif("${isSystemDir}" STREQUAL "-1")

# ----------------------------------------------------------------------------
#  Path for additional contrib modules
# ----------------------------------------------------------------------------
set(PER_CONTRIB_MODULES_PATH "" CACHE PATH "Where to look for additional contrib PER modules")

# Get the OS
set(OS ${CMAKE_SYSTEM_NAME})

#--------------------------------------------------------------------
# Option management
#--------------------------------------------------------------------

# Choose static or shared libraries.
PR_OPTION(BUILD_SHARED_LIBS  "" "" "Build PER shared libraries (.dll/.so) instead of static ones (.lib/.a)" "" ON)
# Build examples as an option.
PR_OPTION(BUILD_EXAMPLES  "" "" "Build PER examples" "" ON)
# Build examples as an option.
PR_OPTION(BUILD_TESTS  "" "" "Build PER tests" "" ON)
PR_OPTION(BUILD_COVERAGE "" "" "Enables test coverage" "" OFF IF (BUILD_TESTS AND CMAKE_COMPILER_IS_GNUCXX AND NOT BUILD_SHARED_LIBS AND CMAKE_BUILD_TYPE MATCHES "Debug"))

# Build demos as an option.
PR_OPTION(BUILD_DEMOS  "" "" "Build PER demos" "" ON)
# Build demos as an option.
PR_OPTION(BUILD_TUTORIALS  "" "" "Build PER tutorials" "" ON)
# Build deprecated functions as an option.
PR_OPTION(BUILD_DEPRECATED_FUNCTIONS  "" "" "Build deprecated functionalities" "" ON)
# Debug and trace cflags
PR_OPTION(ACTIVATE_DEBUG_TRACE  "" "" "Enable debug and trace printings" "" ON)

PR_OPTION(BUILD_WITH_STATIC_CRT  "" "" "Enables use of staticaly linked CRT for staticaly linked PER" "" ON IF MSVC)

# Note that it is better to set MOMENTS_COMBINE_MATRICES to OFF
PR_OPTION(MOMENTS_COMBINE_MATRICES  "" "" "Use linear combination of matrices instead of linear combination of moments to compute interaction matrices." "MOMENTS_COMBINE_MATRICES" OFF)
PR_OPTION(ENABLE_TEST_WITHOUT_DISPLAY  "" "" "Don't use display features when testing" "ENABLE_TEST_WITHOUT_DISPLAY" ON)

if(ENABLE_TEST_WITHOUT_DISPLAY)
  set(OPTION_TO_DESACTIVE_DISPLAY "-d")
endif()

# Check for MIS’s robot drivers
#PR_CHECK_PACKAGE(IDS)
#PR_CHECK_PACKAGE(UPJV)
# Find IsNaN
#PR_CHECK_PACKAGE(IsNaN)
# Find IsInf
#PR_CHECK_PACKAGE(IsInf)
# Find Round
#PR_CHECK_PACKAGE(Round)
# for pioneer
#PR_CHECK_PACKAGE(RT)
#PR_CHECK_PACKAGE(DL)

PR_OPTION(USE_TX60 "" "" "Include Staubli TX60 robot support" "" ON IF UPJV_FOUND)
#PR_OPTION(USE_IDS      IDS      "" "Include IDS support"          "" ON)
PR_OPTION(USE_OPENMP      OpenMP      "" "Include openmp support"          "" ON)
#PR_OPTION(USE_CPP11       CPP11       "" "Include c++11 support"           "" OFF)
PR_OPTION(USE_VISP        VISP        REQUIRED "Include VISP support"         "VISP_DIR" ON)
#PR_OPTION(USE_LIBFREENECT LIBFREENECT ""    "Include libfreenect support"  "" ON)
#PR_OPTION(USE_LIBUSB_1    LIBUSB_1    ""    "Include libusb-1 support"     "" ON)
#PR_OPTION(USE_PTHREAD     PTHREAD     ""    "Include pthread support"      "" ON)
if(WIN32)
    PR_OPTION(USE_XML2        XML2        ""    "Include xml support"          "" ON)
    PR_OPTION(USE_ICONV ICONV "" "Include iconv support" "" ON)
endif()
PR_OPTION(USE_OPENCV      OpenCV      QUIET "Include OpenCV support" "OpenCV_DIR;OpenCV_FOUND;OPENCV_FOUND" ON)
# To control Pioneer mobile robots, under UNIX we need Aria, pthread, rt and dl 3rd party libraries
#PR_OPTION(USE_ARIA        ARIA         ""    "Include aria support"         "" ON)
#PR_OPTION(USE_RT          RT           ""    "Include rt support"           "" ON)
#PR_OPTION(USE_DL          DL           ""    "Include dl support"           "" ON)

#----------------------------------------------------------------------
# Try to find doxygen for documentation generation
# Use "make per_doc" target to generate the documentation
#----------------------------------------------------------------------
find_package(Doxygen)
if(DOXYGEN_FOUND)
  set(PER_HAVE_DOXYGEN "yes")        # for header prConfig.h
  set(PER_HAVE_DOXYGEN_FOUND "yes")  # for PER-third-party.txt
  ## we need latex for doxygen because of the formulas
  find_package(LATEX)
  if(NOT LATEX_COMPILER)
    message(STATUS "latex command LATEX_COMPILER not found but usually required. You will probably get warnings and user interaction on doxy run.")
  endif()
  if(NOT MAKEINDEX_COMPILER)
    message(STATUS "makeindex command MAKEINDEX_COMPILER not found but usually required.")
  endif()
  if(NOT DVIPS_CONVERTER)
    message(STATUS "dvips command DVIPS_CONVERTER not found but usually required.")
  endif()

  configure_file(${PER_SOURCE_DIR}/doc/config-doxygen.in
    ${PER_DOC_DIR}/config-doxygen
    @ONLY )

  configure_file(${PER_SOURCE_DIR}/doc/mainpage.doc.in
    ${PER_DOC_DIR}/mainpage.doc
    @ONLY )
else()
  set(PER_HAVE_DOXYGEN "no")        # for header prConfig.h
  set(PER_HAVE_DOXYGEN_FOUND "no")  # for PER-third-party.txt
endif()

# ----------------------------------------------------------------------------
# Solution folders:
# ----------------------------------------------------------------------------
if(MSVC_IDE OR CMAKE_GENERATOR MATCHES Xcode)
  option(ENABLE_SOLUTION_FOLDERS "Solution folder in Visual Studio or in other IDEs" ON)
endif()

if(ENABLE_SOLUTION_FOLDERS)
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
  set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakeTargets")
endif()

# Extra PER targets: uninstall, etc.
include(cmake/PERExtraTargets.cmake)

include(CheckLibraryExists)

# ----------------------------------------------------------------------------
# Check for system libs
# ----------------------------------------------------------------------------
if(UNIX)
  # try to found -lm requested on some platforms to link with X11
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND PER_LINKER_LIBS ${M_LIBRARY})
  endif()
  # try to found -lsocket -lnsl requested for vpNetwork and vpSickLDMRS
  find_library(SOCKET_LIBRARY NAMES socket)
  find_library(NSL_LIBRARY NAMES nsl)
  mark_as_advanced(SOCKET_LIBRARY NSL_LIBRARY)
  if (SOCKET_LIBRARY)
    list(APPEND PER_LINKER_LIBS ${SOCKET_LIBRARY})
  endif()
  if (NSL_LIBRARY)
    list(APPEND PER_LINKER_LIBS ${NSL_LIBRARY})
  endif()
endif()

#----------------------------------------------------------------------
# Add definitions
#----------------------------------------------------------------------
# With Visual Studio 2005, Microsoft deprecates the standard C library, for
# example fopen() and sprintf(), to non-portable functions fopen_s() and
# sprintf_s(). These functions are considered by Microsoft more secure. This is
# a worthwhile exercise ! The use of these deprecated functions causes a lot of
# warnings. To suppress it, we add the _CRT_SECURE_NO_DEPRECATE preprocessor
# definition
if(WIN32 AND MSVC)
  add_definitions("-D_CRT_SECURE_NO_DEPRECATE")
endif()

#----------------------------------------------------------------------
# Use statically or dynamically linked CRT?
# Default: dynamic
#----------------------------------------------------------------------
#if(MSVC)
#   include(cmake/PERCRTLinkage.cmake)
#endif(MSVC)

#----------------------------------------------------------------------
# Platform specific
#----------------------------------------------------------------------
include(cmake/PERDetectPlatform.cmake)

# Set the path where to install the lib
if(WIN32)
  if(DEFINED PER_RUNTIME AND DEFINED PER_ARCH)
    set(PER_INSTALL_BINARIES_PREFIX "${PER_ARCH}/${PER_RUNTIME}/")
  else()
    message(STATUS "Can't detect runtime and/or arch")
    set(PER_INSTALL_BINARIES_PREFIX "")
  endif()
else()
  set(PER_INSTALL_BINARIES_PREFIX "")
endif()

# where to install the library
if(WIN32)
  if(PER_STATIC)
    set(PER_LIB_INSTALL_PATH   "${PER_INSTALL_BINARIES_PREFIX}static${CMAKE_INSTALL_LIBDIR}")
  else()
    set(PER_LIB_INSTALL_PATH   "${PER_INSTALL_BINARIES_PREFIX}${CMAKE_INSTALL_LIBDIR}")
  endif()
  set(PER_BIN_INSTALL_PATH   "${PER_INSTALL_BINARIES_PREFIX}${CMAKE_INSTALL_BINDIR}")
else()
  set(PER_LIB_INSTALL_PATH     ${CMAKE_INSTALL_LIBDIR})
  set(PER_BIN_INSTALL_PATH     ${CMAKE_INSTALL_BINDIR})
endif()


#----------------------------------------------------------------------
# Configure the files that depend on the build <binary dir> or
# installation <install dir> usage. This is the case of:
# prConfig.h
#----------------------------------------------------------------------

# case 1: when libPR is build with make; files are used in <binary dir>
#--------------
set(PER_ROOT_DIR_DATA_CONFIGCMAKE "${PER_BINARY_DIR}")

#set(PER_SCENES_DIR ${PER_ROOT_DIR_DATA_CONFIGCMAKE}/data/wireframe-simulator)
#set(PER_ROBOT_ARMS_DIR ${PER_ROOT_DIR_DATA_CONFIGCMAKE}/data/robot-simulator)

# case 2: when libPR is build with make install; files are used in <install dir>
#--------------
IF(UNIX)
  set(PER_ROOT_DIR_DATA_CONFIGCMAKE "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATAROOTDIR}/per-${PER_VERSION}")
ELSE()
  set(PER_ROOT_DIR_DATA_CONFIGCMAKE "${CMAKE_INSTALL_PREFIX}")
ENDIF()

#----------------------------------------------------------------------
# Create and install per-config.1.gz man page
#----------------------------------------------------------------------
#if(UNIX)
#  FIND_PROGRAM(GZIP gzip)
#  file(MAKE_DIRECTORY ${PER_BINARY_DIR}/doc/man/man1)
#  ADD_CUSTOM_COMMAND(
#    OUTPUT ${PER_BINARY_DIR}/doc/man/man1/per-config.1.gz
#    COMMAND ${GZIP} --best -c ${CMAKE_CURRENT_SOURCE_DIR}/doc/man/man1/per-config.1 > ${PER_BINARY_DIR}/doc/man/man1/per-config.1.gz
#    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/doc/man/man1/per-config.1
#  )
#  ADD_CUSTOM_TARGET(man ALL
#    DEPENDS ${PER_BINARY_DIR}/doc/man/man1/per-config.1.gz
#  )
#  INSTALL(FILES
#    ${PER_BINARY_DIR}/doc/man/man1/per-config.1.gz
#    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/man/man1
#    PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_WRITE
#    COMPONENT dev
#  )
#  MARK_AS_ADVANCED(GZIP)
#endif()

#----------------------------------------------------------------------
# Modules
#----------------------------------------------------------------------
include(cmake/PERModule.cmake)

# process subdirectories
add_subdirectory(modules)

#-----------------------------------------------------------------------------
# Add extra compilation flags under UNIX
#-----------------------------------------------------------------------------
include(cmake/AddExtraCompilationFlags.cmake) # should be called after FindOpenMP
pr_add_extra_compilation_flags()

#-------------------------------------------------------------------------------
# specific things that need to be updated in prConfig.h and PER-third-party.txt
#-------------------------------------------------------------------------------
PR_SET(PER_HAVE_OPENMP      TRUE IF USE_OPENMP)
PR_SET(PER_HAVE_OPENCV      TRUE IF (BUILD_MODULE_per_core AND USE_OPENCV))
PR_SET(PER_HAVE_VISP        TRUE IF (BUILD_MODULE_per_core AND USE_VISP))
PR_SET(PER_HAVE_PTHREAD     TRUE IF (BUILD_MODULE_per_core AND USE_PTHREAD))
PR_SET(PER_HAVE_XML2        TRUE IF (BUILD_MODULE_per_core AND USE_XML2))
PR_SET(PER_HAVE_XML2        TRUE IF (BUILD_MODULE_per_io AND USE_XML2))
PR_SET(PER_HAVE_ICONV TRUE IF (BUILD_MODULE_per_io AND USE_ICONV))
PR_SET(PER_HAVE_TX60    TRUE IF (BUILD_MODULE_per_robot AND USE_TX60))
#PR_SET(PER_HAVE_PIONEER     TRUE IF (BUILD_MODULE_per_robot AND USE_ARIA))
if(BUILD_MODULE_per_robot AND USE_ARIA)
  if(UNIX AND USE_PTHREAD AND RT_FOUND AND DL_FOUND)
    set(PER_HAVE_PIONEER TRUE)
  elseif(NOT UNIX)
    set(PER_HAVE_PIONEER TRUE)
  endif()
endif()

PR_SET(PER_HAVE_LIBFREENECT TRUE IF (BUILD_MODULE_per_sensor AND USE_LIBFREENECT))
PR_SET(PER_HAVE_LIBUSB_1    TRUE IF (BUILD_MODULE_per_sensor AND USE_LIBUSB_1))

PR_SET(PER_BUILD_SHARED_LIBS TRUE IF BUILD_SHARED_LIBS) # for header prConfig.h

# Check if libfreenect dependencies (ie libusb-1.0 and libpthread) are available
if(USE_LIBFREENECT AND USE_LIBUSB_1 AND USE_PTHREAD)
  if(LIBFREENECT_FOUND AND LIBUSB_1_FOUND AND PTHREAD_FOUND)
    set(PER_HAVE_LIBFREENECT_AND_DEPENDENCIES TRUE)
    set(PER_HAVE_LIBFREENECT_AND_DEPENDENCIES_FOUND "yes")  # for PER-third-party.txt

    # The material is found. Check if libfreenect is an old version
    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_LIBRARIES ${LIBFREENECT_LIBRARIES} ${PTHREAD_LIBRARIES} ${LIBUSB_1_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES ${LIBFREENECT_INCLUDE_DIRS} ${PTHREAD_INCLUDE_DIRS} ${LIBUSB_1_INCLUDE_DIRS})
    CHECK_CXX_SOURCE_COMPILES("
      #include <libfreenect.hpp>

      class prMyKinect : public Freenect::FreenectDevice
      {
      };

      int main()
      {
        Freenect::Freenect<prMyKinect> freenect;
      }
      " LIBFREENECT_IS_OLD_VERSION)
    #MESSAGE("LIBFREENECT_IS_OLD_VERSION: ${LIBFREENECT_IS_OLD_VERSION}")
    if(LIBFREENECT_IS_OLD_VERSION)
      set(PER_HAVE_LIBFREENECT_OLD TRUE)  # for header prConfig.h
    else()
      set(PER_HAVE_LIBFREENECT_OLD FALSE) # for header prConfig.h
    endif()
  endif()
endif()

# check OpenCV nonfree modules and version
if(USE_OPENCV)
  set(PER_HAVE_OPENCV_VERSION "(${OpenCV_VERSION_MAJOR}*65536 + ${OpenCV_VERSION_MINOR}*256 + ${OpenCV_VERSION_PATCH})") # for prConfig.h
  if(OpenCV_VERSION)
    if(OpenCV_VERSION VERSION_LESS "2.4.0")
      message(STATUS "opencv nonfree module found")
       set(PER_HAVE_OPENCV_NONFREE TRUE)  # for header prConfig.h
    elseif(OPENCV_NONFREE_FOUND) # OpenCV < 3.0.0
      message(STATUS "opencv xfeatures2d module found")
      set(PER_HAVE_OPENCV_NONFREE TRUE)  # for header prConfig.h
    elseif(OPENCV_XFEATURES2D_FOUND) # OpenCV >= 3.0.0
      set(PER_HAVE_OPENCV_XFEATURES2D TRUE)  # for header prConfig.h
    else()
      message(STATUS "opencv nonfree or xfeature2d module not found")
    endif()
  else()
    message(STATUS "opencv nonfree not found")
    set(PER_HAVE_OPENCV_VERSION "(0)") # for prConfig.h
  endif()
endif()


# Find isnan macro (C-style)
PR_SET(PER_HAVE_FUNC_ISNAN TRUE IF HAVE_FUNC_ISNAN) # for header prConfig.h
# Find std::isnan function (cmath)
PR_SET(PER_HAVE_FUNC_STD_ISNAN TRUE IF HAVE_FUNC_STD_ISNAN) # for header prConfig.h
# Find _isnan function for MSVC
PR_SET(PER_HAVE_FUNC__ISNAN TRUE IF HAVE_FUNC__ISNAN) # for header prConfig.h
# Find isinf macro (C-style)
PR_SET(PER_HAVE_FUNC_ISINF TRUE IF HAVE_FUNC_ISINF) # for header prConfig.h
# Find std::isinf function (cmath)
PR_SET(PER_HAVE_FUNC_STD_ISINF TRUE IF HAVE_FUNC_STD_ISINF) # for header prConfig.h
# Find _finite function for MSVC
PR_SET(PER_HAVE_FUNC__FINITE TRUE IF HAVE_FUNC__FINITE) # for header prConfig.h
# Find round function (math.h)
PR_SET(PER_HAVE_FUNC_ROUND TRUE IF HAVE_FUNC_ROUND) # for header prConfig.h
# Find std::round function (cmath)
PR_SET(PER_HAVE_FUNC_STD_ROUND TRUE IF HAVE_FUNC_STD_ROUND) # for header prConfig.h

PR_SET(PER_USE_MSVC TRUE IF MSVC) # for header prConfig.h
PR_SET(PER_HAVE_CPP11_COMPATIBILITY TRUE IF USE_CPP11) # for header prConfig.h

if(USE_CPP11)
  set(PER_HAVE_CPP11_COMPATIBILITY_FOUND "yes")  # for PER-third-party.txt
else()
  set(PER_HAVE_CPP11_COMPATIBILITY_FOUND "no")  # for PER-third-party.txt
endif()

# libraries for Pioneer mobile robots
if(USE_ARIA AND UNIX)
  if(ARIA_FOUND AND USE_PTHREAD AND RT_FOUND AND DL_FOUND)
    set(PER_HAVE_PIONEER TRUE) # for header prConfig.h
    set(PER_HAVE_PIONEER_FOUND "yes")  # for PER-third-party.txt
  endif()
elseif(USE_ARIA AND NOT UNIX)
  set(PER_HAVE_PIONEER TRUE) # for header prConfig.h
  set(PER_HAVE_PIONEER_FOUND "yes")  # for PER-third-party.txt
else()
  set(PER_HAVE_PIONEER_FOUND "no")  # for PER-third-party.txt
endif()

if(USE_TX60)
  set(PER_HAVE_TX60_FOUND "yes")  # for PER-third-party.txt
else()
  set(PER_HAVE_TX60_FOUND "no")  # for PER-third-party.txt
endif()

# ----------------------------------------------------------------------------
# Finalization: generate configuration-based files
# ----------------------------------------------------------------------------

# Generate platform-dependent and configuration-dependent headers
include(cmake/PERGenerateHeaders.cmake)

# Configure the file describing how to use PER. PERConfig.cmake
# is the main file configuring a CMake package.
# . Exports build settings and dependencies for projects using PER as a
#   third party project.
# . Create and install files for simple use of find_package(PER)
#   by other cmakified "user" projects and libraries depending on PER.
#   (see "Mastering CMake", pp.72)
# . To use PER in a third party project based on CMake:
#   find_package(PER REQUIRED)
#   include_directories(${PER_INCLUDE_DIRS})
#   target_link_libraries(<target> ${PER_LIBRARIES})
include(cmake/PERGenerateConfig.cmake)


#----------------------------------------------------------------------
# Propagation in sub dirs to build demo, example, test, tutorial
#----------------------------------------------------------------------

set(PER_DIR ${PROJECT_BINARY_DIR})
mark_as_advanced(PER_DIR)
mark_as_advanced(PER_INCLUDE_DIRS)

if(BUILD_EXAMPLES)
  add_subdirectory(example)
endif()
