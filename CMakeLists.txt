# ==============================================================================
# LBFGS CMake File
# ==============================================================================

CMAKE_MINIMUM_REQUIRED(VERSION 2.5 FATAL_ERROR)
#SET(CMAKE_VERBOSE_MAKEFILE ON)

# ==============================================================================
# conditional setup ------------------------------------------------------------
# ==============================================================================

SET(LIBS_SHARED       ON  CACHE BOOL "Should the shared version of the library be built?")
SET(LIBS_STATIC       OFF CACHE BOOL "Should the static version of the library be built?")

# ==============================================================================
# project setup ----------------------------------------------------------------
# ==============================================================================

# define project name ----------------------------------------------------------
PROJECT(LBFGS)
SET(LBFGS_VERSION 1.0)

# what languages we use --------------------------------------------------------
ENABLE_LANGUAGE(Fortran)

# ==============================================================================
# external dependencies --------------------------------------------------------
# ==============================================================================

IF(NOT DEFINED COMPOSITE_PROJECT)
    SET(DEVELOPMENT_ROOT "$ENV{DEVELOPMENT_ROOT}")
ENDIF(NOT DEFINED COMPOSITE_PROJECT)

IF(NOT DEVELOPMENT_ROOT)
    MESSAGE(FATAL_ERROR "DEVELOPMENT_ROOT variable is not set!")
ENDIF(NOT DEVELOPMENT_ROOT)

# ==============================================================================
# external libraries -----------------------------------------------------------
# ==============================================================================

IF(NOT MKL_FOUND)
    # BLAS ======================
    SET(BLAS_ROOT ${DEVELOPMENT_ROOT}/thirdparty/blas/netlib)
    LINK_DIRECTORIES(${BLAS_ROOT}/lib)
    SET(BLAS_LIB_NAME fblasd)
ELSE(NOT MKL_FOUND)
    SET(LAPACK_LIB_NAME ${MKL_LIBS})
ENDIF(NOT MKL_FOUND)

# ==============================================================================
# output directories -----------------------------------------------------------
# ==============================================================================

SET(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)

# ==============================================================================
# rpaths -----------------------------------------------------------------------
# ==============================================================================

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 

# the RPATH to be used when installing
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)

# ==============================================================================
# compiler setups --------------------------------------------------------------
# ==============================================================================

# get the name of fortran compiler
get_filename_component(FCOMP_NAME ${CMAKE_Fortran_COMPILER} NAME_WE)

# default compiler setups - gfortran
SET(CMAKE_Fortran_FLAGS_RELEASE "-O0")
SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0")

# ifort
IF(FCOMP_NAME STREQUAL ifort)
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    SET(CMAKE_Fortran_FLAGS_DEBUG   "-g -traceback -check all -O0")
ENDIF(FCOMP_NAME STREQUAL ifort)

# fix the problem with dynamic libraries and xlf compiler on SP5
IF(FCOMP_NAME STREQUAL xlf)
    SET(CMAKE_Fortran_FLAGS_RELEASE "-O2 -qsuppress=1501-510:cmpmsg")
    SET(CMAKE_Fortran_FLAGS_DEBUG   "-g")
    IF(LIBS_SHARED)
        MESSAGE(FATAL_ERROR "LIBS_SHARED is ON, which is not supported with xlf compiler!")
    ENDIF(LIBS_SHARED)
ENDIF(FCOMP_NAME STREQUAL xlf)

SET(CMAKE_BUILD_TYPE Release)
#SET(CMAKE_BUILD_TYPE Debug)

# ==============================================================================
# project subdirectories  ------------------------------------------------------
# ==============================================================================

ADD_SUBDIRECTORY(src)
