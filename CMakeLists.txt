#####################################################################
#
# CMake file for HYDRA
#
#####################################################################

project(hydra)
cmake_minimum_required(VERSION 3.0)
enable_language(Fortran)

set(HYDRA_MAJOR_VERSION 0)
set(HYDRA_MINOR_VERSION 1)             # then have to give chain to minor
set(HYDRA_VERSION ${HYDRA_MAJOR_VERSION}.${HYDRA_MINOR_VERSION})


#==============================================================
# Uncomment if it is required that Fortran 90 is supported
#==============================================================

IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

#==============================================================
# Add our local modules to the module path
#==============================================================

SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
INCLUDE(${CMAKE_MODULE_PATH}/set_compiler_flags.cmake) 

MESSAGE("HYDRA MAIN: CMAKE FLAGS ARE" ${CMAKE_Fortran_FLAGS_RELEASE})

#==============================================================
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
#==============================================================

IF(${OS} STREQUAL "WINDOWS")
	set(USE_OPENMP FALSE)
ELSE()
	set(USE_OPENMP TRUE)
ENDIF()
MESSAGE("USING OPENMP?" ${USE_OPENMP})
set(USE_MPI FALSE)

set(PAR_TAG (${USE_MPI} OR ${USE_OPENMP}))
IF(${PAR_TAG})
	INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
ENDIF(${PAR_TAG})
MESSAGE("openmp flags are" ${OPENMP_Fortran_FLAGS})
#==============================================================
# There is an error in CMAKE with this flag for pgf90.  Unset it
#==============================================================

#GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
#IF(FCNAME STREQUAL "pgfortran")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
#ENDIF(FCNAME STREQUAL "pgfortran")


#==============================================================
# offer the user the choice of overriding the installation directories
#==============================================================

set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

#==============================================================
# make relative paths absolute
#==============================================================

foreach(p LIB BIN INCLUDE CMAKE)
   set(var INSTALL_${p}_DIR)
   if(NOT IS_ABSOLUTE "${${var}}")
      set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
   endif()
endforeach()

#==============================================================
# Setup the LAPACK libraries.  This also takes care of
# peculiarities, such as the fact the searching for MKL requires
# a C compiler, and that the results are not stored in the cache. 
#==============================================================

#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

#==============================================================
# Define some directories
#==============================================================

SET(SRC ${CMAKE_SOURCE_DIR}/src)
set(BLD "cmake/build")
SET(LIB "cmake/build/lib")
SET(BIN "cmake/build/bin")
set(MOD "cmake/build/modules")
set(MD2 "cmake/build/modules2")
set(PYR ${CMAKE_SOURCE_DIR}/src/Python)
set(EXD "exec")


GET_FILENAME_COMPONENT(LIB ${LIB} ABSOLUTE)
GET_FILENAME_COMPONENT(BIN ${BIN} ABSOLUTE)
GET_FILENAME_COMPONENT(MOD ${MOD} ABSOLUTE)
GET_FILENAME_COMPONENT(MD2 ${MOD} ABSOLUTE)
GET_FILENAME_COMPONENT(PYR ${PYR} ABSOLUTE)
GET_FILENAME_COMPONENT(EXD ${EXD} ABSOLUTE)
GET_FILENAME_COMPONENT(BLD ${BLD} ABSOLUTE)

#==============================================================
# Set directories
#==============================================================

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LIB})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB})
set(CMAKE_BINARY_DIR ${BIN})

MESSAGE("DETERMINED DIRECTORIES")

#==============================================================
# Set cmake policy to ignore nonsense policy: silly fellows
# Tested on CMake 3.2 on Ubuntu 14.04, 3.6 on Mac 10.10.5, 
#==============================================================

MESSAGE("SETTING CMAKE POLICY 0039 BEHAVIOR TO OLD")
cmake_policy(SET CMP0039 OLD)

#==============================================================
# The source for the FOO binary and have it placed in the bin folder
#==============================================================

# add subdirectories
add_subdirectory(src/sizing/)
add_subdirectory(src/fea/)
add_subdirectory(src/bemt/)

#####################################################################
# END OF FILE
#####################################################################
