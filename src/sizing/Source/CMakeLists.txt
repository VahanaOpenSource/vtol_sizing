
#==============================================================
# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
#==============================================================

#==============================================================
# Create empty list
#==============================================================

SET(lib_src "")

#==============================================================
# Add the source files from each directory
#==============================================================

#==============================================================
# Python routines
#==============================================================

file(GLOB pyint_files_src pyint_files/*.F90)

#==============================================================
# component mass/etc routines
#==============================================================

file(GLOB main_src pyinv_files/*.F90)

#==============================================================
# Modules visible to python
#==============================================================

file(GLOB pyint_modules_src pyint_modules/*.F90)

#==============================================================
# add files to list of sources
#==============================================================

LIST(APPEND all_src ${pyint_files_src})
LIST(APPEND all_src ${pyint_modules_src})
LIST(APPEND all_src ${main_src})

#==============================================================
# add files to list of sources to create library
#==============================================================

LIST(APPEND lib_src ${pyint_files_src})
LIST(APPEND lib_src ${pyint_modules_src})
LIST(APPEND lib_src ${main_src})

#==============================================================
# Define the executable in terms of the source files
#==============================================================

#MESSAGE("Before auto-removing O3 flag" ${CMAKE_Fortran_FLAGS})
foreach(flags CMAKE_Fortran_FLAGS_RELEASE)
   string(REPLACE "-O3" "" ${flags} "${${flags}}")    # remove -O3
#   string(REPLACE "-Mipa=fast" "" ${flags} "${${flags}}")    # remove -Mipa=fast
endforeach()
#MESSAGE(STATUS ${CMAKE_Fortran_FLAGS_RELEASE})
#MESSAGE(STATUS ${main_src})
#ADD_EXECUTABLE(${FOOEXE} ${all_src})

#==============================================================
# include preprocessor function from include file
#==============================================================

INCLUDE(${CMAKE_MODULE_PATH}/f90wrap.cmake) 

#==============================================================
# 1. list of source files that need to be wrapped (APIs)
#==============================================================

set(src_list ${pyint_files_src} ${pyint_modules_src})

#==============================================================
# 2. Create user-defined libraries and list libraries & locations
#==============================================================

set(UDL fortran_routines_${prj})
add_library(${UDL} ${lib_src})
MESSAGE("OPENMP FORTRAN FLAGS ARE" ${OpenMP_Fortran_FLAGS})
IF(USE_OPENMP)
  SET_TARGET_PROPERTIES(${UDL} PROPERTIES COMPILE_FLAGS ${OpenMP_Fortran_FLAGS})
ENDIF()
set(UDP ${LIB})

#==============================================================
# 3. Run f90wrap
#==============================================================

set(FC ${CMAKE_Fortran_COMPILER})
MESSAGE("LIBRARIES ARE " ${UDL})
run_f90wrap("${FC}" "${src_list}" "${prj}" "${MOD}" "${UDL}" "${UDP}" "${PYD}")

#==============================================================
# Uncomment if you have parallization
#==============================================================

MESSAGE("OPENMP Fortran Flag used is " ${OpenMP_Fortran_FLAGS})
IF(USE_OPENMP)
    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
ELSEIF(USE_MPI)
    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
    TARGET_LINK_LIBRARIES(${FOOEXE} ${MPI_Fortran_LIBRARIES})
ENDIF(USE_OPENMP)

#==============================================================
# Tell how to install this executable
#==============================================================

IF(WIN32)
    SET(CMAKE_INSTALL_PREFIX "C:\\Program Files")
ELSE()
    SET(CMAKE_INSTALL_PREFIX /usr/local)
ENDIF(WIN32)

