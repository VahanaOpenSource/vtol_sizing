#==============================================================
# Preproecessing for f90wrap: set project name
#==============================================================

set(prj fea)

#==============================================================
# Define some directories
#==============================================================

SET(SRC ${CMAKE_CURRENT_SOURCE_DIR}/Source)
set(PYD ${PYR}/Stage_3)

GET_FILENAME_COMPONENT(PYD ${PYD} ABSOLUTE)

#==============================================================
# Have the .mod files placed in a separate folder
#==============================================================

SET(MOD ${BIN}/modules_${prj})

SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})

#==============================================================
# link the mathops library
#==============================================================

#LINK_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/utils/mathops)
#
#TARGET_LINK_LIBRARIES(${prj} mathops)

#==============================================================
# The source for the FOO binary and have it placed in the bin folder
#==============================================================

ADD_SUBDIRECTORY(${SRC}) 
 
