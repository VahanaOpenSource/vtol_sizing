#==============================================================
# Preproecessing for f90wrap: set project name
#==============================================================

set(prj hydra)

#==============================================================
# Define some directories
#==============================================================

SET(SRC ${CMAKE_CURRENT_SOURCE_DIR}/Source)
set(PYD ${PYR}/Stage_1)

GET_FILENAME_COMPONENT(PYD ${PYD} ABSOLUTE)

#==============================================================
# Have the .mod files placed in a separate folder
#==============================================================

SET(MOD ${BIN}/modules_${prj})

SET(CMAKE_Fortran_MODULE_DIRECTORY ${MOD})

#==============================================================
# The source for the FOO binary and have it placed in the bin folder
#==============================================================

ADD_SUBDIRECTORY(${SRC}) 
 
