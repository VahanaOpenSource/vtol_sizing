#==============================================================
# Preproecessing for f90wrap: set project name
#==============================================================

set(prj bemt)

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
# add source subdirectory
#==============================================================

ADD_SUBDIRECTORY(${SRC}) 
 
