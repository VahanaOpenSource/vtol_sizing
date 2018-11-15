#==============================================================
# Determine and set the Fortran compiler flags we want 
#==============================================================

#==============================================================
# Make sure that the default build type is RELEASE if not specified.
#==============================================================

INCLUDE(${CMAKE_MODULE_PATH}/SetCompileFlag.cmake)

#==============================================================
# Make sure the build type is uppercase
#==============================================================

STRING(TOUPPER "${CMAKE_BUILD_TYPE}" BT)

IF(BT STREQUAL "RELEASE")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE BUILD TYPE RELEASE SELECTED")
ELSEIF(BT STREQUAL "DEBUG")
    SET (CMAKE_BUILD_TYPE DEBUG CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE BUILD TYE DEBUG SELECTED")
ELSEIF(BT STREQUAL "TESTING")
    SET (CMAKE_BUILD_TYPE TESTING CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE BUILD TYPE TESTING SELECTED")
ELSEIF(NOT BT)
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are DEBUG, RELEASE, or TESTING."
      FORCE)
    MESSAGE(STATUS "CMAKE_BUILD_TYPE not given, defaulting to RELEASE")
ELSE()
    MESSAGE(FATAL_ERROR "CMAKE_BUILD_TYPE not valid, choices are DEBUG, RELEASE, or TESTING")
ENDIF(BT STREQUAL "RELEASE")

#==============================================================
# If the compiler flags have already been set, return now
#==============================================================

IF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)
    RETURN ()
ENDIF(CMAKE_Fortran_FLAGS_RELEASE AND CMAKE_Fortran_FLAGS_TESTING AND CMAKE_Fortran_FLAGS_DEBUG)

#==============================================================
# Determine the appropriate flags for this compiler for each build type.
# For each option type, a list of possible flags is given that work
# for various compilers.  The first flag that works is chosen.
# If none of the flags work, nothing is added (unless the REQUIRED 
# flag is given in the call).  This way unknown compiles are supported.
#==============================================================

#==============================================================
### GENERAL FLAGS ###
#==============================================================

#==============================================================
# Don't add underscores in symbols for C-compatability
#==============================================================

#SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
#                 Fortran "-fno-underscoring")

#==============================================================
# There is some bug where -march=native doesn't work on Mac
#==============================================================

IF(APPLE)
    SET(GNUNATIVE "-mtune=native")
ELSE()
    SET(GNUNATIVE "-march=native")
ENDIF()

#==============================================================
# Define the operating system
#==============================================================

SET(OS ${CMAKE_SYSTEM_NAME})
SET(FC ${CMAKE_Fortran_COMPILER})

STRING(TOUPPER "${OS}" OS)
STRING(TOUPPER "${FC}" FC)

SET(Wintel FALSE)
IF(${OS} STREQUAL "WINDOWS")
   IF(${FC} MATCHES "INTEL")
      SET(Wintel TRUE)
   ENDIF()
ENDIF()

MESSAGE("The Operating System Type is " ${OS})
MESSAGE("The Fortran Compiler is      " ${FC})

#==============================================================
# add some definitions
#==============================================================

add_definitions("-fPIC ") #-fbounds=check # for gfortran debug
#add_definitions("-Ddisplay_off")
#add_definitions("-Dsmall_angles")
#add_definitions("-Dperformance_only")

#==============================================================
# Optimize for the host's architecture
#==============================================================

IF(${Wintel}) 
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}"
                 Fortran "/QxHost"       # Intel Windows
                )
ELSEIF(${FC} MATCHES "INTEL" OR ${FC} MATCHES "GFORTRAN")
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}"
                 Fortran "-xHost"        # Intel
                         ${GNUNATIVE}    # GNU
                )
ENDIF()

#==============================================================
# Set GPU suppression flags for non-PGI compilers
#==============================================================

IF(${FC} MATCHES "PGFORTRAN")

  IF(USE_OPENMP)
    add_definitions("-Dopenmp")
  ENDIF()

  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}"
                 Fortran "-Mcuda=cc50 -Minfo=accel -ta=tesla:nordc"       #compute capability 6
                         "-Mcuda=cc60 -Minfo=accel -ta=tesla:nordc"        #compute capability 5x
                         "-Mcuda=cc5x -Minfo=accel -ta=tesla:nordc"        #compute capability 4x
                         "-Mcuda=cc4x -Minfo=accel -ta=tesla:nordc"        #compute capability 3x
                )        # if you have something else, not supported

#  MESSAGE("right after" ${CMAKE_Fortran_FLAGS_RELEASE})
#  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}"
#                 Fortran "-Minfo=accel"        # PGI
#                )

#  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS}"
#                 Fortran "-ta=nvidia"        # PGI/Nvidia
#                )

#  MESSAGE("right after" ${CMAKE_Fortran_FLAGS_RELEASE})

#==============================================================
# Non-PGI compilers: no GPU support
#==============================================================

ELSE()

#==============================================================
# add defnition to suppress GPU part of code
#==============================================================

  add_definitions(-Dnogpu)

#==============================================================
# GNNU compilers: add preprocessing directives
#==============================================================

  IF(${FC} MATCHES "GFORTRAN" OR ${FC} MATCHES "F95" OR ${FC} MATCHES "INTEL")
   
#==============================================================
# Remove parantheses protection
#==============================================================

  SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-fno-protect-parens" # GNU
                )
  ENDIF()
#==============================================================
# Stack arrays
#==============================================================

  IF(${FC} MATCHES "GFORTRAN")
    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-fstack-arrays" # GNU
                )

    IF(USE_OPENMP)
      add_definitions("-cpp -Dopenmp")
    ENDIF()
    
#==============================================================
# Intel fortran compiler
#==============================================================

  ELSEIF(${FC} MATCHES "INTEL")
    IF(USE_OPENMP)
      add_definitions("-fpp -Dopenmp -qopenmp")  
    ENDIF()
  ELSE()
    MESSAGE(STATUS "Unknown compiler: ignoring OpenMP preprocess filter")
  ENDIF()

ENDIF()
  
  
#==============================================================
### DEBUG FLAGS ###
#==============================================================

# NOTE: debugging symbols (-g or /debug:full) are already on by default

#==============================================================
# Disable optimizations
#==============================================================

# IF(${Wintel}) 
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#                  Fortran REQUIRED "/Od" # Intel Windows
#                 )
# ELSE()
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#                  Fortran REQUIRED "-O0" # All compilers not on Windows
#                 )
# ENDIF()

#==============================================================
# Turn on all warnings 
#==============================================================

# IF(OS STREQUAL "WINDOWS") 
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#                  Fortran "/warn:all" # Intel Windows
#                 )
# ELSE()
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#                  Fortran "-warn all" # Intel
#                          "-Wall"     # GNU
#                                      # Portland Group (on by default)
#                 )
# ENDIF()

#==============================================================
# Traceback
#==============================================================

# IF(${Wintel}) 
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#                  Fortran "/traceback"   # Intel Windows
#                 )
# ELSE()
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#                  Fortran "-traceback"   # Intel/Portland Group
#                          "-fbacktrace"  # GNU (gfortran)
#                          "-ftrace=full" # GNU (g95)
#                 )
# ENDIF()

#==============================================================
# Check array bounds
#==============================================================

# IF(${Wintel}) 
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#                  Fortran "/check:bounds"  # Intel Windows
#                 )
# ELSE()
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
#                  Fortran "-check bounds"  # Intel
#                          "-fcheck=bounds" # GNU (New style)
#                          "-fbounds-check" # GNU (Old style)
#                          "-Mbounds"       # Portland Group
#                 )
# ENDIF()

#==============================================================
### TESTING FLAGS ###
#==============================================================

#==============================================================
# Optimizations
#==============================================================

# IF(${Wintel}) 
# SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}"
#                  Fortran REQUIRED "/O2" # Intel Windows
#                 )
# ELSE()
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_TESTING "${CMAKE_Fortran_FLAGS_TESTING}"
#                  Fortran REQUIRED "-O2" # All compilers not on Windows
#                 )
# ENDIF()

#==============================================================
### RELEASE FLAGS ###
#==============================================================

# NOTE: agressive optimizations (-O3) are already turned on by default

#==============================================================
# Unroll loops
#==============================================================

# IF(${Wintel}) 
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                  Fortran "/unroll"        # Intel Windows
#                 )
# ELSE()
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                  Fortran "-funroll-loops" # GNU
#                          "-unroll"        # Intel
#                          "-Munroll"       # Portland Group
#                 )
# ENDIF()

#==============================================================
# Fast flag
#==============================================================

IF(${Wintel}) 
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "/fast"        # Intel Windows
                )
ELSEIF(${FC} MATCHES "GFORTRAN" OR ${FC} MATCHES "INTEL")
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-Ofast" # GNU
                         "-fast"    # Intel and PGI; PGI has some problem..
                )
ENDIF()

#==============================================================
# Inline functions
#==============================================================

# IF(${Wintel}) 
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                  Fortran "/Qinline"           # Intel Windows
#                 )
# ELSE()
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                  Fortran "-inline"            # Intel
#                          "-finline-functions" # GNU
#                          "-Minline"           # Portland Group
#                 )
# ENDIF()

#==============================================================
# Interprocedural (link-time) optimizations
#==============================================================

# IF(${Wintel})
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                  Fortran "/Qipo"    # Intel Windows
#                 )
# ELSE()
#    SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                  Fortran "-ipo"     # Intel
#                          "-flto"    # GNU
#                          "-Mipa"   # Portland Group
#                 )
# ENDIF()

#==============================================================
# Single-file optimizations
#==============================================================

IF(${Wintel})
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "/Qip" # Intel Windows
                )
ELSE()
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-ip"  # Intel
                )
ENDIF()

#==============================================================
# Extend source
#==============================================================

IF(${Wintel})
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "/extend-source" # Intel Windows
                )
ELSE()
   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
                 Fortran "-Mextend"               # PGI
                         "-ffixed-line-length-132"# gnu
                         "-extend_source -132"    # Intel    
                )
ENDIF()

#==============================================================
# Turn off loop vectorization diagnostics; parallelization done by AS
#==============================================================

#IF(${Wintel})
#   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                 Fortran "/Qvec-report0" # Intel Windows
#                )
#ELSE()
#   SET_COMPILE_FLAG(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
#                 Fortran "-vec-report0"  # Intel
#                         "-Mvect"        # Portland Group
#                )
#ENDIF()

