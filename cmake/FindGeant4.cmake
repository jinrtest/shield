# - Try to find Geant4 instalation
# This module sets up Geant4 information: 
# - either from Geant4 CMake configuration file (Geant4Config.cmake), if available
# - or it defines:
# GEANT4_FOUND          If Geant4 is found
# GEANT4_INCLUDE_DIR    PATH to the include directory
# GEANT4_LIBRARIES      Most common libraries
# GEANT4_LIBRARY_DIR    PATH to the library directory 

message(STATUS "Looking for GEANT4 ...")

find_path(Geant4_DIR NAMES geant4-config PATHS
  ${SIMPATH}/transport/geant4/bin
  ${SIMPATH}/bin
  NO_DEFAULT_PATH
)

if(Geant4_DIR)
  Set(PATH ${PATH} ${Geant4_DIR})
endif(Geant4_DIR)

# First search for Geant4Config.cmake on the path defined via user setting 
# Geant4_DIR

if(EXISTS ${Geant4_DIR}/Geant4Config.cmake)
  include(${Geant4_DIR}/Geant4Config.cmake)
  set(GEANT4_INCLUDE_DIR ${Geant4_INCLUDE_DIRS})
  message(STATUS "Looking for GEANT4... - Found Geant4 CMake configuration in ${Geant4_DIR}")
  return()
endif()

# If Geant4_DIR was not set search for geant4-config executable on system path
# to get Geant4 installation directory 

find_program(GEANT4_CONFIG_EXECUTABLE geant4-config PATHS
  ${PATH}
  )

if(GEANT4_CONFIG_EXECUTABLE)
  execute_process(
    COMMAND ${GEANT4_CONFIG_EXECUTABLE} --prefix 
    OUTPUT_VARIABLE G4PREFIX 
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(
    COMMAND ${GEANT4_CONFIG_EXECUTABLE} --version 
    OUTPUT_VARIABLE GEANT4_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(EXISTS ${G4PREFIX}/lib/Geant4-${GEANT4_VERSION}/Geant4Config.cmake)
    set(Geant4_DIR ${G4PREFIX}/lib/Geant4-${GEANT4_VERSION})
    include(${Geant4_DIR}/Geant4Config.cmake)
    set(GEANT4_INCLUDE_DIR ${Geant4_INCLUDE_DIRS})
    message(STATUS "Looking for GEANT4... - Found Geant4 CMake configuration in ${Geant4_DIR}")
    return()
  endif()

  if(EXISTS ${G4PREFIX}/lib64/Geant4-${GEANT4_VERSION}/Geant4Config.cmake)
    set(Geant4_DIR ${G4PREFIX}/lib64/Geant4-${GEANT4_VERSION})
    include(${Geant4_DIR}/Geant4Config.cmake)
    set(GEANT4_INCLUDE_DIR ${Geant4_INCLUDE_DIRS})
    message(STATUS "Looking for GEANT4... - Found Geant4 CMake configuration in ${Geant4_DIR}")
    return()
  endif()

endif()

# If search for geant4-config failed try to use directly user paths if set
# or environment variables 
# We can add Geant4Config.cmake  and UseGeant4.cmake to local cmake directory and include them. 
# Now liblist and libname not exist in Geant4 at FairRoot.
#if (NOT GEANT4_FOUND)
#  find_path(GEANT4_INCLUDE_DIR NAMES G4Event.hh PATHS
#    ${SIMPATH}/transport/geant4/include
#    ${SIMPATH}/transport/geant4/include/geant4
#    ${SIMPATH}/transport/geant4/include/Geant4
#    ${SIMPATH}/include/geant4
#    ${SIMPATH}/include/Geant4
#    NO_DEFAULT_PATH
#  )
#  set(GEANT4_INCLUDE_DIR
#    ${GEANT4_INCLUDE_DIR}
#    ${SIMPATH}/transport/geant4/source/interfaces/common/include 
#    ${SIMPATH}/transport/geant4/physics_lists/hadronic/Packaging/include   
#    ${SIMPATH}/transport/geant4/physics_lists/hadronic/QGSP/include
#  )
#  find_path(GEANT4_LIBRARY_DIR NAMES libG3toG4.so PATHS
#    ${SIMPATH}/transport/geant4/lib/Linux-g++
#    ${SIMPATH}/transport/geant4/lib/Linux-icc
#    ${SIMPATH}/transport/geant4/lib
#    ${SIMPATH}/lib
#    NO_DEFAULT_PATH
#  )
#  if (GEANT4_INCLUDE_DIR AND GEANT4_LIBRARY_DIR)
#    execute_process(
#      COMMAND ${GEANT4_LIBRARY_DIR}/liblist -m ${GEANT4_LIBRARY_DIR}                  
#      INPUT_FILE ${GEANT4_LIBRARY_DIR}/libname.map 
#      OUTPUT_VARIABLE GEANT4_LIBRARIES
#      OUTPUT_STRIP_TRAILING_WHITESPACE
#      TIMEOUT 2)
#  endif()
#  set(GEANT4_LIBRARIES "-L${GEANT4_LIBRARY_DIR} ${GEANT4_LIBRARIES} -lexpat -lz")
#endif()

if (GEANT4_INCLUDE_DIR AND GEANT4_LIBRARY_DIR AND GEANT4_LIBRARIES)
  set (GEANT4_FOUND TRUE)
endif()

if (GEANT4_FOUND)
  if (NOT GEANT4_FIND_QUIETLY)
    MESSAGE(STATUS "Looking for GEANT4... - found ${GEANT4_LIBRARY_DIR}")
#    message(STATUS "Found ${GEANT4_LIBRARY_DIR}")
  endif (NOT GEANT4_FIND_QUIETLY)
  SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${GEANT4_LIBRARY_DIR})
else (GEANT4_FOUND)
  if (GEANT4_FIND_REQUIRED)
    message(FATAL_ERROR "Looking for GEANT4... - Not found")
  endif (GEANT4_FIND_REQUIRED)
endif (GEANT4_FOUND)

# Make variables changeble to the advanced user
mark_as_advanced(GEANT4_INCLUDE_DIR GEANT4_LIBRARY_DIR GEANT4_LIBRARIES)


