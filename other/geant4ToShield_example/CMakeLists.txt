#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(N01)

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
#include_directories(${PROJECT_SOURCE_DIR}/include 
#                    ${Geant4_INCLUDE_DIR})
set(INCLUDE_DIRECTORIES
            ${Geant4_INCLUDE_DIR}
	    ${ROOT_INCLUDE_DIR} 
	    ${TSHIELD_ROOT_PATH}/shield/inc
	    ${TSHIELD_ROOT_PATH}/hadgen/inc
	    include 
    )
include_directories(${INCLUDE_DIRECTORIES})
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)
set(sources src/convertGeant4ToShield.cpp  src/ExN01DetectorConstruction.cpp  src/ExN01PhysicsList.cpp  src/ExN01PrimaryGeneratorAction.cpp)
#set(sources src/ExN01DetectorConstruction.cpp  src/ExN01PhysicsList.cpp  src/ExN01PrimaryGeneratorAction.cpp)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleN01 exampleN01.cpp ${sources})
#add_executable(exampleN01 exampleN01.cc ${sources})
#add_executable(exampleN01 exampleN01.cpp ${sources} ${headers})
target_link_libraries(exampleN01 ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build N01. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
#set(N01_SCRIPTS
#    exampleN01.in exampleN01.out
#  )
#
#foreach(_script ${N01_SCRIPTS})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/${_script}
#    ${PROJECT_BINARY_DIR}/${_script}
#    COPYONLY
#    )
#endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(N01 DEPENDS exampleN01)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS exampleN01 DESTINATION bin)

