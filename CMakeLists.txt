cmake_minimum_required (VERSION 2.8 FATAL_ERROR)
project (TShield)

set(TSHIELD_ROOT_PATH ${CMAKE_CURRENT_SOURCE_DIR})
include(cmake/TShieldConfig.cmake)

add_subdirectory (hadgen)
add_subdirectory (shield)
add_subdirectory (THadgen)
add_subdirectory (TShield)
