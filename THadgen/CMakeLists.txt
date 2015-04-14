enable_language (C CXX Fortran)

set(INCLUDE_DIRECTORIES ${ROOT_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/inc ${TSHIELD_ROOT_PATH}/hadgen/inc)
include_directories(${INCLUDE_DIRECTORIES})
link_directories(${ROOT_LIBRARY_DIR})

configure_file(inc/THadgen.h inc/THadgen.h COPYONLY)
configure_file(inc/LinkDef.h inc/LinkDef.h COPYONLY)
set(G__Hadgen_HEADERS inc/THadgen.h)
set(G__Hadgen_LINKDEF inc/LinkDef.h)
set(G__Hadgen_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/G__Hadgen.cxx)
ROOT_GENERATE_DICTIONARY("${G__Hadgen_HEADERS}" "${G__Hadgen_LINKDEF}" "${G__Hadgen_DICTIONARY}" "${INCLUDE_DIRECTORIES}")

set(THadgen_SRCS
    src/THadgen.cxx
    src/THadgenRoutines.cxx
    ${G__Hadgen_DICTIONARY}
)

add_library(THadgen SHARED ${THadgen_SRCS})
target_link_libraries(THadgen HADGEN ${ROOT_LIBRARIES})

