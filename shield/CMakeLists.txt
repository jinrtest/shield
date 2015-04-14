enable_language (C Fortran)

include_directories(inc ${TSHIELD_ROOT_PATH}/hadgen/inc)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -w")

set(SHIELD_SRCS
    src/shield_decode_tree.c 
    src/shield_geo.c 
    src/shield_init.c 
    src/shield_iso.c 
    src/shield_media.c 
    src/shield_setget.c 
     
    src/Fragments.f
    src/Inspar.f
    src/Loent28.f
    src/Medium.f
    src/Microd.f
    src/Shield_so.f
    src/Tree.f
    src/decay1.f
    src/elscat.f
    src/gemca.f
    src/files.f
    src/gentree.f
    src/lndata28.f
    src/range.f
    src/Output.f
)

add_library(SHIELD SHARED ${SHIELD_SRCS})
target_link_libraries(SHIELD HADGEN)

