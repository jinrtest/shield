enable_language (C Fortran)

add_subdirectory (parse)
include_directories(inc ${CMAKE_CURRENT_BINARY_DIR}/parse)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")

set(HADGEN_SRCS
    src/hadgen_common_blocks.c
    src/hadgen_main.c
    src/hadgen_options.c
    src/hadgen_prev.c
    src/hadgen_tabnuc.c
    src/hadgen_time.c
    src/hadgen_setdky.c

    src/HistAGT.f
    src/Qgsmn1.f
    src/Qgsmn2.f
    src/casevp.f
    src/coales.f
    src/deex.f
    src/files.f
    src/genagt.f
    src/gqaar.f
    src/gqamn.f
    src/gqstn.f
    src/hadgen_so.f
    src/preco.f
    src/ranlux.f
)

add_library(HADGEN SHARED ${HADGEN_SRCS})
add_dependencies(HADGEN parse_run)

