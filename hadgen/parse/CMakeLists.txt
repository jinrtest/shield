enable_language (C Fortran)
set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR})
link_directories(${LINK_DIRECTORIES})

set(HADGEN_parse_SRCS
    src/fparse.f
    src/main.c
)
add_executable(parse ${HADGEN_parse_SRCS})

#Copy files for run generated executable.
configure_file(TABNUC.DAT TABNUC.DAT COPYONLY)
configure_file(ATAB.DAT   ATAB.DAT   COPYONLY)
#Run executable after build (always outdated target)
add_custom_target(parse_run DEPENDS parse TABNUC.DAT ATAB.DAT VERBATIM)
add_custom_command(
    TARGET parse_run
    POST_BUILD
    COMMAND parse TABNUC.DAT
    DEPENDS parse TABNUC.DAT ATAB.DAT
    COMMENT "Parsing TABNUC.DAT and ATAB.DAT for HADGEN"
    VERBATIM
)
