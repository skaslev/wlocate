set(SADA_SOURCES 386.c comparray4.c mmap.c qsufsort.c suftest3.test.c)
add_library(sada_csa STATIC ${SADA_SOURCES})
set_target_properties(sada_csa PROPERTIES LINKER_LANGUAGE C)
