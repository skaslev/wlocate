cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(wlocate VERSION 1.0 LANGUAGES C)

set(WLOCATE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/wlocate)

set(CMAKE_C_FLAGS "-std=c99 -Wall -pedantic-errors")
add_definitions(-D_XOPEN_SOURCE=700)

add_subdirectory(sada_csa)

add_executable(wupdatedb wupdatedb.c buf.c die.c)
target_link_libraries(wupdatedb sada_csa)
set_target_properties(wupdatedb PROPERTIES LINKER_LANGUAGE C)
install(TARGETS wupdatedb DESTINATION ${WLOCATE_INSTALL_DIR})

add_executable(wlocate wlocate.c buf.c die.c)
target_link_libraries(wlocate sada_csa)
set_target_properties(wlocate PROPERTIES LINKER_LANGUAGE C)
install(TARGETS wlocate DESTINATION ${WLOCATE_INSTALL_DIR})
