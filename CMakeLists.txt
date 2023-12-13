cmake_minimum_required(VERSION 3.12)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(FXLT VERSION 1.0.0 LANGUAGES C)
cmake_policy(SET CMP0076 NEW)
include(CTest)

include_directories(${CMAKE_INSTALL_PREFIX}/include)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

set(lib "fuse_XORier_lookup_table")
add_library(${lib} SHARED)
target_compile_options(${lib} PUBLIC -Wall -Wextra -g -Wno-implicit-fallthrough -O3  -lm)
add_subdirectory(src)
install(TARGETS ${lib})

add_custom_target(fxlt_install ALL make fuse_XORier_lookup_table && sudo cmake --install .)
add_custom_target(force COMMAND make clean && make)

add_subdirectory(test)