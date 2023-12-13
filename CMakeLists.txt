cmake_minimum_required(VERSION 3.12)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(FXLT VERSION 1.0.0 LANGUAGES C)
cmake_policy(SET CMP0076 NEW)
include(CTest)

include_directories(${CMAKE_INSTALL_PREFIX}/include)
link_directories(${CMAKE_INSTALL_PREFIX}/lib)

set(lib "fused_XORier_lookup_table")
add_library(${lib} SHARED)
target_include_directories(${lib} PUBLIC .)
target_compile_options(${lib} PUBLIC -Wall -Wextra -g -Wno-implicit-fallthrough -O3  -lm)
add_subdirectory(${lib})
install(TARGETS ${lib})


set(testlibs "gnuplot_i;hash_table;linked_hash_table;rbtree")
foreach(lib IN LISTS testlibs)
add_library(${lib} STATIC)
target_include_directories(${lib} PUBLIC .)
target_compile_options(${lib} PUBLIC -Wall -Wextra -g -Wno-implicit-fallthrough -O3  -lm)
add_subdirectory(${lib})
endforeach()

add_executable(main ./main.c)

target_compile_options(main PUBLIC -O3 -Wall -Wextra -g -Wno-implicit-fallthrough -lm)
target_include_directories(main PUBLIC ${testlibs})
target_link_libraries(main ${testlibs} -lpthread -fopenmp -lfused_XORier_lookup_table)


add_custom_target(force COMMAND make clean && make)