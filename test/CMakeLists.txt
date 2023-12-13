set(testlibs "gnuplot_i;hash_table;linked_hash_table;rbtree")
foreach(lib IN LISTS testlibs)
add_library(${lib} STATIC)
target_include_directories(${lib} PUBLIC .)
target_compile_options(${lib} PUBLIC -Wall -Wextra -g -Wno-implicit-fallthrough -O3  -lm)
add_subdirectory(${lib})
endforeach()

add_executable(main ./main.c)

target_compile_options(main PUBLIC -O3 -Wall -Wextra -g -Wno-implicit-fallthrough -lm -fopenmp)
target_include_directories(main PUBLIC ${testlibs})
target_link_libraries(main ${testlibs} -lpthread -fopenmp -lfuse_XORier_lookup_table)