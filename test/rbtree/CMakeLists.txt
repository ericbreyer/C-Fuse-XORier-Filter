target_include_directories(rbtree PUBLIC . ./rbnode)
target_sources(rbtree PRIVATE rbtree.c rbtree.h ./rbnode/rbnode.c ./rbnode/rbnode.h)