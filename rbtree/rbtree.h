#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct rbtree;

void rbtree_print(struct rbtree * mapThis, FILE * out,void (*printKey)(void * thing,FILE * out),void (*printVal)(void * thing,FILE * out));
int rbtree_getKeys(struct rbtree *mapThis, void ***keys);
void rbtree_clear(struct rbtree *mapThis);
struct rbtree *rbtree_copy(struct rbtree *mapThis);
bool rbtree_contains(struct rbtree *mapThis, void * key);
void* rbtree_getValue(struct rbtree *mapThis, void * key);
bool rbtree_tryGetValue(struct rbtree *mapThis, void * key, void **out);
bool rbtree_setValue(struct rbtree *mapThis, void * key, void * value);
bool rbtree_insert(struct rbtree *mapThis, void * key, void * value);
bool rbtree_remove(struct rbtree *mapThis, void * key);
void destroy_rbtree(struct rbtree *mapThis);
struct rbtree *construct_rbtree();
size_t rbtree_getSize(struct rbtree *mapThis);