#include "./rbtree.h"
#include "./rbnode/rbnode.h"

struct rbnode;

struct rbtree {
    struct rbnode *root; // = NULL;
};


void rbtreePrintBookends(struct rbtree *this, FILE * out) {
    fprintf(out, "---------bh:%d rp: %s---------\n", (this->root ? this->root->vtable->bh(this->root) : 0), (this->root ? (this->root->vtable->rp(this->root) ? "✅" : "❌ ") : "✅"));
}

int rbtree_checkValid(struct rbtree *this) {
    if (!this->root) {
        return true;
    }
    return this->root->vtable->rp(this->root) && (this->root->vtable->bh(this->root) != -1);
}

void rbtree_print(struct rbtree * mapThis, FILE * out,void (*printKey)(void * thing,FILE * out),void (*printVal)(void * thing,FILE * out)) {
    struct rbtree * this = (struct rbtree *)mapThis;
    rbtreePrintBookends(this, out);
    if (!this->root) {
        return;
    }
    this->root->vtable->printTree(this->root, 0, out, printKey, printVal);
    rbtreePrintBookends(this, out);
}

int rbtree_getKeys(struct rbtree *mapThis, void ***keys) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        return 0;
    }
    return this->root->vtable->getKeys(this->root, keys);
}

void rbtree_clear(struct rbtree *mapThis) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        return;
    }
    destroy_rbnode(this->root);
    free(this->root);
    this->root = NULL;
}

struct rbtree *rbtree_copy(struct rbtree *mapThis) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        return construct_rbtree();
    }
    struct rbtree *ret = (struct rbtree *)construct_rbtree();
    ret->root = this->root->vtable->copy(this->root, NULL);
    return (struct rbtree *)ret;
}

bool rbtree_contains(struct rbtree *mapThis, void * key) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        return false;
    }
    return this->root->vtable->contains(this->root, key);
};

void* rbtree_getValue(struct rbtree *mapThis, void * key) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        return 0;
    }
    return this->root->vtable->getValue(this->root, key);
};

bool rbtree_tryGetValue(struct rbtree *mapThis, void * key, void **out) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        return false;
    }
    return this->root->vtable->tryGetValue(this->root, key, out);
};

bool rbtree_setValue(struct rbtree *mapThis, void * key, void * value) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        return false;
    }
    return this->root->vtable->setValue(this->root, key, value);
}

bool rbtree_insert(struct rbtree *mapThis, void * key, void * value) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        this->root = construct_rbnode(key, value, this->root, &this->root);
        return true;
    }
    int ret = this->root->vtable->insert(this->root, key, value);
    assert(rbtree_checkValid(this));
    return ret;
};

bool rbtree_remove(struct rbtree *mapThis, void * key) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        return false;
    }
    int ret = this->root->vtable->remove(this->root, key, &this->root);
    assert(rbtree_checkValid(this));
    return ret;
}

void destroy_rbtree(struct rbtree *mapThis) {
    rbtree_clear(mapThis);
}
struct rbtree *construct_rbtree() {
    struct rbtree *new = calloc(1, sizeof *new);
    new->root = NULL;

    return (struct rbtree *)new;
}

size_t rbtree_getSize(struct rbtree *mapThis) {
    struct rbtree * this = (struct rbtree *)mapThis;
    if (!this->root) {
        return 0;
    }
    void ** temp;
    size_t nodes = this->root->vtable->getKeys(this->root, &temp);
    free(temp);
    return nodes * sizeof(struct rbnode) + sizeof(struct rbtree);
}