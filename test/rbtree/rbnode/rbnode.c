#include "./rbnode.h"


void rbnode_static_doubleBlackFixup(struct rbnode *dbnode, struct rbnode *parent);

int rbnode_compare(struct rbnode *this, void * other) {
    return this->key - other;
}

void rbnode_commitApoptosis(struct rbnode *this, struct rbnode **parentsChild) {
    // if there are two children (hardest case)
    //  make this node's key the key of it's successor, then delete the successor node
    if (this->l && this->r) {
        struct rbnode *prev = this;
        struct rbnode *cur = this->r;
        while (cur->l) {
            prev = cur;
            cur = cur->l;
        }
        this->key = cur->key;
        if (prev->r == cur) {
            // prev->r = NULL;
            cur->vtable->commitApoptosis(cur, &(prev->r));
        } else {
            // prev->l = NULL;
            cur->vtable->commitApoptosis(cur, &(prev->l));
        }
        return;
    }
    // if there is a left or right child
    // replace this node with it, and fixup if need be
    else if (this->l) {
        *parentsChild = this->l;
        this->l->p = this->p;
        if (this->c == red) {
            // nothing, all good
        } else if (this->l->c == red) {
            this->l->c = black;
        } else if (this->c == black) {
            rbnode_static_doubleBlackFixup(this->l, this->l->p);
        }
    } else if (this->r) {
        *parentsChild = this->r;
        this->r->p = this->p;
        if (this->c == red) {
            // nothing, all good
        } else if (this->r->c == red) {
            this->r->c = black;
        } else if (this->c == black) {
            rbnode_static_doubleBlackFixup(this->r, this->r->p);
        }
    }
    // if it is a leaf, just fixup if need be
    else {
        *parentsChild = NULL;
        if (this->c == black) {
            rbnode_static_doubleBlackFixup(NULL, this->p);
        }
    }
    // clear the pointers so it dosen't think it has children when it destructs
    this->l = NULL;
    this->r = NULL;
    // destruct
    destroy_rbnode(this);
    free(this);
}

struct rbnode *rbnode_getUncle(struct rbnode *this) {
    if (!this->p->p)
        return NULL;
    int leftBranch = this->p == this->p->p->l;
    return leftBranch ? (this->p->p->r ? this->p->p->r : NULL) : (this->p->p->l ? this->p->p->l : NULL);
}

struct insertFixupInfo rbnode_getInsertCase(struct rbnode *this) {
    struct insertFixupInfo info = {0, 0};
    if (!this->p) {
        info.fc = root;
    } else if (!this->p->p || this->p->c == black) {
        info.fc = good;
    } else if (this->vtable->getUncle(this) && this->vtable->getUncle(this)->c == red) {
        info.fc = uncleRed;
        info.fo = (this->p == this->p->p->l) ? case_right : case_left;
    } else if (this == this->p->l && this->p == this->p->p->l) {
        info.fc = line;
        info.fo = case_left;
    } else if (this == this->p->r && this->p == this->p->p->l) {
        info.fc = triangle;
        info.fo = case_left;
    } else if (this == this->p->r && this->p == this->p->p->r) {
        info.fc = line;
        info.fo = case_right;
    } else if (this == this->p->l && this->p == this->p->p->r) {
        info.fc = triangle;
        info.fo = case_right;
    }
    return info;
};

void rbnode_insertFixup(struct rbnode *this) {
    int left = false;
    switch (this->vtable->getInsertCase(this).fc) {
    case root:
        this->c = black;
        break;
    case uncleRed:
        this->p->c = 1 - this->p->c;
        this->p->p->c = 1 - this->p->p->c;
        this->vtable->getUncle(this)->c = 1 - this->vtable->getUncle(this)->c;
        this->p->p->vtable->insertFixup(this->p->p);
        break;
    case line:
        left = (this->vtable->getInsertCase(this).fo == case_left);
        if (!left) {
            this->p->p->vtable->rotateLeft(this->p->p);
            this->p->c = black;
            this->p->l->c = red;
        } else {
            this->p->p->vtable->rotateRight(this->p->p);
            this->p->c = black;
            this->p->r->c = red;
        }
        break;
    case triangle:
        left = (this->vtable->getInsertCase(this).fo == case_left);
        if (!left) {
            this->p->vtable->rotateRight(this->p);
            this->r->vtable->insertFixup(this->r);
        } else {
            this->p->vtable->rotateLeft(this->p);
            this->l->vtable->insertFixup(this->l);
        }
        break;
    default:
        break;
    }
};

void rbnode_rotateLeft(struct rbnode *this) {
    struct rbnode *x = this;
    struct rbnode *y = this->r;
    this->r = y->l;
    if (y->l) {
        y->l->p = x;
    }
    y->p = x->p;
    if (!x->p) {
        this->vtable->setRoot(this, y);
    } else if (x == x->p->l) {
        x->p->l = y;
    } else {
        x->p->r = y;
    }
    y->l = x;
    x->p = y;
}

void rbnode_rotateRight(struct rbnode *this) {
    struct rbnode *x = this;
    struct rbnode *y = this->l;
    this->l = y->r;
    if (y->r) {
        y->r->p = x;
    }
    y->p = x->p;
    if (!x->p) {
        this->vtable->setRoot(this, y);
    } else if (x == x->p->r) {
        x->p->r = y;
    } else {
        x->p->l = y;
    }
    y->r = x;
    x->p = y;
}

int rbnode_static_isBlack(struct rbnode *n) {
    return !n || n->c == black;
}

int rbnode_static_isRed(struct rbnode *n) {
    return n && n->c == red;
}

struct rbnode *rbnode_static_getSibling(struct rbnode *dbnode, struct rbnode *parent) {
    if (!parent)
        return NULL;
    int left = dbnode == parent->l;
    return left ? (parent->r ? parent->r : NULL) : (parent->l ? parent->l : NULL);
}

struct rbnode *rbnode_static_getSiblingFarChild(struct rbnode *dbnode, struct rbnode *parent) {
    if (!parent)
        return NULL;
    int left = dbnode == parent->l;
    return left ? (parent->r ? parent->r->r : NULL) : (parent->l ? parent->l->l : NULL);
}

struct rbnode *rbnode_static_getSiblingNearChild(struct rbnode *dbnode, struct rbnode *parent) {
    if (!parent)
        return NULL;
    int left = dbnode == parent->l;
    return left ? (parent->r ? parent->r->l : NULL) : (parent->l ? parent->l->r : NULL);
}

void rbnode_static_doubleBlackFixup(struct rbnode *dbnode, struct rbnode *parent) {
    if (!dbnode && !parent) {
        return;
    } else if (dbnode && !dbnode->p) {
        return;
    }
    struct rbnode *sib = rbnode_static_getSibling(dbnode, parent);
    if (!sib) {
        printf("⚠️⚠️⚠️ shouldn't be here ⚠️⚠️⚠️\n");
        assert(false);
        return;
    } else if (rbnode_static_isBlack(sib) && rbnode_static_isBlack(sib->l) && rbnode_static_isBlack(sib->r)) {
        sib->c = red;
        if (rbnode_static_isBlack(parent)) {
            rbnode_static_doubleBlackFixup(parent, parent->p);
        } else {
            parent->c = black;
        }
    } else if (rbnode_static_isRed(sib)) {
        enum color t = sib->c;
        sib->c = parent->c;
        parent->c = t;
        (parent->r == dbnode ? parent->vtable->rotateRight(parent) : parent->vtable->rotateLeft(parent));
        rbnode_static_doubleBlackFixup(dbnode, parent);
    } else if (rbnode_static_isBlack(sib) && rbnode_static_isBlack(rbnode_static_getSiblingFarChild(dbnode, parent)) && rbnode_static_isRed(rbnode_static_getSiblingNearChild(dbnode, parent))) {
        rbnode_static_getSiblingNearChild(dbnode, parent)->c = black;
        sib->c = red;
        (sib == parent->r) ? sib->vtable->rotateRight(sib) : sib->vtable->rotateLeft(sib);
        rbnode_static_doubleBlackFixup(dbnode, parent);
    } else if (rbnode_static_isBlack(sib) && rbnode_static_isRed(rbnode_static_getSiblingFarChild(dbnode, parent))) {
        enum color t = sib->c;
        sib->c = parent->c;
        parent->c = t;
        rbnode_static_getSiblingFarChild(dbnode, parent)->c = black;
        (parent->r == dbnode ? parent->vtable->rotateRight(parent) : parent->vtable->rotateLeft(parent));
    }
}

void rbnode_setRoot(struct rbnode *this, struct rbnode *root) {
    *(this->root) = root;
}

void rbnode_printTree(struct rbnode *this, int depth, FILE * out, void (*printKey)(void * thing,FILE * out),void (*printVal)(void * thing,FILE * out)) {
    fprintf(out, " ");
    for (int i = 0; i < depth; i++) {
        fprintf(out, "│ ");
    }
    char *color;
    switch (this->c) {
    case red:
        color = "🟥";
        break;
    case black:
        color = "⬛";
        break;
    }
    fprintf(out, "%s", color);
    printKey(this->key, out);
    fprintf(out, " -> ");
    printVal(this->value, out);
    fprintf(out, "\n");

    if (this->r) {
        this->r->vtable->printTree(this->r, depth + 1, out, printKey, printVal);
    }
    if (this->l)
        this->l->vtable->printTree(this->l, depth + 1, out, printKey, printVal);
}

int rbnode_getKeys(struct rbnode *this, void * **keys) {

    int numRKeys = 0;
    int numLKeys = 0;
    void * *rkeys = NULL;
    void * *lkeys = NULL;

    if (this->r) {
        numRKeys = this->r->vtable->getKeys(this->r, &rkeys);
    }
    if (this->l) {
        numLKeys = this->l->vtable->getKeys(this->l, &lkeys);
    }

    *keys = calloc((numRKeys + numLKeys + 1), sizeof(void *));
    if (!(*keys)) {
        assert(false);
    }
    **keys = this->key;
    if (rkeys) {
        memcpy((*keys) + 1, rkeys, sizeof *rkeys * numRKeys);
        free(rkeys);
    }
    if (lkeys) {
        memcpy((*keys) + 1 + numRKeys, lkeys, sizeof *lkeys * numLKeys);
        free(lkeys);
    }

    return numRKeys + numLKeys + 1;
}

struct rbnode *rbnode_copy(struct rbnode *this, struct rbnode *parent) {
    // make a copy of this node
    struct rbnode *newme = construct_rbnode(this->key, this->value, parent, this->root);
    assert(newme);
    newme->c = this->c;
    // copy its subtrees, if they exist
    struct rbnode *newl = this->l ? this->l->vtable->copy(this->l, newme) : NULL;
    struct rbnode *newr = this->r ? this->r->vtable->copy(this->r, newme) : NULL;
    // update the pointers
    newme->l = newl;
    newme->r = newr;
    return newme;
}

int rbnode_contains(struct rbnode *this, void * key) {
    // if this is the node we are looking for, it is
    if (this->key == key) {
        return true;
    }
    // see if we need to look in the left or right subtree
    struct rbnode **recurseOn = (this->vtable->compare(this, key) >= 0) ? &(this->l) : &(this->r);
    // if that child dosen't exist, the node isn't in the tree,
    // otherwise, see if the subtree contains it
    if (!*recurseOn) {
        return false;
    }
    return (*recurseOn)->vtable->contains((*recurseOn), key);
}

void* rbnode_getValue(struct rbnode *this, void * key) {
    // if this is the node we are looking for, get it
    if (this->key == key) {
        return this->value;
    }
    // see if we need to look in the left or right subtree
    struct rbnode **recurseOn = (this->vtable->compare(this, key) >= 0) ? &this->l : &this->r;
    // if that child dosen't exist, the node isn't in the tree,
    // otherwise, get the value of the node in the subtree
    if (!*recurseOn) {
        return 0;
    }
    return (*recurseOn)->vtable->getValue((*recurseOn), key);
}

int rbnode_tryGetValue(struct rbnode *this, void * key, void * *out) {
    // if this is the node we are looking for, get it
    if (this->key == key) {
        *out = this->value;
        return true;
    }
    // see if we need to look in the left or right subtree
    struct rbnode **recurseOn = (this->vtable->compare(this, key) >= 0) ? &this->l : &this->r;
    // if that child dosen't exist, the node isn't in the tree,
    // otherwise, get the value of the node in the subtree
    if (!*recurseOn) {
        return false;
    }
    return (*recurseOn)->vtable->tryGetValue((*recurseOn), key, out);
}

int rbnode_setValue(struct rbnode *this, void * key, void * value) {
    // if this is the node we are looking for, update it
    if (this->key == key) {
        this->value = value;
        return true;
    }
    // see if we need to look in the left or right subtree
    struct rbnode **recurseOn = (this->vtable->compare(this, key) >= 0) ? &this->l : &this->r;
    // if that child dosen't exist, the node isn't in the tree,
    // otherwise, set the value of the node in the subtree
    if (!*recurseOn) {
        return false;
    }
    return (*recurseOn)->vtable->setValue((*recurseOn), key, value);
}

int rbnode_insert(struct rbnode *this, void * key, void * value) {
    // if we are they key, it is already in the tree
    if (this->vtable->compare(this, key) == 0)
        return false;
    // see if we need to look in the left or right subtree
    struct rbnode **recurseOn = (this->vtable->compare(this, key) >= 0) ? &(this->l) : &(this->r);
    // if that child dosen't exist, the node isn't in the tree, so insert it
    // otherwise, insert the node in the subtree
    if (!*recurseOn) {
        *recurseOn = construct_rbnode(key, value, this, this->root);
        assert(*recurseOn);
        (*recurseOn)->vtable->insertFixup((*recurseOn));
        return true;
    } else {
        return (*recurseOn)->vtable->insert((*recurseOn), key, value);
    }
}

int rbnode_remove(struct rbnode *this, void * key, struct rbnode **parentsChild) {
    // if this is the node we are looking for, delete
    if (this->key == key) {
        this->vtable->commitApoptosis(this, parentsChild);
        return true;
    }
    // see if we need to look in the left or right subtree
    struct rbnode **recurseOn = (this->vtable->compare(this, key) >= 0) ? &this->l : &this->r;
    // if that child dosen't exist, the node isn't in the tree
    // otherwise, remove the node from the subtree
    if (!*recurseOn) {
        return false;
    } else {
        return (*recurseOn)->vtable->remove((*recurseOn), key, recurseOn);
    }
}

int rbnode_bh(struct rbnode *this) {
    // if this is a leaf its bh is its color
    if (!this->l && !this->r) {
        return (this->c == red) ? 0 : 1;
    }
    // check left subtree, if it exists
    int lbh = this->l ? this->l->vtable->bh(this->l) : 0;
    // check right subtree, if it exists
    int rbh = this->r ? this->r->vtable->bh(this->r) : 0;
    // if one or both is broken, this is broken by default
    if (lbh == -1 || rbh == -1) {
        return -1;
    }
    // if they are equal, return their bh plus this node's color
    // else invalid
    return (lbh == rbh ? rbh + (this->c == red ? 0 : 1) : -1);
}

int rbnode_rp(struct rbnode *this) {
    // check left subtree, if it exists
    int lrp = this->l ? this->l->vtable->rp(this->l) : true;
    // check right subtree, if it exists
    int rrp = this->r ? this->r->vtable->rp(this->r) : true;
    // make sure this and its children don't make a new violation
    // then and everything together
    return (lrp && rrp && (this->c == red ? ((this->l ? this->l->c != red : true) && (this->r ? this->r->c != red : true)) : true));
}

void destroy_rbnode(struct rbnode *this) {
    if (this->r) {
        destroy_rbnode(this->r);
        free(this->r);
    }
    if (this->l) {
        destroy_rbnode(this->l);
        free(this->l);
    }
}

struct rbnode_VTABLE_s rbnode_VTABLE = {
    .compare = rbnode_compare,
    .commitApoptosis = rbnode_commitApoptosis,
    .getUncle = rbnode_getUncle,
    .getInsertCase = rbnode_getInsertCase,
    .insertFixup = rbnode_insertFixup,
    .rotateLeft = rbnode_rotateLeft,
    .rotateRight = rbnode_rotateRight,
    .static_isBlack = rbnode_static_isBlack,
    .static_isRed = rbnode_static_isRed,
    .static_getSibling = rbnode_static_getSibling,
    .static_getSiblingFarChild = rbnode_static_getSiblingFarChild,
    .static_getSiblingNearChild = rbnode_static_getSiblingNearChild,
    .static_doubleBlackFixup = rbnode_static_doubleBlackFixup,
    .setRoot = rbnode_setRoot,
    .printTree = rbnode_printTree,
    .getKeys = rbnode_getKeys,
    .copy = rbnode_copy,
    .contains = rbnode_contains,
    .getValue = rbnode_getValue,
    .tryGetValue = rbnode_tryGetValue,
    .setValue = rbnode_setValue,
    .insert = rbnode_insert,
    .remove = rbnode_remove,
    .bh = rbnode_bh,
    .rp = rbnode_rp,
};

struct rbnode *construct_rbnode(void * key, void * value, struct rbnode *parent, struct rbnode **r) {
    struct rbnode *node = calloc(1, sizeof *node);
    node->l = NULL;
    node->r = NULL;

    node->root = r;
    node->p = parent;
    node->key = key;
    node->value = value;
    if (!parent) {
        node->c = black;
    }

    node->vtable = &rbnode_VTABLE;
    // printf("%p\n", &rbnode_VTABLE);

    return node;
};
