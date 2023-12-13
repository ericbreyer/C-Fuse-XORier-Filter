#include "chaining_list.h"
#include <stdio.h>
struct ChainingNode * construct_ChainingNode(void * key, void * value) {
    struct ChainingNode * new = malloc(sizeof *new);
    new->key = key;
    new->value = value;
    new->next = NULL;
    return new;
}

bool ChainingList_insert(struct ChainingList * this, void * key, void * value) {
    struct ChainingNode *prev = NULL;
    struct ChainingNode *curr = this->root;
    while (curr != NULL) {
        if (key == curr->key) {
            return false;
        }
        prev = curr;
        curr = curr->next;
    }
    struct ChainingNode *newChainingNode = construct_ChainingNode(key, value);
    if (!this->root) { // root is null
        this->root = newChainingNode;
    } else {
        prev->next = newChainingNode;
    }
    return true;
}

bool ChainingList_setValue(struct ChainingList * this, void * key, void * value) {
    // struct ChainingNode *prev = NULL;
    struct ChainingNode *curr = this->root;
    while (curr != NULL) {
        if (key == curr->key) {
            curr->value = value;
            return true;
        }
        // prev = curr;
        curr = curr->next;
    }
    return false;
}

bool ChainingList_remove(struct ChainingList * this, void * key) {
    struct ChainingNode *prev = NULL;
    struct ChainingNode *curr = this->root;
    while (curr != NULL) {
        if (key == curr->key) {
            break;
        }
        prev = curr;
        curr = curr->next;
    }
    if (prev != NULL && curr != NULL) { // key was found and is not root
        prev->next = curr->next;
        free(curr);
        curr = NULL;
    } else if (prev == NULL && curr != NULL) { // key was root
        this->root = curr->next;
        free(curr);
    } else { // key was not found
        return false;
    }
    return true;
}

inline bool ChainingList_tryGet(struct ChainingList * this, void * key, void * *out) {
    struct ChainingNode *curr = this->root;
    while (curr != NULL) {
        if (key == curr->key) {
            *out = curr->value;
            return true;
        }
        curr = curr->next;
    }
    return false;
}

int ChainingList_keys(struct ChainingList * this, void * *keys) {
    struct ChainingNode *curr = this->root;
    int numKeys = 0;
    while (curr != NULL) {
        *(keys + numKeys) = curr->key;
        ++numKeys;
        curr = curr->next;
    }
    return numKeys;
}

struct ChainingList * construct_ChainingList() {
    struct ChainingList * new = malloc(sizeof *new);
    new->root = NULL;
    return new;
    
}

void destroy_ChainingList(struct ChainingList * this) {
    struct ChainingNode *prev = NULL;
    struct ChainingNode *curr = this->root;
    while (curr != NULL) {
        prev = curr;
        curr = curr->next;
        free(prev);
    }
    this->root = NULL;
}

size_t ChainingList_size(struct ChainingList * this) {
    // memory footprint
    int numKeys = 0;
    struct ChainingNode *curr = this->root;
    while (curr != NULL) {
        ++numKeys;
        curr = curr->next;
    }
    return sizeof *this + numKeys * sizeof * curr;

    
}