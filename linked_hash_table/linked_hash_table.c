#include "./linked_hash_table.h"
#include "./chaining_list/chaining_list.h"
#include <stdint.h>

struct LinkedHashTable * construct_LinkedHashTable() {
    struct LinkedHashTable * new = malloc(sizeof *new);
    new->numSlots = 10;
    new->numKeys = 0;
    new->alphaMax = .75;
    new->underlyingArray = malloc(sizeof *new->underlyingArray * new->numSlots);
    for(int i = 0; i < new->numSlots; ++i) {
        new->underlyingArray[i] = construct_ChainingList();
    }
    return(struct LinkedHashTable *)new;
}

unsigned int LinkedHashTable_hash(struct LinkedHashTable * this, unsigned int x) {
    (void)this;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

float LinkedHashTable_loadFactor(struct LinkedHashTable * this) {
    return (float)this->numKeys / (float)this->numSlots;
}

void LinkedHashTable_resize(struct LinkedHashTable * this, float factor) {
    if (this->numSlots <= 10 && factor < 1) {
        return;
    }
    int prevSlots = this->numSlots;
    this->numSlots = this->numSlots * factor;
    this->numKeys = 0;
    struct ChainingList **oldUnderlyingArray = this->underlyingArray;
    this->underlyingArray = malloc(sizeof *this->underlyingArray * this->numSlots);
    for(int i = 0; i < this->numSlots; ++i) {
        this->underlyingArray[i] = construct_ChainingList();
    }
    for (int i = 0; i < prevSlots; i++) {
        struct ChainingList **temp = oldUnderlyingArray + i;
        struct ChainingNode *curr = (*temp)->root;

        while (curr != NULL) {
            LinkedHashTable_insert((struct LinkedHashTable *)this, curr->key, curr->value);
            curr = curr->next;
        }
    }
    for(int i = 0; i < prevSlots; ++i) {
        destroy_ChainingList(oldUnderlyingArray[i]);
        free(oldUnderlyingArray[i]);
    }
    free(oldUnderlyingArray);
}

void LinkedHashTable_print(struct LinkedHashTable * mapThis, FILE * out,void (*printKey)(void * thing,FILE * out),void (*printVal)(void * thing,FILE * out)) {
    struct LinkedHashTable * this = (struct LinkedHashTable *)mapThis;
    fprintf(out, "{");
    for (int i = 0; i < this->numSlots; i++) {
        void **keys = malloc(sizeof(void *) * this->numKeys);
        int len = ChainingList_keys(this->underlyingArray[i], keys);
        for(int j = 0; j < len; ++j) {
            void * val;
            ChainingList_tryGet(this->underlyingArray[i], keys[j], &val);
            printKey(keys[j], out);
            fprintf(out, ": ");
            printVal(val, out);
            fprintf(out, ", ");
        }
    }
    fprintf(out, "}\n");
}

int LinkedHashTable_getKeys(struct LinkedHashTable * mapThis, void **keys) {
    struct LinkedHashTable * this = (struct LinkedHashTable *)mapThis;

    int inPtr = 0;
    for (int i = 0; i < this->numSlots; i++) {
        inPtr += ChainingList_keys(this->underlyingArray[i], keys + inPtr);
    }
    return this->numKeys;
};

void LinkedHashTable_clear(struct LinkedHashTable * mapThis) {
    struct LinkedHashTable * this = (struct LinkedHashTable *)mapThis;
    for (int i = 0; i < this->numSlots; i++) {
        destroy_ChainingList(this->underlyingArray[i]);
    }
    this->numKeys = 0;
}

struct LinkedHashTable * LinkedHashTable_copy(struct LinkedHashTable * mapThis) {

    struct LinkedHashTable * newTable = construct_LinkedHashTable();

    int numKeys = mapThis->numKeys;
    void **keys = malloc(numKeys * sizeof * keys);
    LinkedHashTable_getKeys(mapThis, keys);

    for (int i = 0; i < numKeys; ++i) {
        LinkedHashTable_insert(newTable, keys[i], LinkedHashTable_getValue(mapThis, keys[i]));
    }

    return newTable;
}

bool LinkedHashTable_contains(struct LinkedHashTable * mapThis, void * key) {
    struct LinkedHashTable * this = (struct LinkedHashTable *)mapThis;

    void * temp;
    int slot = LinkedHashTable_hash(this, (uintptr_t)key) % this->numSlots;
    return ChainingList_tryGet(this->underlyingArray[slot],key, &temp);
}

void * LinkedHashTable_getValue(struct LinkedHashTable * mapThis, void * key) {
    struct LinkedHashTable * this = (struct LinkedHashTable *)mapThis;

    void * out;
    int slot = LinkedHashTable_hash(this, (uintptr_t)key) % this->numSlots;
    ChainingList_tryGet(this->underlyingArray[slot], key, &out);
    return out;
}

bool LinkedHashTable_tryGetValue(struct LinkedHashTable * mapThis, void * key, void **out) {
    struct LinkedHashTable * this = (struct LinkedHashTable *)mapThis;

    int slot = LinkedHashTable_hash(this, (uintptr_t)key) % this->numSlots;
    return ChainingList_tryGet(this->underlyingArray[slot], key, out);
}

bool LinkedHashTable_setValue(struct LinkedHashTable * mapThis, void * key, void * value) {
    struct LinkedHashTable * this = (struct LinkedHashTable *)mapThis;

    int slot = LinkedHashTable_hash(this, (uintptr_t)key) % this->numSlots;
    return ChainingList_setValue(this->underlyingArray[slot], key, value);
}

bool LinkedHashTable_insert(struct LinkedHashTable * mapThis, void * key, void * value) {
    struct LinkedHashTable * this = (struct LinkedHashTable *)mapThis;
    int slot = LinkedHashTable_hash(this, (uintptr_t)key) % this->numSlots;
    bool ret = ChainingList_insert(this->underlyingArray[slot], key, value);
    if (ret) {
        this->numKeys++;
    }
    if (LinkedHashTable_loadFactor(this) - this->alphaMax > 0) {
        LinkedHashTable_resize(this, 2);
    }
    return ret;
}

bool LinkedHashTable_remove(struct LinkedHashTable * mapThis, void * key) {
    struct LinkedHashTable * this = (struct LinkedHashTable *)mapThis;

    int slot = LinkedHashTable_hash(this, (uintptr_t)key) % this->numSlots;
    bool ret = ChainingList_remove(this->underlyingArray[slot], key);
    if (ret) {
        this->numKeys--;
    }
    if (LinkedHashTable_loadFactor(this) - (this->alphaMax / 4) < 0) {
        LinkedHashTable_resize(this, .5);
    }
    return ret;
}

void destroy_LinkedHashTable(struct LinkedHashTable * this) {
    for (int i = 0; i < this->numSlots; i++) {
        destroy_ChainingList(this->underlyingArray[i]);
        free(this->underlyingArray[i]);
    }
    free(this->underlyingArray);
}
size_t LinkedHashTable_getSize(struct LinkedHashTable * this) {
    // mem footprint
    int mem = sizeof *this + sizeof *this->underlyingArray * this->numSlots;
    for(int i = 0; i < this->numSlots; ++i) {
        mem += ChainingList_size(this->underlyingArray[i]);
    }
    return mem;
}


// static struct LinkedHashTable vTable = {
//     .print = LinkedHashTable_print,
//     .clear = LinkedHashTable_clear,
//     .contains = LinkedHashTable_contains,
//     .destroy = destroy_LinkedHashTable,
//     .getKeys = LinkedHashTable_getKeys,
//     .getValue = LinkedHashTable_getValue,
//     .tryGetValue = LinkedHashTable_tryGetValue,
//     .insert = LinkedHashTable_insert,
//     .remove = LinkedHashTable_remove,
//     .setValue = LinkedHashTable_setValue,
//     .copy = LinkedHashTable_copy
// };
