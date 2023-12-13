#pragma once
#include <string.h>
#include <stdio.h>
#include <stdbool.h>

struct LinkedHashTable {
    int numSlots;
    int numKeys;
    float alphaMax;
    struct ChainingList ** underlyingArray;
};

unsigned int LinkedHashTable_hash(struct LinkedHashTable * this, unsigned int x);
float LinkedHashTable_loadFactor(struct LinkedHashTable * this);
void LinkedHashTable_resize(struct LinkedHashTable * this, float factor);
void LinkedHashTable_print(struct LinkedHashTable * mapThis, FILE * out,void (*printKey)(void * thing,FILE * out),void (*printVal)(void * thing,FILE * out));
int LinkedHashTable_getKeys(struct LinkedHashTable * mapThis, void **keys);
void LinkedHashTable_clear(struct LinkedHashTable * mapThis);
struct LinkedHashTable * LinkedHashTable_copy(struct LinkedHashTable * mapThis);
bool LinkedHashTable_contains(struct LinkedHashTable * mapThis, void * key);
void * LinkedHashTable_getValue(struct LinkedHashTable * mapThis, void * key);
bool LinkedHashTable_tryGetValue(struct LinkedHashTable * mapThis, void * key, void **out);
bool LinkedHashTable_setValue(struct LinkedHashTable * mapThis, void * key, void * value);
bool LinkedHashTable_insert(struct LinkedHashTable * mapThis, void * key, void * value);
bool LinkedHashTable_remove(struct LinkedHashTable * mapThis, void * key);
void destroy_LinkedHashTable(struct LinkedHashTable * this);
struct LinkedHashTable * construct_LinkedHashTable();
size_t LinkedHashTable_getSize(struct LinkedHashTable * this);