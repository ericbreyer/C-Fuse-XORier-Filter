#pragma once
#include <string.h>
#include <stdio.h>
#include <stdbool.h>
#include <sized_pointer.h>

enum SlotState {
    kEmpty,
    kDeleted,
    kOccupied
};

struct Node {
    enum SlotState state;
    SizedPointer key;
    void * value;
};

struct HashTable {
    size_t numSlots;
    size_t numKeys;
    float alphaMax;
    struct Node * underlyingArray;
};

// unsigned int HashTable_hash(struct HashTable * this, unsigned int x);
float HashTable_loadFactor(struct HashTable * this);
void HashTable_resize(struct HashTable * this, float factor);
void HashTable_print(struct HashTable * mapThis, FILE * out,void (*printKey)(SizedPointer thing,FILE * out),void (*printVal)(void * thing,FILE * out));
size_t HashTable_getKeys(struct HashTable * mapThis, SizedPointer *keys);
void HashTable_clear(struct HashTable * mapThis);
struct HashTable * HashTable_copy(struct HashTable * mapThis);
bool HashTable_contains(struct HashTable * mapThis, SizedPointer key);
void * HashTable_getValue(struct HashTable * mapThis, SizedPointer key);
bool HashTable_tryGetValue(struct HashTable * mapThis, SizedPointer key, void **out);
bool HashTable_setValue(struct HashTable * mapThis, SizedPointer key, void * value);
bool HashTable_insert(struct HashTable * mapThis, SizedPointer key, void * value);
bool HashTable_remove(struct HashTable * mapThis, SizedPointer key);
void destroy_HashTable(struct HashTable * this);
struct HashTable * construct_HashTable();
struct HashTable *construct_HashTableSize(size_t size);
size_t HashTable_getSize(struct HashTable * this);
