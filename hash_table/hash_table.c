#include "./hash_table.h"
#include <stdint.h>
#include <stdlib.h>

// hash table using open addressing


static inline uint32_t murmur_32_scramble(uint32_t k) {
    k *= 0xcc9e2d51;
    k = (k << 15) | (k >> 17);
    k *= 0x1b873593;
    return k;
}
static uint32_t murmur3_32(const uint8_t *key, size_t len, uint32_t seed) {
    uint32_t h = seed;
    uint32_t k;
    /* Read in groups of 4. */
    for (size_t i = len >> 2; i; i--) {
        // Here is a source of differing results across endiannesses.
        // A swap here has no effects on hash properties though.
        memcpy(&k, key, sizeof(uint32_t));
        key += sizeof(uint32_t);
        h ^= murmur_32_scramble(k);
        h = (h << 13) | (h >> 19);
        h = h * 5 + 0xe6546b64;
    }
    /* Read the rest. */
    k = 0;
    for (size_t i = len & 3; i; i--) {
        k <<= 8;
        k |= key[i - 1];
    }
    // A swap is *not* necessary here because the preceding loop already
    // places the low bytes in the low places according to whatever endianness
    // we use. Swaps only apply when the memory is copied in a chunk.
    h ^= murmur_32_scramble(k);
    /* Finalize. */
    h ^= len & ((uint32_t)-1);
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
}

static uint64_t roundUp2(uint64_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v |= v >> 32;
    v++;
    return v;
}

size_t getSlot(size_t hash, size_t i, size_t mod) {
    // return (int)(hash + .5 * (float)i + .5 * (float)(i * i)) % mod;
    return (size_t)(hash + i) % mod;
}

float HashTable_loadFactor(struct HashTable *this) {
    return (float)this->numKeys / (float)this->numSlots;
}

void HashTable_resize(struct HashTable *this, float factor) {
    if (this->numSlots <= 8 && factor < 1) {
        return;
    }
    size_t prevSlots = this->numSlots;
    this->numSlots = (size_t)((float)this->numSlots * factor);
    this->numKeys = 0;
    struct Node *oldUnderlyingArray = this->underlyingArray;
    this->underlyingArray =
        malloc(sizeof *this->underlyingArray * this->numSlots);
    if (this->underlyingArray == NULL) {
        perror("malloc");
        exit(1);
    }
    for (size_t i = 0; i < this->numSlots; i++) {
        this->underlyingArray[i].state = kEmpty;
    }
    for (size_t i = 0; i < prevSlots; i++) {
        struct Node curr = oldUnderlyingArray[i];
        if (curr.state == kOccupied) {
            HashTable_insert((struct HashTable *)this, curr.key, curr.value);
        }
    }
    free(oldUnderlyingArray);
}

void HashTable_print(struct HashTable *mapThis, FILE *out,
                     void (*printKey)(SizedPointer thing, FILE *out),
                     void (*printVal)(void *thing, FILE *out)) {
    struct HashTable *this = (struct HashTable *)mapThis;
    fprintf(out, "{");
    for (size_t i = 0; i < this->numSlots; i++) {
        if (this->underlyingArray[i].state == kOccupied) {
            printKey(this->underlyingArray[i].key, out);
            fprintf(out, ": ");
            printVal(this->underlyingArray[i].value, out);
            fprintf(out, ", ");
        }
    }
    fprintf(out, "}\n");
}

size_t HashTable_getKeys(struct HashTable *mapThis, SizedPointer *keys) {
    struct HashTable *this = (struct HashTable *)mapThis;

    int inPtr = 0;
    for (size_t i = 0; i < this->numSlots; i++) {
        if (this->underlyingArray[i].state == kOccupied) {
            keys[inPtr] = this->underlyingArray[i].key;
            inPtr++;
        }
    }
    return this->numKeys;
};

void HashTable_clear(struct HashTable *mapThis) {
    struct HashTable *this = (struct HashTable *)mapThis;
    for (size_t i = 0; i < this->numSlots; i++) {
        this->underlyingArray[i].state = kEmpty;
    }
    this->numKeys = 0;
}

struct HashTable *HashTable_copy(struct HashTable *mapThis) {
    struct HashTable *this = (struct HashTable *)mapThis;

    struct HashTable *newTable = construct_HashTable(this->numSlots);

    for (size_t i = 0; i < this->numSlots; ++i) {
        if (this->underlyingArray[i].state == kOccupied) {
            HashTable_insert(newTable, this->underlyingArray[i].key,
                             this->underlyingArray[i].value);
        }
    }
    return newTable;
}

bool HashTable_contains(struct HashTable *mapThis, SizedPointer key) {
    struct HashTable *this = (struct HashTable *)mapThis;

    // size_t slot = HashTable_hash(this, (uintptr_t)key) % this->numSlots;
    size_t slot = murmur3_32(key.ptr, key.size, 0) % this->numSlots;

    size_t firstDeletedSlot = (size_t)-1;

    for (size_t i = 0; i < this->numSlots; ++i) {
        size_t tSlot = getSlot(slot, i, this->numSlots);

        if (this->underlyingArray[tSlot].state == kEmpty) {
            return false;
        }

        if (this->underlyingArray[tSlot].state == kDeleted) {
            if (firstDeletedSlot == (size_t)-1) {
                firstDeletedSlot = tSlot;
            }
            continue;
        }

        if (this->underlyingArray[tSlot].key.ptr == key.ptr &&
            this->underlyingArray[tSlot].state == kOccupied) {
            if (firstDeletedSlot != (size_t)-1) {
                this->underlyingArray[firstDeletedSlot].key =
                    this->underlyingArray[tSlot].key;
                this->underlyingArray[firstDeletedSlot].value =
                    this->underlyingArray[tSlot].value;
                this->underlyingArray[firstDeletedSlot].state = kOccupied;
                this->underlyingArray[tSlot].state = kDeleted;
            }
            return true;
        }
    }

    return false;
}

void *HashTable_getValue(struct HashTable *mapThis, SizedPointer key) {
    struct HashTable *this = (struct HashTable *)mapThis;

    // size_t slot = HashTable_hash(this, (uintptr_t)key) % this->numSlots;
    size_t slot = murmur3_32(key.ptr, key.size, 0) % this->numSlots;

    size_t firstDeletedSlot = (size_t)-1;

    for (size_t i = 0; i < this->numSlots; ++i) {
        size_t tSlot = getSlot(slot, i, this->numSlots);
        if (this->underlyingArray[tSlot].state == kEmpty) {
            return NULL;
        }
        if (this->underlyingArray[tSlot].state == kDeleted &&
            firstDeletedSlot == (size_t)-1) {
            firstDeletedSlot = tSlot;
            continue;
        }
        if (this->underlyingArray[tSlot].key.ptr == key.ptr &&
            this->underlyingArray[tSlot].state == kOccupied) {
            if (firstDeletedSlot != (size_t)-1) {
                this->underlyingArray[firstDeletedSlot].key =
                    this->underlyingArray[tSlot].key;
                this->underlyingArray[firstDeletedSlot].value =
                    this->underlyingArray[tSlot].value;
                this->underlyingArray[firstDeletedSlot].state = kOccupied;
                this->underlyingArray[tSlot].state = kDeleted;
            }
            return this->underlyingArray[tSlot].value;
        }
    }

    return NULL;
}

bool HashTable_tryGetValue(struct HashTable *mapThis, SizedPointer key, void **out) {
    struct HashTable *this = (struct HashTable *)mapThis;

    // size_t slot = HashTable_hash(this, (uintptr_t)key) % this->numSlots;
    size_t slot = murmur3_32(key.ptr, key.size, 0) % this->numSlots;

    size_t firstDeletedSlot = (size_t)-1;

    for (size_t i = 0; i < this->numSlots; ++i) {
        size_t tSlot = getSlot(slot, i, this->numSlots);

        if (this->underlyingArray[tSlot].state == kEmpty) {
            return false;
        }

        if (this->underlyingArray[tSlot].state == kDeleted) {
            if (firstDeletedSlot == (size_t)-1) {
                firstDeletedSlot = tSlot;
            }
            continue;
        }

        if (this->underlyingArray[tSlot].key.ptr == key.ptr &&
            this->underlyingArray[tSlot].state == kOccupied) {
            if (firstDeletedSlot != (size_t)-1) {
                this->underlyingArray[firstDeletedSlot].key =
                    this->underlyingArray[tSlot].key;
                this->underlyingArray[firstDeletedSlot].value =
                    this->underlyingArray[tSlot].value;
                this->underlyingArray[firstDeletedSlot].state = kOccupied;
                this->underlyingArray[tSlot].state = kDeleted;
            }
            *out = this->underlyingArray[tSlot].value;
            return true;
        }
    }

    return false;
}
bool HashTable_setValue(struct HashTable *mapThis, SizedPointer key, void *value) {
    struct HashTable *this = (struct HashTable *)mapThis;

    // size_t slot = HashTable_hash(this, (uintptr_t)key) % this->numSlots;
    size_t slot = murmur3_32(key.ptr, key.size, 0) % this->numSlots;

    size_t firstDeletedSlot = (size_t)-1;

    for (size_t i = 0; i < this->numSlots; ++i) {
        size_t tSlot = getSlot(slot, i, this->numSlots);
        if (this->underlyingArray[tSlot].state == kEmpty) {
            return false;
        }
        if (this->underlyingArray[tSlot].state == kDeleted &&
            firstDeletedSlot == (size_t)-1) {
            firstDeletedSlot = tSlot;
        }
        if (this->underlyingArray[tSlot].key.ptr == key.ptr &&
            this->underlyingArray[tSlot].state == kOccupied) {
            if (firstDeletedSlot == (size_t)-1) {
                this->underlyingArray[tSlot].value = value;
                return true;
            }
            this->underlyingArray[firstDeletedSlot].key = key;
            this->underlyingArray[firstDeletedSlot].value = value;
            this->underlyingArray[firstDeletedSlot].state = kOccupied;
            this->underlyingArray[tSlot].state = kDeleted;
            return true;
        }
    }
    return false;
}

bool HashTable_insert(struct HashTable *mapThis, SizedPointer key, void *value) {
    struct HashTable *this = (struct HashTable *)mapThis;
    // size_t slot = HashTable_hash(this, (uintptr_t)key) % this->numSlots;
    size_t slot = murmur3_32(key.ptr, key.size, 0) % this->numSlots;

    for (size_t i = 0; i < this->numSlots; ++i) {
        size_t tSlot = getSlot(slot, i, this->numSlots);
        if (this->underlyingArray[tSlot].state == kEmpty ||
            this->underlyingArray[tSlot].state == kDeleted) {
            this->underlyingArray[tSlot].key = key;
            this->underlyingArray[tSlot].value = value;
            this->underlyingArray[tSlot].state = kOccupied;
            this->numKeys++;
            if (HashTable_loadFactor(this) - this->alphaMax > 0) {
                HashTable_resize(this, 2);
            }
            return true;
        }
        if (this->underlyingArray[tSlot].key.ptr == key.ptr &&
            this->underlyingArray[tSlot].state == kOccupied) {
            return false;
        }
    }

    return false;
}

bool HashTable_remove(struct HashTable *mapThis, SizedPointer key) {
    struct HashTable *this = (struct HashTable *)mapThis;

    // size_t slot = HashTable_hash(this, (uintptr_t)key) % this->numSlots;
    size_t slot = murmur3_32(key.ptr, key.size, 0) % this->numSlots;
    for (size_t i = 0; i < this->numSlots; ++i) {
        size_t tSlot = getSlot(slot, i, this->numSlots);
        if (this->underlyingArray[tSlot].state == kEmpty) {
            return false;
        }
        if (this->underlyingArray[tSlot].key.ptr == key.ptr &&
            this->underlyingArray[tSlot].state == kOccupied) {
            this->underlyingArray[tSlot].state = kDeleted;
            this->numKeys--;
            if (HashTable_loadFactor(this) - (this->alphaMax / 4) < 0) {
                HashTable_resize(this, .5);
            }
            return true;
        }
    }
    return false;
}

void destroy_HashTable(struct HashTable *this) { free(this->underlyingArray); }

struct HashTable *construct_HashTable() {
    struct HashTable *new = malloc(sizeof *new);
    new->numSlots = 8;
    new->numKeys = 0;
    new->alphaMax = .75;
    new->underlyingArray =
        calloc(1, sizeof *new->underlyingArray *new->numSlots);
    return (struct HashTable *)new;
}

struct HashTable *construct_HashTableSize(size_t size) {
    struct HashTable *new = malloc(sizeof *new);
    new->numSlots = roundUp2(size);
    // fprintf(stderr,"numSlots: %d\n", new->numSlots);
    new->numKeys = 0;
    new->alphaMax = .75;
    new->underlyingArray =
        calloc(1, sizeof *new->underlyingArray *new->numSlots);
    return (struct HashTable *)new;
}

size_t HashTable_getSize(struct HashTable *this) { 
    // memory footprint
    return sizeof *this + sizeof *this->underlyingArray * this->numSlots;
 }