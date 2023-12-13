//------------------------------------------------------------------------------
// Information
//------------------------------------------------------------------------------

/**
 * @file fuse_xorier_lookup_table.h
 * @author Eric Breyer (eric.breyer@gmail.com) [https://github.com/ericbreyer]
 * @brief Declaration of a fuse XORier lookup table
 * @version 1.0
 * @date 2023-12-13
 *
 * @copyright Copyright (c) 2023, released under GNU Public Licence 3.0 or later
 *
 */
// SPDX-License-Identifier: GNU-3.0-or-later

//------------------------------------------------------------------------------
// Includes
//------------------------------------------------------------------------------

#include <math.h>   // pow
#include <stdio.h>  // logging
#include <stdlib.h> // free and malloc
#include <string.h> // memcpy

#include "config.h"                    // configuration options
#include "fused_XORier_lookup_table.h" // interface

//------------------------------------------------------------------------------
// Types
//------------------------------------------------------------------------------

struct fuseXORierLookupTable {
    const size_t k;         // number of hash functions
    uint32_t *hash_seeds;   // the k + 2 hash seeds the filter needs
    const size_t m;         // number of slots in the filter
    const int q;            // number of bits in a slot
    slot_t *const table1;   // the underlying array structure of the filter used
                            // to test membership
    void **const table2;    // the stored values
    const bool verbose;     // verbose printing (for debugging)
    const bool cacheHashes; // whether or not to cache hashes

    const size_t numSegments; // the number of segments, or 0 if this filter not
                              // fused (just a bloomier filter)
    const size_t segmentSize; // the segment size

    const size_t
        n; // number of elements inserted (for bookkeeping/debugging purposes)
};

/**
 * @brief A set of elements, represented as the xor of said elements and the
 * cardinality of the set. Useful for finding singletons.
 *
 */
struct xorSet {
    size_t elems; // the xor of the elements in the set
    size_t cardinality; // the number of elements in the set
};

/**
 * @brief A circular buffer of size_t's. Can grow dynamically if needed. Used as
 * a queue.
 *
 */
struct circularBuffer {
    size_t *elems; // the elements in the queue
    size_t size;  // the size of the queue
    size_t insert; // the index to insert the next element
    size_t numElems; // the number of elements in the queue
};

/**
 * @brief A pair struct to hold the hashes of an element
 *
 */
struct itemHash {
    bool hashesValid; // whether or not the hashes are valid or need to be (re)computed
    size_t *hashes; // the hashes of the element
};

/**
 * @brief A pair struct to hold the result of fuseXORierLT_findMatch. Namely the
 * index of the key and the index of the hash function that maps it to
 * its location.
 *
 */
struct match {
    size_t keyIdx; // the index of the key
    size_t l; // the index of the hash function that maps the key to its location
};

/**
 * @brief The result of fuseXORierLT_findMatch. Namely the ordering of the keys
 * and the mapping from keys to locations.
 *
 */
struct mappingResult {
    struct match *ordering; // the ordering of the keys (and their mappings)
    bool success; // whether or not the matching succeeded
};

//------------------------------------------------------------------------------
// Forward Declarations of Helper Functions
//------------------------------------------------------------------------------

//-------------
//---Hashing---
//-------------

/**
 * @brief implements the murmur3_32 hash algorithm
 * 
 * @param[in] key a pointer to the key
 * @param[in] len the length of the key
 * @param[in] seed the seed for the hash
 * @return uint32_t representing the hash
 */
static uint32_t murmur3_32(const uint8_t *key, size_t len, uint32_t seed);

//-------------
//---XOR Sets--
//-------------

/**
 * @brief add an element to a xorSet
 * 
 * @param[in] set said xorSet
 * @param[in] elem the element to add
 */
void addElem(struct xorSet *set, size_t elem);

/**
 * @brief remove an element from a xorSet
 * 
 * @param[in] set said xorSet
 * @param[in] elem the element to remove
 */
void removeElem(struct xorSet *set, size_t elem);

/**
 * @brief check if a xorSet is a singleton
 * 
 * @param[in] set said xorSet
 * @return true if the xorSet is a singleton
 * @return false otherwise
 */
bool singleton(struct xorSet *set);

/**
 * @brief get the singleton element from a xorSet
 * 
 * @param[in] set said xorSet
 * @return size_t the singleton element
 */
size_t getSingleton(struct xorSet *set);

//---------------------
//---Circular Buffer---
//---------------------

/**
 * @brief construct a circular buffer
 * 
 * @param[in] self the circular buffer
 * @param[in] size the initial size of the circular buffer
 */
void construct_circularBuffer(struct circularBuffer *self, size_t size);

/**
 * @brief destroy a circular buffer
 * 
 * @param[in] self the circular buffer
 */
void destroy_circularBuffer(struct circularBuffer *self);

/**
 * @brief check if a circular buffer is empty
 * 
 * @param[in] self the circular buffer
 * @return true if the circular buffer is empty
 * @return false otherwise
 */
bool circularBuffer_isEmpty(struct circularBuffer *self);

/**
 * @brief get the number of elements in a circular buffer
 * 
 * @param[in] self the circular buffer
 * @return size_t the number of elements in the circular buffer
 */
size_t circularBuffer_numElems(struct circularBuffer *self);

/**
 * @brief grow a circular buffer
 * 
 * @param[in] self the circular buffer
 */
void circularBuffer_grow(struct circularBuffer *self);

/**
 * @brief enqueue an element into a circular buffer
 * 
 * @param[in] self the circular buffer
 * @param[in] elem the element to enqueue
 */
void circularBuffer_enqueue(struct circularBuffer *self, size_t elem);

/**
 * @brief dequeue an element from a circular buffer
 * 
 * @param[in] self the circular buffer
 * @return size_t the dequeued element
 */
size_t circularBuffer_dequeue(struct circularBuffer *self);

//----------
//---Misc---
//----------

/**
 * @brief round a number up to the next power of 2
 * 
 * @param[in] v said number
 * @return uint64_t the next power of 2
 */
static uint64_t roundUp2(uint64_t v);

//------------------------------------------------------------------------------
// FXLT Helper Functions
//------------------------------------------------------------------------------

/**
 * @brief Hashes an element k + 1 times and stores the results in neighborhood
 * given the hash seeds in self
 *
 * @param[in] self the fxlt
 * @param[in] elem the element to hash
 * @param[out] neighborhood an array of size k + 1 to store the results
 */
void fuseXORierLT_hashAll(struct fuseXORierLookupTable *self, SizedPointer elem,
                          size_t *neighborhood) {
    if (self->numSegments == 0) {
        for (size_t i = 0; i < self->k; ++i) {
            neighborhood[i] =
                murmur3_32(elem.ptr, elem.size, self->hash_seeds[i]) % self->m;
        }
    } else {

        size_t firstSegmentStart =
            (murmur3_32(elem.ptr, elem.size, self->hash_seeds[self->k + 1]) %
             (self->numSegments - self->k)) *
            self->segmentSize;

        for (size_t i = 0; i < self->k; ++i) {
            size_t segmentStart = firstSegmentStart + i * self->segmentSize;
            size_t temp = murmur3_32(elem.ptr, elem.size, self->hash_seeds[i]);
            neighborhood[i] = temp % self->segmentSize + segmentStart;
        }
    }
    neighborhood[self->k] =
        murmur3_32(elem.ptr, elem.size, self->hash_seeds[self->k]) &
        (((uint64_t)1 << self->q) - 1);
}

/**
 * @brief hashes an element k + 1 times and stores the results in neighborhood
 * given the hash seeds in self, reusing the hashes in elemHash if they are
 * valid
 *
 * @param[in] self the fxlt
 * @param[in] elem the element to hash
 * @param[in] elemHash the cached hashes for elem
 * @param[out] neighborhood if self->cacheHashes is false, neighborhood should
 * point to an array of size k + 1. Otherwise, it does not need to point to
 * valid memory (it will be overwritten). However it is best practice for the
 * passed pointer to initially point to a stack-allocated array, to ensure
 * correctness whether or not caching is enabled.
 *
 * eg.
 * size_t hashes[self->k + 1];
 * size_t *neighborhood = hashes;
 * fuseXORierLT_hashAllCache(self, ... , ... , &neighborhood);
 *
 * whether or not caching is enabled, neighborhood will point to valid memory,
 * and hashes will not cause memory issues
 */
void fuseXORierLT_hashAllCache(struct fuseXORierLookupTable *self,
                               SizedPointer elem, struct itemHash *elemHash,
                               size_t **neighborhood) {
    if (self->cacheHashes == false) {
        fuseXORierLT_hashAll(self, elem, *neighborhood);
        return;
    }
    if (!elemHash->hashesValid) {
        fuseXORierLT_hashAll(self, elem, elemHash->hashes);
        elemHash->hashesValid = true;
    }
    *neighborhood = elemHash->hashes;
}

/**
 * @brief given a set of elements S, finds an ordering on S and a mapping from S
 * to [0, m) such that every element in S is mapped to a unique location in [0,
 * m) and the mapping is consistent with the ordering
 *
 * @param[in] self the fxlt
 * @param[in] S the set of keys to map
 * @param[in] SHashes the cached hashes for S
 * @param[in] SLen the length of S
 * @return struct mappingResult
 */
struct mappingResult fuseXORierLT_findMatch(struct fuseXORierLookupTable *self,
                                            SizedPointer *S,
                                            struct itemHash *SHashes,
                                            size_t SLen) {

    size_t printInterval = roundUp2(SLen / 10000);

    // make a hash table of the frequencies of each location
    struct xorSet *H = calloc(self->m, sizeof *H);
    size_t hashes[self->k + 1], *neighborhood = hashes;
    for (size_t i = 0; i < SLen; ++i) {
        if (self->verbose && (i & (printInterval - 1)) == 0)
            fprintf(
                stderr,
                "counting frequencies - |Counted|: %ld/%ld = %.1f%%"
                "                                                              "
                "  \r",
                i, SLen, 100.0 * i / SLen);

        SizedPointer t = S[i];
        struct itemHash *tHash = SHashes + i;

        fuseXORierLT_hashAllCache(self, t, tHash, &neighborhood);

        for (size_t l = 0; l < self->k; ++l) {
            size_t L = neighborhood[l];
            addElem(&H[L], i);
        }
    }

    struct mappingResult res;
    struct match *ordering = malloc(sizeof *ordering * SLen);
    size_t orderingInsert = SLen - 1;

    res.success = true;
    res.ordering = ordering;

    struct circularBuffer Q;
    construct_circularBuffer(&Q, 1);

    for (size_t i = 0; i < self->m; ++i) {
        if (singleton(&H[i])) {
            circularBuffer_enqueue(&Q, i);
        }
    }

    while (!circularBuffer_isEmpty(&Q)) {

        if (self->verbose &&
            ((SLen - 1 - orderingInsert) & (printInterval - 1)) == 0)
            fprintf(stderr,
                    "find match - |Singletons On Deck|: %-10ld |Matched|: "
                    "%ld/%ld = %.1f%%"
                    "                         "
                    "  \r",
                    circularBuffer_numElems(&Q), SLen - 1 - orderingInsert,
                    SLen, 100.0 * (SLen - 1 - orderingInsert) / SLen);

        size_t i = circularBuffer_dequeue(&Q);

        if (singleton(&H[i])) {
            size_t keyIdx = getSingleton(&H[i]);
            SizedPointer x = S[keyIdx];
            struct itemHash *xHash = SHashes + keyIdx;

            fuseXORierLT_hashAllCache(self, x, xHash, &neighborhood);
            int l = -1;
            for (size_t j = 0; j < self->k; ++j) {
                if (neighborhood[j] == i) {
                    l = j;
                }
                removeElem(&H[neighborhood[j]], keyIdx);
                if (singleton(&H[neighborhood[j]])) {
                    circularBuffer_enqueue(&Q, neighborhood[j]);
                }
            }
            if (l == -1) {
                fprintf(stderr, "l == -1 (how?)\n");
                exit(1);
            }

            // HashTable_insert(matching, x, (void *)(uintptr_t)l);
            ordering[orderingInsert--] = (struct match){keyIdx, l};
        }
    }

    if (orderingInsert != (size_t)-1) {
        free(ordering);
        res.success = false;
    }

    free(H);
    destroy_circularBuffer(&Q);

    return res;
}

size_t fuseXORierLT_findPlace(struct fuseXORierLookupTable *self,
                              SizedPointer t) {
    size_t hashes[self->k + 1];
    fuseXORierLT_hashAll(self, t, hashes);
    slot_t M = hashes[self->k];

    slot_t l = M;
    for (size_t j = 0; j < self->k; ++j) {
        l ^= self->table1[hashes[j]];
    }
    if (l < self->k) {
        size_t L = hashes[l];
        return L;
    }
    return (size_t)-1;
}

//------------------------------------------------------------------------------
// FXLT Interface Function Definitions
//------------------------------------------------------------------------------

struct fuseXORierLookupTable *build_fuseXORierLT(size_t n, SizedPointer keys[n],
                                                 void *values[n], float c,
                                                 size_t k, size_t flags) {
    struct fuseXORierLookupTable *self =
        malloc(sizeof(struct fuseXORierLookupTable));
    size_t m = (size_t)(c * (float)n);

    if (k != 3 && k != 4) {
        fprintf(stderr, "k must be 3 or 4\n");
        return NULL;
    }

    size_t segmentSize =
        (flags & FXLT_FLAG_NO_SPATIAL_COUPLING)
            ? 0
            : ((k == 3) ? 4.8 * pow(n, .58) : .7 * pow(n, .65));
    segmentSize = roundUp2(segmentSize);

    struct fuseXORierLookupTable init = {
        .k = k,
        .q = sizeof(slot_t) * 8,
        .segmentSize = segmentSize,
        .n = n,
        .m = (segmentSize == 0) ? m : segmentSize * (m / segmentSize + 1),
        .numSegments = (segmentSize == 0) ? 0 : m / segmentSize + 1,
        .table1 = malloc(sizeof *self->table1 * m),
        .table2 = malloc(sizeof *self->table2 * m),
        .verbose = flags & FXLT_FLAG_PRINT_STATS,
        .cacheHashes = flags & FXLT_FLAG_CACHE_HASHES,
    };
    memcpy(self, &init, sizeof init);

    self->hash_seeds = malloc(sizeof *self->hash_seeds * (k + 2));

    struct itemHash *itemHashes = NULL;
    size_t(*hashScratchpad)[self->k + 1] = NULL;
    if (self->cacheHashes) {
        itemHashes = malloc(sizeof *itemHashes * n);
        hashScratchpad = malloc(sizeof *hashScratchpad * (k + 1) * n);
        for (size_t i = 0; i < n; ++i) {
            itemHashes[i].hashes = hashScratchpad[i];
            itemHashes[i].hashesValid = false;
        }
    }

    struct mappingResult res;
    res.success = false;
    size_t tries = 0;
    for (;;) {
        if (((!(flags & FXLT_FLAG_ALLOW_RETRY_BUILD)) && tries != 0) ||  tries > 10  ) {
            delete_fuseXORierLT(self);
            return NULL;
        }

        for (size_t i = 0; i < k + 2; ++i) {
            self->hash_seeds[i] = random();
        }

        res = fuseXORierLT_findMatch(self, keys, itemHashes, n);
        ++tries;
        if (self->verbose)
            fprintf(stderr,
                    "                                                      "
                    "                                                      "
                    "                       \r");

        if (res.success) {
            break;
        }
        if (self->cacheHashes) {
            for (size_t i = 0; i < n; ++i) {
                itemHashes[i].hashesValid = false;
            }
        }
    }

    int printInterval = roundUp2(n / 10000);

    size_t hashes[self->k + 1];
    size_t *neighborhood = hashes;

    for (size_t i = 0; i < n; ++i) {
        if (self->verbose && ((i & (printInterval - 1)) == 0))
            fprintf(stderr,
                    "building - |Built|: %ld/%ld = %.1f%%                      "
                    "                                     \r",
                    i, n, 100.0 * i / n);
        struct match pair = res.ordering[i];
        SizedPointer key = keys[pair.keyIdx];
        struct itemHash *keyHash = itemHashes + pair.keyIdx;
        void *val = values[pair.keyIdx];
        slot_t l = pair.l;

        fuseXORierLT_hashAllCache(self, key, keyHash, &neighborhood);
        slot_t M = neighborhood[self->k];

        size_t L = neighborhood[l];

        self->table1[L] = M ^ l;
        for (size_t j = 0; j < self->k; ++j) {
            if (j == l) {
                continue;
            }
            self->table1[L] ^= self->table1[neighborhood[j]];
        }
        self->table2[L] = val;
    }

    if (self->cacheHashes) {
        free(hashScratchpad);
        free(itemHashes);
    }

    free(res.ordering);
    if (self->verbose)
        fprintf(stderr, "                                                      "
                        "                                         \r");
    return self;
}

void *fuseXORierLT_lookup(struct fuseXORierLookupTable *self,
                          SizedPointer elem) {
    size_t L = fuseXORierLT_findPlace(self, elem);
    if (L == (size_t)-1) {
        return NULL;
    }
    return self->table2[L];
}

bool fuseXORierLT_set(struct fuseXORierLookupTable *self, SizedPointer elem,
                      void *value) {
    size_t L = fuseXORierLT_findPlace(self, elem);
    if (L == (size_t)-1) {
        return false;
    }
    self->table2[L] = value;
    return true;
}

void delete_fuseXORierLT(struct fuseXORierLookupTable *self) {
    free(self->hash_seeds);
    free(self->table1);
    free(self->table2);
    free(self);
}

void fuseXORierLT_printStats(struct fuseXORierLookupTable *this) {
    if (this->cacheHashes) {
        size_t mem = sizeof *this + sizeof *this->hash_seeds * (this->k + 2) +
                     sizeof *this->table1 * this->m +
                     sizeof *this->table2 * this->m;
        size_t peakMem =
            mem + sizeof(struct itemHash) * this->n +
            sizeof(size_t[this->k + 1]) * this->n +
            sizeof(struct mappingResult) + sizeof(struct itemHash *) * this->n +
            sizeof(struct xorSet) * this->m + sizeof(size_t) * this->n / 5;

        fprintf(stderr, "--Cached Bloomier Filter--\n");
        fprintf(stderr, "k: %ld ", this->k);
        fprintf(stderr, "m: %ld ", this->m);
        fprintf(stderr, "q: %d ", this->q);
        fprintf(stderr, "segmentSize: %ld ", this->segmentSize);
        fprintf(stderr, "\n");
        fprintf(stderr, "mem:     %ld\n", mem);
        fprintf(stderr, "peakMem: %ld\n", peakMem);
        fprintf(stderr, "--------------------------\n");
    } else {
        size_t mem = sizeof *this + sizeof *this->hash_seeds * (this->k + 2) +
                     sizeof *this->table1 * this->m +
                     sizeof *this->table2 * this->m;
        size_t peakMem = mem + sizeof(struct mappingResult) +
                         sizeof(struct xorSet) * this->m +
                         sizeof(size_t) * this->n / 5;

        fprintf(stderr, "--Bloomier Filter--\n");
        fprintf(stderr, "k: %ld ", this->k);
        fprintf(stderr, "m: %ld ", this->m);
        fprintf(stderr, "q: %d ", this->q);
        fprintf(stderr, "segmentSize: %ld ", this->segmentSize);
        fprintf(stderr, "\n");
        fprintf(stderr, "mem:     %ld\n", mem);
        fprintf(stderr, "peakMem: %ld\n", peakMem);
        fprintf(stderr, "-------------------\n");
    }
}

struct footprint fuseXORierLT_getSize(struct fuseXORierLookupTable *self) {
    if (self->cacheHashes) {
        size_t mem = sizeof *self + sizeof *self->hash_seeds * (self->k + 2) +
                     sizeof *self->table1 * self->m +
                     sizeof *self->table2 * self->m;
        size_t peakMem =
            mem + sizeof(struct itemHash) * self->n +
            sizeof(size_t[self->k + 1]) * self->n +
            sizeof(struct mappingResult) + sizeof(struct itemHash *) * self->n +
            sizeof(struct xorSet) * self->m + sizeof(size_t) * self->n / 5;
        return (struct footprint){mem, peakMem};
    } else {
        size_t mem = sizeof *self + sizeof *self->hash_seeds * (self->k + 2) +
                     sizeof *self->table1 * self->m +
                     sizeof *self->table2 * self->m;
        size_t peakMem = mem + sizeof(struct mappingResult) +
                         sizeof(struct xorSet) * self->m +
                         sizeof(size_t) * self->n / 5;
        return (struct footprint){mem, peakMem};
    }
}

//------------------------------------------------------------------------------
// Helper Function Definitions
//------------------------------------------------------------------------------

void addElem(struct xorSet *set, size_t elem) {
    set->elems ^= elem;
    ++set->cardinality;
}
void removeElem(struct xorSet *set, size_t elem) {
    set->elems ^= elem;
    --set->cardinality;
}
bool singleton(struct xorSet *set) { return set->cardinality == 1; }
size_t getSingleton(struct xorSet *set) { return set->elems; }

void construct_circularBuffer(struct circularBuffer *self, size_t size) {
    size = roundUp2(size);
    self->elems = malloc(sizeof *self->elems * size);
    self->size = size;
    self->insert = 0;
    self->numElems = 0;
}

void destroy_circularBuffer(struct circularBuffer *self) { free(self->elems); }

bool circularBuffer_isEmpty(struct circularBuffer *self) {
    return self->numElems == 0;
}

size_t circularBuffer_numElems(struct circularBuffer *self) {
    return self->numElems;
}

void circularBuffer_grow(struct circularBuffer *self) {
    size_t *newElems = malloc(sizeof *newElems * self->size * 2);
    for (size_t i = 0; i < self->numElems; ++i) {
        newElems[i] = self->elems[(self->insert + i) & (self->size - 1)];
    }
    free(self->elems);
    self->elems = newElems;
    self->insert = self->numElems;
    self->size *= 2;
}

void circularBuffer_enqueue(struct circularBuffer *self, size_t elem) {
    if (self->numElems == self->size) {
        circularBuffer_grow(self);
    }
    self->elems[(self->insert) & (self->size - 1)] = elem;
    self->insert = (self->insert + 1) & (self->size - 1);
    ++self->numElems;
}

size_t circularBuffer_dequeue(struct circularBuffer *self) {
    size_t removeFrom =
        (self->insert - self->numElems + self->size) & (self->size - 1);
    size_t elem = self->elems[removeFrom];
    --self->numElems;
    return elem;
}


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
        // Here is a source of differing results across endianness.
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