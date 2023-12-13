//------------------------------------------------------------------------------
// Information
//------------------------------------------------------------------------------

/**
 * @file fuse_xorier_lookup_table.h
 * @author Eric Breyer (eric.breyer@gmail.com) [https://github.com/ericbreyer]
 * @brief Declaration of a fuse XORier lookup table
 * @version 1.0
 * @date 2023-12-12
 *
 * @copyright Copyright (c) 2023, released under GNU Public Licence 3.0 or later
 *
 */
// SPDX-License-Identifier: GNU-3.0-or-later

//------------------------------------------------------------------------------
// Includes
//------------------------------------------------------------------------------

#include "sized_pointer.h" // SizedPointer
#include <stdbool.h>       // bool
#include <stddef.h>        // size_t

#pragma once

//------------------------------------------------------------------------------
// Types
//------------------------------------------------------------------------------

/**
 * @brief Encapsulation of a fuse XORier lookup table
 *
 */
struct fuseXORierLookupTable;

/**
 * @brief Pair type for fuse XORier lookup table statistics
 *
 */
struct footprint {
    size_t size;
    size_t peak_size;
};

/**
 * @brief Build option flags for fuse XORier lookup table construction
 *
 */
enum fuseXORierLTFlags {
    FXLT_FLAG_PRINT_STATS = 1 << 0,
    FXLT_FLAG_CACHE_HASHES = 1 << 1,
    FXLT_FLAG_NO_SPATIAL_COUPLING = 1 << 2,
    FXLT_FLAG_ALLOW_RETRY_BUILD = 1 << 3,
};

//------------------------------------------------------------------------------
// Functions
//------------------------------------------------------------------------------

/**
 * @brief Retrieve a value corrsponding to a key from a fused XORier lookup
 * table. This function has a small probability of returning a false positive,
 * ie. returning a value for a key that was not inserted. If the key is not in
 * the table, the function will return NULL. As such, NULL must not be used as a
 * meaningful value in the table. There is no probability of a false negative,
 * ie. if the function returns NULL, the key is not in the table.
 *
 * @param[in] self
 * @param[in] elem
 * @return void*
 */
void *fuseXORierLT_lookup(struct fuseXORierLookupTable *self,
                          SizedPointer elem);

/**
 * @brief Change the value associated with a key in a fuse XORier lookup table.
 * NOTE: The key must have been inserted into the table in the construction. If
 * the key is not in the table, but the lookup is a false positive, the set
 * could appear to succeed while overwriting a different key's value.
 *
 * @param[in] self The fuse XORier lookup table
 * @param[in] elem The key
 * @param[in] value The new value
 * @return bool Whether the set (seemingly) succeeded
 */
bool fuseXORierLT_set(struct fuseXORierLookupTable *self, SizedPointer elem,
                      void *value);

/**
 * @brief Construct a fuse XORier lookup table
 *
 * @param[in] n The number of keys and values to insert
 * @param[in] keys The keys to insert
 * @param[in] values The values to insert
 * @param[in] c The load factor of the table
 * @param[in] k The number of XORier hash functions to use
 * @param[in] flags Bitmask of flags to control the construction (multiple flags
 * can be combined with bitwise OR)
 * @return The constructed table
 */
struct fuseXORierLookupTable *build_fuseXORierLT(size_t n, SizedPointer keys[n],
                                                 void *values[n], float c,
                                                 size_t k, size_t flags);

/**
 * @brief Deconstruct the fuse XORier lookup table and free its memory
 *
 * @param[in] self the Fuse XORier lookup table
 */
void delete_fuseXORierLT(struct fuseXORierLookupTable *self);

//------------------------------------------------------------------------------
// Debugging functions
//------------------------------------------------------------------------------

/**
 * @brief get the size and peak footprint of a fuse XORier lookup table
 *
 * @param[in] self The fuse XORier lookup table
 * @return struct footprint
 */
struct footprint fuseXORierLT_getSize(struct fuseXORierLookupTable *self);

/**
 * @brief Debugging function to print the statistics of a fuse XORier lookup
 * table
 *
 * @param[in] self The fuse XORier lookup table
 */
void fuseXORierLT_printStats(struct fuseXORierLookupTable *self);