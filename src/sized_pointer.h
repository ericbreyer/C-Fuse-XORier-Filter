/**
 * @file sized_pointer.h
 * @author Eric Breyer (eric.breyer@gmail.com) [https://github.com/ericbreyer]
 * @brief Declaration of a sized pointer type. Used to be able to hash arbitrary data.
 * @version 1.0
 * @date 2023-12-13
 *
 * @copyright Copyright (c) 2023, released under GNU Public Licence 3.0 or later
 *
 */
// SPDX-License-Identifier: GNU-3.0-or-later

#pragma once
#include <stdint.h>

typedef struct SizedPointer_s {
    void * ptr;
    size_t size;
} SizedPointer;