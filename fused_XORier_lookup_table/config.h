/**
 * @file config.h
 * @author Eric Breyer (eric.breyer@gmail.com)
 * @brief Configuration file for fuse XORier lookup table
 * @version 1.0
 * @date 2023-12-12
 *
 * @copyright Copyright (c) 2023, released under GNU Public Licence 3.0 or later
 *
 */
// SPDX-License-Identifier: GNU-3.0-or-later

#pragma once

/**
 * @brief The type of slot used in the lookup table. False positive rate is
 * equal to 2^(-sizeof(slot_t))
 *
 */
typedef unsigned char slot_t;