# SPDX-License-Identifier: GNU-3.0-or-later

target_include_directories(fuse_XORier_lookup_table PUBLIC .)
target_sources(fuse_XORier_lookup_table PRIVATE fuse_XORier_lookup_table.c fuse_XORier_lookup_table.h)
set_target_properties(fuse_XORier_lookup_table PROPERTIES PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/fuse_XORier_lookup_table.h;${CMAKE_CURRENT_SOURCE_DIR}/sized_pointer.h")
