#pragma once
#include <stdint.h>

typedef struct SizedPointer_s {
    void * ptr;
    size_t size;
} SizedPointer;