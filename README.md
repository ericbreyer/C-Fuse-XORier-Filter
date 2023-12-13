# C Fuse XORier Lookup Table

C Fuse XORier Lookup Table is a C library that implements the [Fuse XORier Lookup Table](about:blank). 
An FXLT is a static, probabilistic, associate array (map) in sublinear space. Similar to [Bloom](https://en.wikipedia.org/wiki/Bloom_filter) and [XOR filters](https://arxiv.org/pdf/1912.08258.pdf), an extension of the [Bloomier Filter](https://www.cs.princeton.edu/~chazelle/pubs/soda-rev04.pdf).

## Installation

### Build and Install From Source

```bash
git clone https://github.com/ericbreyer/C-Fuse-XORier-Filter
cd C-Fuse-XORier-Filter
mkdir build
cd build
cmake -B . -S ..
make fxlt_install
```

### Run Example/Graphs
```bash
make main
./test/main
```

## Usage

```c
#include <fuse_XORier_lookup_table.h>

int n = 10000 /* the number of keys to be inserted */;
SizedPointer keys[n] = /* <the keys> */;
void *values[n] = /* <the values the keys map to> */;
float c = 1.075 /* the constant multiple of n which is the size of the filter */;
size_t k = 4 /* the number of hash functions to use */;

struct fuseXORierLookupTable * FXLT = build_fuseXORierLT(n, keys, values, c, k, 0);
void * val = fuseXORierLT_lookup(FXLT, /* SizedPointer elem */);
if(val == NULL) {
    // elem is definitely not in the map
} else {
    // elem is vary likely in the map with value val
}
```

see [fuse_XORier_lookup_table.h](https://github.com/ericbreyer/C-Fuse-XORier-Filter/blob/main/src/fused_XORier_lookup_table.h) and [main.c](https://github.com/ericbreyer/C-Fuse-XORier-Filter/blob/main/test/main.c) for more examples and functionality.

## Contributing

Pull requests are welcome. Feel free to fork this repo and make a pull request when ready. Don't forget to add a star!

## Authors

- [Eric Breyer](https://github.com/ericbreyer) - code and research
- [Alan Liu](https://github.com/alanyliu) - research (python code in separate repo)

## License

[GNU 3.0 or Later](https://www.gnu.org/licenses/gpl-3.0.html#license-text)