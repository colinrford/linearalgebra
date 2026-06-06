# `lam.linearalgebra`

After taking some time away from battling these dragons, I have returned to
this project: `lam.linearalgebra` is now a c++ module and part of
[LAM](https://www.github.com/colinrford/lam). It is currently in a state of flux.

## linearalgebra, a [LAM](https://www.github.com/colinrford/lam) library

First things first, a quick disclaimer, if you will:
There is no beating `BLAS` or `LAPACK` it would seem, so I have decided to offer the user the choice to use `BLAS` or `LAPACK` via a `cmake` `LAM_USE_BLAS` option, which `cmake` will try and set to `ON` if it can find the `BLAS` or `LAPACK` libraries.

## building

What you'll need:
- [`lam.concepts`](https://www.github.com/colinrford/concepts)
- `cmake`
- `Ninja`
- compiler with support for `import std;`

Build along the lines of:
```bash
git clone https://www.github.com/colinrford/linearalgebra
cd linearalgebra
cmake -S . -B build -G Ninja
cmake --build build
```

Alternatively, consume with `conan` using an `import std`-friendly profile:
```bash
conan create . --profile <my-import-std-profile>
```

## usage
For example, see [`lam.lebesgue`](https://www.github.com/colinrford/lebesgue).
