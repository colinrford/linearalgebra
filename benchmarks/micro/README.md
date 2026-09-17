# benchmarks/micro

Slide-grade throughput benchmarks built on [google/benchmark](https://github.com/google/benchmark)
(warmup, auto-iteration to statistical stability, mean/median/stddev). This is the
only external dependency in `lam.linearalgebra` outside the BLAS backends.

## Why it's isolated here

LAM is deliberately minimal-dependency — `lam.concepts` and the standard library
only. google/benchmark is therefore **optional and quarantined**:

- The parent `benchmarks/CMakeLists.txt` gates this folder behind
  `find_package(benchmark QUIET)`. Not installed → `micro/` is silently skipped;
  everything else still builds. (We use `find_package`, **not** `FetchContent`,
  which would always pull the dep and defeat the point.)
- Nothing here touches the library or its install/export surface.

## Layout

- `bench_support.cppm` — module `linalg.bench.support`. Reusable, **dependency-clean**
  scaffolding (allocation-profile defs, deterministic fillers, the RAII monotonic
  arena). Imports `std` + `lam.linearalgebra` only — **no** google/benchmark. The
  `DoNotOptimize`/`ClobberMemory` barrier deliberately lives in the `.cpp` TUs
  below, where `<benchmark/benchmark.h>` is already included.
- `bm_*.cpp` — the benchmark translation units (`linalg_gb_*` targets). Each
  `#include <benchmark/benchmark.h>` and imports the support module.

## Conventions

- Targets are prefixed `linalg_gb_` — the `gb` flags "needs google/benchmark".
- Binaries build into `bin/benchmarks/`, not the flat `bin/` (deliberate).

## Running

```sh
brew install google-benchmark      # or your platform's package
# configure/build the linearalgebra dev build, then:
./bin/benchmarks/linalg_gb_smoke --benchmark_out=results.json --benchmark_out_format=json
```
