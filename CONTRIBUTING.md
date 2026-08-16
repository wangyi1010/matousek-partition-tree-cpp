# Contributing

## Setup

Requirements: CMake 3.20+, a C++20 compiler, and Boost headers.

```bash
git clone https://github.com/wangyi1010/matousek-partition-tree-cpp
cd matousek-partition-tree-cpp
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

## Workflow

Before pushing:

```bash
cmake --build build --parallel
ctest --test-dir build --output-on-failure
./build/matousek-demo 400 42
./build/measure-performance --smoke
```

## Conventions

- Use C++20 and keep the public API in `include/matousek_partition_tree/core.hpp`.
- All geometric coordinates and predicates use `boost::rational<boost::multiprecision::cpp_int>`.
- Floating point is allowed only for randomized cutting control parameters, SVG coordinates, and timing output.
- Verified construction failures throw `CuttingError` or `TestSetError`; never weaken a runtime check silently.
- Avoid mutable global state. Pass `std::mt19937_64` explicitly when deterministic random-state continuity matters.
- Add property tests for behavioral changes.
