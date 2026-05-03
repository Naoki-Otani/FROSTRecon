# FixBMpln

A temporary tool to fix incorrect `pln` values in ROOT files.

## Purpose

This program reads all `.root` files in an input directory and fixes entries in `BMBasicRecon` where:

- `mod == 5`
- `view == 0`

For those elements, it applies:

```cpp
pln -= 1;
```

The output file name is:

- `input.root` -> `input_fixed.root`

## What it copies

For each input file, the program reads and writes:

- `BMBasicRecon`
- `BMBeaminfo`
- `BMBSD`

Only `BMBasicRecon::pln` is fixed.

## Directory layout

```text
FixBMpln/
├── CMakeLists.txt
├── include/
│   ├── BMConst.hpp
│   └── project/
│       ├── BMBasicRecon.hpp
│       ├── BMBeaminfo.hpp
│       ├── BMBSD.hpp
│       └── PMRecon.hpp
├── src/
│   ├── BMBasicRecon.cpp
│   ├── BMBeaminfo.cpp
│   ├── BMBSD.cpp
│   ├── PMRecon.cpp
│   ├── BMBasicReconLinkDef.hh
│   ├── BMBeaminfoLinkDef.hh
│   ├── BMBSDLinkDef.hh
│   └── PMReconLinkDef.hh
└── app/
    └── fix_bmbasicrecon_class.cpp
```

## Notes

- `BMConst.hpp` should be placed in `include/BMConst.hpp`.
- `BMBSD.hpp` needs `#include <cstdint>`.
- For this temporary build, `spdlog` usage in `BMConst.hpp` was commented out.

## Build

```bash
mkdir -p build
cd build
cmake ..
make -j"$(nproc)"
```

## Run

```bash
./fix_bmbasicrecon_class INPUT_DIR OUTPUT_DIR
```

Example:

```bash
./fix_bmbasicrecon_class \
  /group/nu/ninja/work/otani/FROSTReconData/BMdata/2-BMBSD \
  /group/nu/ninja/work/otani/FROSTReconData/BMdata/2-BMBSD_fixed
```

## Root cause

The original problem came from the new `PickSignal` code (PickSignal.cpp) storing BM horizontal hits with:

```cpp
pln.push_back(imod + 1);
```

instead of:

```cpp
pln.push_back(imod);
```

This tool is only a temporary fix for already-produced ROOT files.
