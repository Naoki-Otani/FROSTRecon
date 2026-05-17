# FixBMpln

Temporary tools to modify `BMBasicRecon` content in ROOT files.

## Apps

This project currently provides two apps:

- `fix_bmbasicrecon_class`
- `remove_deadpln_topview`

## 1. `fix_bmbasicrecon_class`

Reads all `.root` files in an input directory and fixes entries in `BMBasicRecon` where:

- `mod == 5`
- `view == 0`

For those elements, it applies:

```cpp
pln -= 1;
```

Output file name:

- `input.root` -> `input_fixed.root`

## 2. `remove_deadpln_topview`

Reads all `.root` files in an input directory and removes entries in `BMBasicRecon` where:

- `mod == 5`
- `view == 1`
- `pln == 9`, `10`, or `17`

This is used to intentionally create dead planes in Baby MIND top-view.

Output file name:

- `input.root` -> `input.root`

## What they copy

For each input file, both apps read and write:

- `BMBasicRecon`
- `BMBeaminfo`
- `BMBSD`

Only `BMBasicRecon` content is modified.

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
    ├── fix_bmbasicrecon_class.cpp
    └── remove_deadpln_topview.cpp
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

### Fix wrong `pln`

```bash
./fix_bmbasicrecon_class INPUT_DIR OUTPUT_DIR
```

Example:

```bash
./fix_bmbasicrecon_class \
  /group/nu/ninja/work/otani/FROSTReconData/BMdata/2-BMBSD \
  /group/nu/ninja/work/otani/FROSTReconData/BMdata/2-BMBSD_fixed
```

### Remove top-view dead planes

```bash
./remove_deadpln_topview INPUT_DIR OUTPUT_DIR
```

Example:

```bash
./remove_deadpln_topview \
  /group/nu/ninja/work/otani/FROSTReconData/BMdata/2-BMBSD \
  /group/nu/ninja/work/otani/FROSTReconData/BMdata/2-BMBSD_deadpln
```

## Root cause of the original `pln` bug

The original problem came from the new `PickSignal` code storing BM horizontal hits with:

```cpp
pln.push_back(imod + 1);
```

instead of:

```cpp
pln.push_back(imod);
```

`fix_bmbasicrecon_class` is only a temporary fix for already-produced ROOT files.
