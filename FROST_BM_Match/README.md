# FROST-Baby MIND Matching Tools

This package contains two executables used in sequence:

1. **`HitConverter`**: aligns FROST events with Baby MIND/WAGASCI spills and adds FROST-related trees to a B2 ROOT file.
2. **`TrackMatch`**: reconstructs Baby MIND track quantities, matches them to FROST hits, and writes analysis summary trees.

Typical workflow:

```text
B2/Baby MIND file + FROST recon ROOT directory
        -> HitConverter
        -> ROOT file with frost_match and match_info
        -> TrackMatch
        -> ROOT file with ntbm and enriched match_info
```

## Directory structure

```text
.
├── CMakeLists.txt
├── lib/
│   ├── CMakeLists.txt
│   ├── NTBMConst.hh
│   ├── NTBMLinkDef.hh
│   ├── NTBMSummary.hh
│   └── NTBMSummary.cc
├── src/
│   ├── HitConverter/
│   │   ├── CMakeLists.txt
│   │   ├── HitConverter.cc
│   │   ├── HitConverterDict.hh
│   │   └── HitConverterLinkDef.h
│   └── TrackMatch/
│       ├── CMakeLists.txt
│       ├── TrackMatch.cpp
│       └── TrackMatch.hpp
└── cmake/modules/
    ├── FindB2MC.cmake
    ├── FindEventDisplay.cmake
    └── FindNTBM.cmake
```

The files under `cmake/modules/` are custom CMake find modules. They are source files, not generated build products. Generated files should appear under your build directory.

## Build Environment

This package is intended to be built inside the Apptainer environment used by the NINJA and WAGASCI groups. The example build commands below assume that ROOT, Geant4, Boost, and the WAGASCI-Baby MIND/B2MC software are available in that environment.

Native builds outside the Apptainer environment may also work, but all required dependencies must be installed and discoverable by CMake.

## Requirements

The following dependencies are expected to be available in the NINJA/WAGASCI Apptainer environment:

- CMake
- C++17 compiler
- ROOT 6.18 or newer
- Geant4
- Boost, including `system`, `filesystem`, `program_options`, `unit_test_framework`, `log_setup`, and `log`
- B2MC / WAGASCI-Baby MIND libraries

## Build and install

```bash
# Set this to the installation directory of wagasci-babymind-monte-carlo.
export B2MCDIR=/path/to/wagasci-babymind-monte-carlo

# Set this to the installation directory of wagasci-babymind-event-display
# if your CMake configuration requires it.
export EVENT_DISPLAY_DIR=/path/to/wagasci-babymind-event-display

mkdir -p build
cd build

cmake .. \
  -DB2MC_PATH="${B2MCDIR}" \
  -DB2MC_INCLUDE_DIR="${B2MCDIR}/include/wagasci/b2" \
  -DB2MC_LIBRARY="${B2MCDIR}/lib/wagasci/b2/libB2MC.so" \
  -DEVENT_DISPLAY_PATH="${EVENT_DISPLAY_DIR}" \
  -DCMAKE_INSTALL_PREFIX=..

make -j"$(nproc)"
make install
```

Installed executables are placed under:

```text
bin/HitConverter/HitConverter
bin/TrackMatch/TrackMatch
```

`libNTBM.so`, ROOT dictionary files, and public headers are installed under the configured install prefix, mainly in `lib/ninja/recon` and `include/ninja/recon`.

## HitConverter

### What it does

`HitConverter` reads a B2/Baby MIND file and a directory of FROST reconstructed ROOT files. It matches B2 spills to FROST events using spill number and timestamp information, then writes an output B2 ROOT file with two additional trees:

- `frost_match`: FROST event data aligned one-to-one with B2/Baby MIND spills.
- `match_info`: spill-level metadata describing the B2-FROST match.

If a FROST event is missing, an empty `frost_match` entry is still written so that the row order remains aligned with `match_info` and the B2 spill entries.

### FROST input file names

FROST ROOT files are collected from one directory and sorted by run/event range. Supported patterns are:

```text
run00006_recon.root
run00006_0_9999_recon.root
```

### Usage

```bash
HitConverter <input B2/BabyMIND file> <input FROST directory> <output ROOT file> [trace|debug|info|warning|error|fatal]
```

Example:

```bash
bin/HitConverter/HitConverter \
  input_bm.root \
  /path/to/frost_recon_root_files \
  afterHitConverter.root \
  info
```

Log verbosity is ordered as:

```text
trace < debug < info < warning < error < fatal
```

`fatal` gives the least output.

### HitConverter output

`frost_match` stores the full FROST event schema, including raw arrays and reconstructed branches such as `is_hit`, `x_rec`, `y_rec`, `chi2`, `x_rectype`, `y_rectype`, `groupx`, and `groupy`.

`match_info` stores one entry per B2 spill, including:

- `bm_index`
- `wagasci_spill`, `wagasci_unixtime`
- `matched`
- `frost_file`, `frost_entry`, `frost_spillnum`, `frost_unixtime`
- `frost_run_number`, `frost_start_event`, `frost_end_event`
- `time_diff`

## TrackMatch

### What it does

`TrackMatch` reads the `HitConverter` output. It reconstructs Baby MIND track position/tangent information, extrapolates each Baby MIND track to the FROST plane, and stores BM-FROST matching information.

It writes:

- `ntbm`: a tree with one `NTBMSummary` object per spill.
- updated `match_info`: the original spill-level matching information plus track-wise vectors.
- updated `frost_match`: copied from the input file.

### Usage

```bash
TrackMatch <input B2 ROOT file> <output ROOT file> <z shift> <MC(0)/data(1)> [trace|debug|info|warning|error|fatal]
```

Example:

```bash
bin/TrackMatch/TrackMatch \
  afterHitConverter.root \
  afterTrackMatch.root \
  0.0 \
  1 \
  info
```

Use `0` for MC and `1` for real data.

### Matching rules

Current matching behavior is:

- Baby MIND bunch is `1..8`; FROST arrays are indexed as `0..7`.
- A FROST hit is accepted only when `is_hit[bunch - 1] == 1`.
- Position allowance cuts are not used as the BM-FROST matching criterion.
- A matched track requires both x-view and y-view FROST reconstructed candidates to exist.
- For real data, FROST reconstructed x/y positions are scaled by constants in `NTBMConst.hh`; MC positions are not scaled.
- Charge handling uses the sign of `BsdGoodSpillFlag`: positive means FHC and negative means RHC.

### Main output content

#### `ntbm`

The `ntbm` tree contains one branch:

```text
NTBMSummary
```

`NTBMSummary` stores:

- file and beam information
- Baby MIND track quantities: track type, momentum, position, tangent, maximum plane, track length, likelihoods, charge, bunch
- BM-FROST candidate information: candidate counts, FROST x/y candidates, expected positions, differences, and candidate tangents
- MC truth information when running on MC

#### updated `match_info`

`match_info` remains one entry per spill. Track-wise values are stored as vectors, with one vector element per Baby MIND track in that spill.

Important added branches:

- `trackmatch_has_match`
- `trackmatch_bm_track_id`
- `trackmatch_frost_match_bunch`
- `trackmatch_ninja_track_type`
- `trackmatch_baby_mind_tangent_x`, `trackmatch_baby_mind_tangent_y`
- `trackmatch_expected_x`, `trackmatch_expected_y`
- `trackmatch_frost_nearest_x`, `trackmatch_frost_nearest_y`
- `trackmatch_dx`, `trackmatch_dy`
- `trackmatch_tangent_x`, `trackmatch_tangent_y`
- `trackmatch_dtanx`, `trackmatch_dtany`
- `trackmatch_frost_is_hit`
- beam analysis branches: `spill_pot`, `bunch_pot`, `bsd_spill_number`, `timestamp`, `bsd_good_spill_flag`, `wagasci_good_spill_flag`, `detector_flags`

Definitions used in `match_info`:

- `trackmatch_expected_x/y`: Baby MIND track extrapolated to the FROST plane.
- `trackmatch_frost_nearest_x/y`: nearest FROST reconstructed position selected for compact output.
- `trackmatch_dx/dy`: `FROST position - expected position`.
- `trackmatch_tangent_x/y`: tangent from the Baby MIND second-layer position to the FROST position.
- `trackmatch_dtanx/dtany`: difference between BM-FROST tangent and Baby MIND fitted tangent.
- `trackmatch_frost_is_hit`: `is_hit[bunch - 1]`, stored in the same row order as the other `trackmatch_*` vectors.

## Developer notes

- `lib/NTBMSummary.*` defines the custom summary object and its ROOT dictionary.
- `src/HitConverter/HitConverterDict.hh` and `HitConverterLinkDef.h` provide ROOT dictionaries for nested STL containers used by FROST branches.
- `NTBMConst.hh` collects matching constants, geometry constants, FROST scale factors, and common non-initialized values.
- The project rejects in-tree builds; always configure in a separate `build/` directory.
