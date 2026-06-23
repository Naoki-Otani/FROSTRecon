# FROST-Baby MIND Matching Tools

This package contains two executables:

1. **`HitConverter`**: aligns FROST events with Baby MIND/WAGASCI spills and adds FROST-related trees to a B2 ROOT file. This is used for real data.
2. **`TrackMatch`**: reconstructs Baby MIND track quantities, matches them to FROST hits, and writes analysis summary trees. It supports both real data and MC input.

Typical real-data workflow:

```text
B2/Baby MIND file + FROST recon ROOT directory
        -> HitConverter
        -> ROOT file with frost_match and match_info trees
        -> TrackMatch in data mode
        -> ROOT file with ntbm, frost_match, and enriched match_info
```

Typical MC workflow:

```text
FROSTRecon MC ROOT file
  containing tree, frost, frostmc, and optionally nRooTracker
        -> TrackMatch in MC mode
        -> ROOT file with ntbm, frost_match, match_info, frostmc, and nRooTracker
```

For MC, `HitConverter` is not used. The FROSTRecon MC file already contains the WAGASCI/BabyMIND B2 tree and the reconstructed FROST tree:

```text
tree + frost + frostmc + nRooTracker
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

`TrackMatch` reconstructs Baby MIND track position/tangent information, extrapolates each Baby MIND track to the FROST plane, and stores BM-FROST matching information.

The input tree layout depends on the input mode:

- **Data mode**: reads `frost_match` and `match_info` from the `HitConverter` output.
- **MC mode**: reads `frost` directly from the FROSTRecon MC output. It also reads `frostmc` to fill FROST truth-particle information. `nRooTracker` is copied to the output if present.

It writes:

- `ntbm`: a tree with one `NTBMSummary` object per spill.
- updated `match_info`: spill-level metadata plus track-wise vectors.
- updated `frost_match`: copied from the input FROST tree.

For data input, `frost_match` is cloned from the input `frost_match` tree. For MC input, the input `frost` tree is cloned and written as `frost_match` so that downstream jobs can use a common output tree name.

For MC input, `frostmc` and `nRooTracker` are also copied to the output if they exist in the input file.

### Usage

```bash
TrackMatch <input B2 ROOT file> <output ROOT file> <z shift> <MC(0)/data(1)> [fhc|rhc] [trace|debug|info|warning|error|fatal]
```

Data example:

```bash
bin/TrackMatch/TrackMatch \
  afterHitConverter.root \
  afterTrackMatch.root \
  0.0 \
  1 \
  info
```
MC example:

```bash
bin/TrackMatch/TrackMatch \
  frost_recon_mc.root \
  afterTrackMatchMC.root \
  0.0 \
  0 \
  fhc \
  info
```

Use `0` for MC and `1` for real data.

For MC input, the beam mode must be specified explicitly with `fhc` or `rhc`. The MC input currently does not provide a meaningful `BsdGoodSpillFlag`, so TrackMatch does not infer the beam mode from that flag in MC mode.

For real data, the optional `fhc`/`rhc` argument is ignored. The beam mode is determined from `BsdGoodSpillFlag` in the B2 beam summary.

### Matching rules

Current matching behavior is:

- Position allowance cuts are not used as the BM-FROST matching criterion.
- A matched track requires both x-view and y-view FROST reconstructed candidates to exist.
- For real data, FROST reconstructed x/y positions are scaled by constants in `NTBMConst.hh`; MC positions are not scaled.
- Charge handling uses the beam mode. For real data, the beam mode is read from `BsdGoodSpillFlag`. For MC, the beam mode is taken from the explicit `fhc` or `rhc` command-line argument.

Data matching:

- Baby MIND bunch is `1..8`; FROST arrays are indexed as `0..7`.
- A FROST hit is accepted only when `is_hit[bunch - 1] == 1`.
- The matching candidates are taken from `x_rec[bunch - 1]` and `y_rec[bunch - 1]`.
- A matched track requires both `x_rec[bunch - 1]` and `y_rec[bunch - 1]` to contain at least one candidate.

MC matching:

- MC has no physical bunch structure in the TrackMatch input.
- FROSTRecon stores the meaningful MC hit information in index `0`.
- A FROST hit is accepted only when `is_hit[0] == 1`.
- The matching candidates are taken from `x_rec[0]` and `y_rec[0]`.
- A matched track requires both `x_rec[0]` and `y_rec[0]` to contain at least one candidate.
- The MC charge assignment uses the explicitly specified beam mode, `fhc` or `rhc`, rather than `BsdGoodSpillFlag`.

For both data and MC:

- `trackmatch_expected_x/y` is the Baby MIND track extrapolated to the FROST plane.
- `trackmatch_frost_nearest_x/y` is the nearest reconstructed FROST position to the Baby MIND extrapolated position.
- `trackmatch_dx/dy` is `FROST position - expected position`.
- `trackmatch_tangent_x/y` is the tangent from the Baby MIND second-layer position to the selected FROST position.
- `trackmatch_dtanx/dtany` is the difference between the BM-FROST tangent and the Baby MIND fitted tangent.
- Position allowance cuts are not used as the BM-FROST matching criterion.

Beam-mode handling:

- Data: positive `BsdGoodSpillFlag` is treated as FHC, and negative `BsdGoodSpillFlag` is treated as RHC.
- MC: `BsdGoodSpillFlag` is not used for beam-mode determination. The fifth argument after `<MC(0)/data(1)>` must be `fhc` or `rhc`.
- Running MC without `fhc` or `rhc` is treated as an invalid command-line argument.

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
- MC normalization information when running on MC
- FROST truth-particle information when running on MC

The FROST truth-particle information is stored as spill-level vectors:

- `true_frost_particle_id_`: PDG particle id of particles that hit FROST.
- `true_frost_position_x_`: x position at `z = 0 mm` in FROST local coordinates.
- `true_frost_position_y_`: y position at `z = 0 mm` in FROST local coordinates.
- `true_frost_tangent_x_`: x-z tangent, computed as `px / pz`.
- `true_frost_tangent_y_`: y-z tangent, computed as `py / pz`.

The same vector index corresponds to the same truth particle in all five vectors. The number of FROST truth particles is independent of `number_of_frost_match_entries_`.

The truth positions are calculated from `frostmc` using `particle_local_first_*_mm` and `particle_local_last_*_mm`, linearly extrapolated to `z = 0 mm`. The truth tangents are calculated from `particle_initial_px_mev_c / particle_initial_pz_mev_c` and `particle_initial_py_mev_c / particle_initial_pz_mev_c`.

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
- `trackmatch_true_frost_nearest_particle_id`
- `trackmatch_true_frost_nearest_position_x`
- `trackmatch_true_frost_nearest_position_y`
- `trackmatch_true_frost_nearest_tangent_x`
- `trackmatch_true_frost_nearest_tangent_y`

Definitions used in `match_info`:

- `trackmatch_expected_x/y`: Baby MIND track extrapolated to the FROST plane.
- `trackmatch_frost_nearest_x/y`: nearest FROST reconstructed position selected for compact output.
- `trackmatch_dx/dy`: `FROST position - expected position`.
- `trackmatch_tangent_x/y`: tangent from the Baby MIND second-layer position to the FROST position.
- `trackmatch_dtanx/dtany`: difference between BM-FROST tangent and Baby MIND fitted tangent.
- `trackmatch_frost_is_hit`:
  - Data: `is_hit[bunch - 1]`
  - MC: `is_hit[0]`
  - stored in the same row order as the other `trackmatch_*` vectors.
- `bsd_good_spill_flag`: copied from the B2 beam summary. In MC this may remain non-initialized and should not be used to infer FHC/RHC.
- `trackmatch_true_frost_nearest_particle_id`: PDG particle id of the nearest FROST truth particle.
- `trackmatch_true_frost_nearest_position_x/y`: nearest truth-particle position at `z = 0 mm` in FROST local coordinates.
- `trackmatch_true_frost_nearest_tangent_x/y`: nearest truth-particle tangent, computed as `px / pz` and `py / pz`.

For the `trackmatch_true_frost_nearest_*` branches, TrackMatch selects the truth particle closest to the Baby MIND extrapolated position. If one or more FROST-hit muons are present, only muons with PDG id `+13` or `-13` are considered. If no FROST-hit muon is present, all FROST-hit truth particles are considered.

For data, the `trackmatch_true_frost_nearest_*` values are filled with non-initialized defaults.

## Developer notes

- `lib/NTBMSummary.*` defines the custom summary object and its ROOT dictionary.
- `src/HitConverter/HitConverterDict.hh` and `HitConverterLinkDef.h` provide ROOT dictionaries for nested STL containers used by FROST branches.
- `NTBMConst.hh` collects matching constants, geometry constants, FROST scale factors, and common non-initialized values.
- The project rejects in-tree builds; always configure in a separate `build/` directory.
