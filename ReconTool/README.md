# FROST ReconTools

This package provides two command-line tools for generating mapping functions and running event-by-event reconstruction.

- `mapfunction_tool` generates mapping graphs from MC samples.
- `FROST_reconstruction` reconstructs hit positions event-by-event using the mapping graphs.

The code is refactored so that:

1. Common multi-hit utilities (`DivideGroup1/2`, `MultiReconstructionX/Y`) live in a shared library source (`src/MultiHitAlgo.cc`).
2. Input `TTree` / branch-name dependencies are isolated in `FROST::TreeReader` (`src/TreeReader.cc`).
3. Mapping graphs are written by `mapfunction_tool` and directly read by `FROST_reconstruction`.

---

## Build

You need a working ROOT environment (`root-config` must be available).

```bash
make
```

This builds:

- `./mapfunction_tool`
- `./FROST_reconstruction`

### ROOT dictionary note (required for groupx/groupy branches)

`FROST_reconstruction` writes `groupx`/`groupy` as `std::vector<std::vector<int>>` branches. This requires a ROOT dictionary
(CollectionProxy). The build generates a dictionary using:

- `dict/LinkDef.h`
- `include/DictTypes.h`
- auto-generated: `src/dict.cxx`, `src/dict_rdict.pcm`

If you see errors about `vector<vector<int>>` CollectionProxy, rebuild from scratch:

```bash
make clean
make
```

---

## mapfunction_tool

### Single-hit map (gx/gy + per-fiber response)

This mode produces:

- `gx` : `TGraphErrors` mapping **xg -> x_true**
- `gy` : `TGraphErrors` mapping **yg -> y_true**
- `gmcx_i` (i=0..131): `TGraph` mapping **x_true -> expected p.e.** for each X fiber
- `gmcy_i` (i=0..139): `TGraph` mapping **y_true -> expected p.e.** for each Y fiber

The graphs are completed on a 1-mm grid using symmetry.

Example:

```bash
./mapfunction_tool \
  --mode single \
  --alpha 1.0 \
  --in-x-pattern "/path/to/singlehitmuon...x{label}mmy0mm....root" \
  --in-y-pattern "/path/to/singlehitmuon...y{label}mmx0mm....root" \
  --out single_map.root
```

Optional:

- `--max-events N` (use only the first N entries per file)

### Two-hit map (merged-peak calibration)

This mode produces:

- `g1xm`, `g1xp` : `TGraphErrors` mapping **xg_candidate -> x_true** for the merged-peak case
- `g1ym`, `g1yp` : `TGraphErrors` mapping **yg_candidate -> y_true** for the merged-peak case

The symmetry used matches the original reconstruction logic:

- `xm(-x) = -xp(+x)`
- `xp(-x) = -xm(+x)`

Example:

```bash
./mapfunction_tool \
  --mode two \
  --alpha 1.0 \
  --in-x-pattern "/path/to/twohit...x{label}mmy-5mmto5mm....root" \
  --in-y-pattern "/path/to/twohit...y{label}mmx-5mmto5mm....root" \
  --out two_map.root
```

Optional:

- `--threshold1 10` (default 10)
- `--threshold2 20` (default 20)
- `--max-events N`

Parameters written to output:

- `TParameter<double> alpha`
- `TParameter<double> threshold1` (two-hit map only)
- `TParameter<double> threshold2` (two-hit map only)

---

## FROST_reconstruction

This tool reads an input ROOT file and reconstructs `(x,y)` event-by-event.

### Input modes

You must specify exactly one of:

- `--mc`
- `--data`

The default input tree name depends on the mode:

- MC: `wls`
- Data: `tree`

You can override this with:

- `--tree <tree_name>`

Photon counts are read via `FROST::TreeReader`.

#### MC input format

Expected branches:

- `lightyieldx[132]`
- `lightyieldy[140]`

#### Data input format

Expected branch:

- `lightyield[272][8]`

Internally, this is converted by `TreeReader` into MC-like X/Y views for each bunch:

- X side: `lightyield[140]` to `lightyield[271]` correspond to MC-like `lightyieldx[131]` to `lightyieldx[0]`
- Y side: `lightyield[0]` to `lightyield[139]` correspond to MC-like `lightyieldy[139]` to `lightyieldy[0]`

That is, both X and Y channel orders are reversed when converting from Data format.

The second index `[8]` is the bunch index:

- bunch = `0,1,2,3,4,5,6,7`

Also, negative light-yield values are clamped to zero inside `TreeReader` before reconstruction is performed.


### Reconstruction algorithm

For each event and for each bunch used by the selected input mode:

1. Compute single-hit `xg`, `yg` and the corresponding `(x,y)` using the single-hit map.
2. Compute Poisson reduced chi-square.
3. If reduced reduced chi-square exceeds `--chi2-threshold`, generate multi-hit candidates using:
   - `DivideGroup1/DivideGroup2` (grouping)
   - `MultiReconstructionX/MultiReconstructionY` (candidate xg/yg)
   - two-hit maps (`g1xm/g1xp/g1ym/g1yp`) to convert candidates to position when `rectype` indicates +/- mapping

The map parameters are taken from the map ROOT files:

- `alpha` is read from both single-map and two-map and **must match** (mismatch => error).
- `threshold1`/`threshold2` are read from the **two-map** file.

For Data, reconstruction is performed independently for all 8 bunches.
For MC, reconstruction is performed only once, and the results are stored in bunch slot 0.

Example:

```bash
./FROST_reconstruction \
  --mc \
  --in input.root \
  --out out.root \
  --single-map single_map.root \
  --two-map two_map.root
```

```bash
./FROST_reconstruction \
  --data \
  --in input.root \
  --out out.root \
  --single-map single_map.root \
  --two-map two_map.root
```

Important options:

- `--chi2-threshold 1.26` (default 1.26)
- `--max-events N`
- `--tree <tree_name>`

### Output

The output ROOT file:

- Copies all trees and objects from the input file, except for the source input tree.
- Writes a single output tree named **`frost`**, which is a clone of the source input tree (all original active branches preserved) plus additional reconstruction branches.
- If the input mode is Data and the source tree has a `waveform` branch, that branch is not copied into `frost` to reduce output size.

Added branches:

- `is_hit` (`std::vector<int>`)
  - size = 8
  - `is_hit[b] = 1` if `max(lightyield_x) >= 10` and `max(lightyield_y) >= 10` in bunch `b`
  - otherwise `0`
  - for MC, only `is_hit[0]` is meaningful

- `xg` (`std::vector<std::vector<double>>`)
  - `xg[bunch][candidate]`

- `yg` (`std::vector<std::vector<double>>`)
  - `yg[bunch][candidate]`

- `x_rec` (`std::vector<std::vector<double>>`)
  - `x_rec[bunch][candidate]`
  - aligned with `xg[bunch]`

- `y_rec` (`std::vector<std::vector<double>>`)
  - `y_rec[bunch][candidate]`
  - aligned with `yg[bunch]`

- `chi2` (`std::vector<double>`)
  - reduced chi-square used for multihit decision
  - one value per bunch

- `is_multihit` (`std::vector<int>`)
  - one value per bunch
  - `1` if multi-hit candidate generation was triggered, else `0`

- `x_rectype` (`std::vector<std::vector<int>>`)
  - `x_rectype[bunch][candidate]`

- `y_rectype` (`std::vector<std::vector<int>>`)
  - `y_rectype[bunch][candidate]`

- `groupx` (`std::vector<std::vector<std::vector<int>>>`)
  - `groupx[bunch][group][bin_index]`

- `groupy` (`std::vector<std::vector<std::vector<int>>>`)
  - `groupy[bunch][group][bin_index]`

For single-hit events, the candidate vectors for that bunch contain a single entry and `rectype = 2`.

### Definition of `rectype`

`rectype` indicates which mapping function was used to convert reconstructed `xg` / `yg` to position:

- `rectype = 0`
  - use the “minus-side” two-hit mapping
  - X: `g1xm`
  - Y: `g1ym`

- `rectype = 1`
  - use the “plus-side” two-hit mapping
  - X: `g1xp`
  - Y: `g1yp`

- `rectype = 2`
  - use the single-hit mapping
  - X: `gx`
  - Y: `gy`

In other words:

- `0/1` mean the candidate came from the multi-hit logic and is mapped with the two-hit map
- `2` means the candidate is treated with the single-hit map

Parameters written to output:

- `TParameter<double> alpha`
- `TParameter<double> threshold1`
- `TParameter<double> threshold2`
- `TParameter<double> chi2_threshold`

---

## Notes

- `FROST_reconstruction` reports progress (percentage of processed events) to stderr.
- Data input is handled bunch-by-bunch.
- Negative light-yield values from Data are floored to zero in `TreeReader`.
- The ROOT dictionary is required because `groupx` / `groupy` branches use nested STL containers.
