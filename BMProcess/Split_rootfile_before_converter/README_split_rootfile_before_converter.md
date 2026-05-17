# ROOT Time Split Utilities

This directory contains two small C++ utilities for splitting ROOT files by timestamp.

## Programs

- `split_daily.cpp`: split BM or PM ROOT files in an input directory into daily ROOT files.
- `split_at_time.cpp`: split one BM/PM ROOT file, or one WG run directory, into two parts at a specified time.

Both programs use local time when converting between date strings and Unix time. Set `TZ=Asia/Tokyo` before running if you need JST explicitly.

```bash
export TZ=Asia/Tokyo
```

## Build

Example executable names are assumed to be:

```text
split_daily
split_at_time
```

Build with your Makefile, for example:

```bash
make
```

The programs require ROOT and `root-config`.

## split_daily

Split all ROOT files in an input directory into one output ROOT file per day.

### Usage

```bash
./split_daily -bm <input_dir> <output_dir>
./split_daily -pm <input_dir> <output_dir>
```

### Options

| Option | Input | Timestamp branch | Output name |
|---|---|---|---|
| `-bm` | directory containing BM ROOT files | `unixtime` | `BMBSD_YYYY-MM-DD_00-00-00.root` |
| `-pm` | directory containing PM ROOT files | `BMBSD.unixtime` | `PMBSD_YYYY-MM-DD_00-00-00.root` |

Entries with timestamp `<= 0` are ignored.

### Examples

```bash
./split_daily -bm ./BM_input ./BM_daily
./split_daily -pm ./PM_input ./PM_daily
```

Example outputs:

```text
BMBSD_2026-01-24_00-00-00.root
PMBSD_2026-01-24_00-00-00.root
```

## split_at_time

Split data into two outputs at a specified timestamp.

The split time format must be:

```text
YYYY-MM-DD_HH-MM-SS
```

The following format is also accepted:

```text
YYYY-MM-DD_HH-MM:SS
```

The first output contains entries before the split time.  
The second output contains entries at or after the split time.

## BM / PM mode

### Usage

```bash
./split_at_time -bm <input.root> <split_time> <output_dir>
./split_at_time -pm <input.root> <split_time> <output_dir>
```

### Options

| Option | Input | Timestamp branch | Output name |
|---|---|---|---|
| `-bm` | one BM ROOT file | `unixtime` | `BMBSD_YYYY-MM-DD_HH-MM-SS.root` |
| `-pm` | one PM ROOT file | `BMBSD.unixtime` | `PMBSD_YYYY-MM-DD_HH-MM-SS.root` |

The input file name must match the selected mode, for example:

```text
BMBSD_2026-01-24_00-00-00.root
PMBSD_2026-01-24_00-00-00.root
```

### Examples

```bash
./split_at_time -bm BMBSD_2026-01-24_00-00-00.root 2026-01-24_12-00-00 ./output
./split_at_time -pm PMBSD_2026-01-24_00-00-00.root 2026-01-24_12-00-00 ./output
```

Example BM outputs:

```text
output/BMBSD_2026-01-24_00-00-00.root
output/BMBSD_2026-01-24_12-00-00.root
```

## WG mode

WG data are handled as a run directory, not as a single ROOT file.

### Usage

```bash
./split_at_time -wg <input_dir> <split_time> <output_dir>
```

### Input directory format

The input directory must be named like:

```text
physics_run_YYYY-MM-DD_HH-MM-SS_0
```

Example:

```text
physics_run_2026-01-24_00-00-00_0
```

Inside the directory, the program looks for existing DIF files from `dif_0` to `dif_7`:

```text
physics_run_2026-01-24_00-00-00_0_ecal_dif_0_tree.root
physics_run_2026-01-24_00-00-00_0_ecal_dif_1_tree.root
...
physics_run_2026-01-24_00-00-00_0_ecal_dif_7_tree.root
```

Missing DIF files are skipped.

### Timestamp branch

WG splitting uses:

```text
bsd.timestamp
```

The `bsd` tree is used to build the selected entry list. The same selected entries are copied for all TTrees in each DIF ROOT file.

Entries with `timestamp <= 0` are ignored.

### Example

```bash
./split_at_time -wg physics_run_2026-01-24_00-00-00_0 2026-01-24_12-30-03 ./output
```

Example outputs:

```text
output/physics_run_2026-01-24_00-00-00_0/
  physics_run_2026-01-24_00-00-00_0_ecal_dif_0_tree.root
  physics_run_2026-01-24_00-00-00_0_ecal_dif_1_tree.root
  ...

output/physics_run_2026-01-24_12-30-03_0/
  physics_run_2026-01-24_12-30-03_0_ecal_dif_0_tree.root
  physics_run_2026-01-24_12-30-03_0_ecal_dif_1_tree.root
  ...
```

The first output directory contains entries with:

```text
timestamp < split_time
```

The second output directory contains entries with:

```text
timestamp >= split_time
```

## Notes

- Use an output directory different from the input location to avoid overwriting original data.
- The tools print the number of processed files, entries, invalid timestamps, and output entries.
- For BM/PM, only the TTree named `tree` is copied.
- For WG, all TTrees in each DIF ROOT file are copied using the `bsd.timestamp` selection.
