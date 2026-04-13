# match_pm_to_bm.cpp

`match_pm_to_bm` matches Proton Module events to Baby MIND files using `spill/spillnum` and `unixtime`, then writes Proton Module output files with Baby MIND-aligned file names.

## Build

```bash
make
```

## Run

```bash
./match_pm_to_bm \
  /path/to/babymind/2-BMBSD \
  /path/to/protonmodule/4-PMBSD \
  /path/to/output_dir
```

## Input

- Baby MIND files: `BMBSD_*.root`
- Proton Module files: `PMBSD_*.root`
- Both must contain a TTree named `tree`

Used branches:

- Baby MIND: `spillnum`, `unixtime`
- Proton Module: `spill`, `unixtime`

## Matching rule

A Proton Module event is matched to a Baby MIND event when:

- `spillnum % 32768 == spill % 32768`
- `abs(BM.unixtime - PM.unixtime) <= 3000`
- `BM.unixtime != -1`

Each Proton Module event is used at most once.

## Output

For each Baby MIND file:

- input: `BMBSD_2025-11-29_13-46-59_Run0.root`
- output: `PMBSD_2025-11-29_13-46-59_Run0.root`

The output tree keeps the original Proton Module branch structure.
