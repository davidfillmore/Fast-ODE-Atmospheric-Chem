# The eleven TS1 paper scenarios

Each scenario contains cells **0-63**, the first 64 cells used in the paper,
in their original order. The original datasets had 1,000 cells. CSV values
are copied verbatim; no concentrations or reference solutions were recomputed
during the import. Scenario comments and `references.log` retain the original
1,000-cell descriptions and validation results.

| Scenario | Host interval (s) | Intervals | Snapshot every (intervals) | Original group |
|---|---:|---:|---:|---|
| clean-10min | 600 | 144 | 36 | disclosed |
| clean-48h-5min | 300 | 576 | 72 | held out |
| diurnal-10min | 600 | 144 | 36 | disclosed |
| diurnal-30min | 1800 | 48 | 12 | disclosed |
| diurnal-5min | 300 | 288 | 72 | disclosed |
| polluted-10min | 600 | 144 | 36 | disclosed |
| polluted-48h-30min | 1800 | 96 | 12 | held out |
| summer-15min | 900 | 96 | 24 | held out |
| tier1-30s | 30 | 101 | final only | disclosed |
| tropics-20min | 1200 | 72 | 18 | held out |
| winter60-10min | 600 | 144 | 36 | held out |


`scenario.txt` contains `key = value` settings, with paths relative to its
directory. `cells.csv` sets each cell's temperature, pressure, concentrations,
and rate parameters. `photolysis.csv` contains `time_s,factor`; each host
interval holds every photolysis parameter at its cell value times the factor
at the interval start. Other rates and the cell conditions stay fixed.

`reference_final.csv` has columns `temperature,pressure,<species...>` and one
row per cell. Where present, `reference_snapshots.csv` has columns
`step,cell,temperature,pressure,<species...>` and contains all 64 cells at
each recorded time, including the final time. The references were generated
with MICM's six-stage Rosenbrock method at `rtol = 1e-10`,
`atol = 1e-23 mol m-3`, and checked at `1e-9` and `1e-22`.

The scoring floor is `1e-14 mol m-3`. The drivers score all saved snapshots
and match reference columns by species name. See `fastode.ts1.score` for
the significant-digits calculation and the finite, nonnegative-state gate.

The former held-out scenarios are included for public reproduction. Their
original grouping is retained only to explain the paper's development protocol.

See `../../provenance/import-manifest.json` for original and imported hashes.
