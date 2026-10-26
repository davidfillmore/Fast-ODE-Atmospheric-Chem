# TS1 summary

528 runs on 64 cells; counts per interval and per cell over the eleven scenarios.

| config | rtol | atol | L2 | L3 | mean D | lowest D (scenario) | lowest species | steps | rej | LU | solves | rhs | newton | failed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| alg | 0.03 | 3e-16 | 10 | 1 | 2.80 | 1.98 (diurnal-5min) | BRO 1.10 (diurnal-5min) | 1.72 | 0.06 | 1.72 | 5.15 | 3.37 | 0.00 | 0 |
| alg | 0.01 | 1e-16 | 11 | 8 | 3.41 | 2.70 (diurnal-30min) | ISOP 1.80 (polluted-48h-30min) | 2.67 | 0.25 | 2.67 | 8.00 | 5.08 | 0.00 | 0 |
| alg | 0.005 | 5e-17 | 11 | 11 | 3.65 | 3.21 (tropics-20min) | ISOP 2.12 (polluted-48h-30min) | 3.96 | 0.53 | 3.96 | 11.88 | 7.39 | 0.00 | 0 |
| alg | 0.002 | 2e-17 | 11 | 11 | 4.12 | 3.68 (tropics-20min) | ISOP 2.47 (diurnal-30min) | 6.90 | 1.13 | 6.90 | 20.70 | 12.67 | 0.00 | 0 |
| alg | 0.001 | 1e-17 | 11 | 11 | 4.38 | 3.90 (diurnal-30min) | ISOP 2.71 (diurnal-30min) | 10.27 | 1.81 | 10.27 | 30.81 | 18.73 | 0.00 | 0 |
| alg | 0.0005 | 5e-18 | 11 | 11 | 4.66 | 4.12 (diurnal-30min) | ISOP 2.91 (diurnal-30min) | 14.27 | 2.54 | 14.27 | 42.82 | 26.00 | 0.00 | 0 |
| alg-a15 | 0.03 | 1e-15 | 8 | 1 | 2.54 | 1.81 (clean-10min) | BRO 1.03 (diurnal-10min) | 1.63 | 0.05 | 1.63 | 4.90 | 3.22 | 0.00 | 0 |
| alg-a15 | 0.01 | 1e-15 | 11 | 6 | 3.24 | 2.58 (diurnal-30min) | ISOP 1.76 (polluted-48h-30min) | 2.28 | 0.17 | 2.28 | 6.85 | 4.40 | 0.00 | 0 |
| alg-a15 | 0.005 | 1e-15 | 11 | 10 | 3.51 | 2.95 (diurnal-30min) | MTERP 2.09 (polluted-48h-30min) | 2.98 | 0.31 | 2.98 | 8.95 | 5.66 | 0.00 | 0 |
| alg-a15 | 0.002 | 1e-15 | 11 | 11 | 3.81 | 3.49 (tropics-20min) | ISOP 2.20 (polluted-48h-30min) | 4.34 | 0.58 | 4.34 | 13.03 | 8.10 | 0.00 | 0 |
| alg-a15 | 0.001 | 1e-15 | 11 | 11 | 4.13 | 3.72 (diurnal-30min) | ISOP 2.44 (polluted-48h-30min) | 5.77 | 0.85 | 5.77 | 17.30 | 10.68 | 0.00 | 0 |
| alg-a15 | 0.0005 | 1e-15 | 11 | 11 | 4.37 | 3.85 (diurnal-30min) | ISOP 2.65 (diurnal-30min) | 7.38 | 1.12 | 7.38 | 22.15 | 13.65 | 0.00 | 0 |
| alg-a20 | 0.03 | 1e-20 | 11 | 4 | 3.11 | 2.44 (diurnal-30min) | PHENO 1.67 (clean-48h-5min) | 2.46 | 0.19 | 2.46 | 7.39 | 4.73 | 0.00 | 0 |
| alg-a20 | 0.01 | 1e-20 | 11 | 11 | 3.55 | 3.04 (tropics-20min) | ISOP 2.14 (diurnal-30min) | 4.28 | 0.56 | 4.28 | 12.85 | 8.00 | 0.00 | 0 |
| alg-a20 | 0.005 | 1e-20 | 11 | 11 | 3.82 | 3.43 (summer-15min) | ISOP 2.40 (diurnal-30min) | 6.30 | 0.96 | 6.30 | 18.89 | 11.63 | 0.00 | 0 |
| alg-a20 | 0.002 | 1e-20 | 11 | 11 | 4.30 | 3.87 (diurnal-30min) | ISOP 2.63 (diurnal-30min) | 10.25 | 1.73 | 10.25 | 30.76 | 18.78 | 0.00 | 0 |
| alg-a20 | 0.001 | 1e-20 | 11 | 11 | 4.59 | 4.07 (diurnal-30min) | ISOP 2.84 (diurnal-30min) | 13.92 | 2.38 | 13.92 | 41.75 | 25.45 | 0.00 | 0 |
| alg-a20 | 0.0005 | 1e-20 | 11 | 11 | 4.78 | 4.32 (diurnal-30min) | ISOP 3.09 (diurnal-30min) | 17.74 | 2.96 | 17.74 | 53.23 | 32.52 | 0.00 | 0 |
| alg-elementary | 0.01 | 1e-16 | 11 | 3 | 2.96 | 2.09 (diurnal-30min) | HYDRALD 1.31 (diurnal-30min) | 2.06 | 0.28 | 2.06 | 6.19 | 3.85 | 0.00 | 0 |
| alg-elementary | 0.002 | 2e-17 | 11 | 10 | 3.63 | 3.14 (diurnal-30min) | ISOP 1.95 (diurnal-30min) | 5.16 | 1.20 | 5.16 | 15.47 | 9.12 | 0.00 | 0 |
| alg-elementary | 0.0005 | 5e-18 | 11 | 11 | 4.17 | 3.62 (diurnal-30min) | ISOP 2.45 (diurnal-30min) | 10.30 | 2.61 | 10.30 | 30.89 | 17.98 | 0.00 | 0 |
| alg-nomemory | 0.01 | 1e-16 | 11 | 3 | 3.08 | 2.29 (tropics-20min) | HYDRALD 1.55 (tropics-20min) | 2.54 | 0.25 | 2.54 | 7.62 | 4.84 | 0.00 | 0 |
| alg-nomemory | 0.002 | 2e-17 | 11 | 11 | 3.84 | 3.44 (clean-10min) | DMS 2.32 (winter60-10min) | 6.81 | 1.13 | 6.81 | 20.43 | 12.49 | 0.00 | 0 |
| alg-nomemory | 0.0005 | 5e-18 | 11 | 11 | 4.46 | 3.92 (diurnal-30min) | ISOP 2.91 (diurnal-30min) | 14.05 | 2.54 | 14.05 | 42.15 | 25.57 | 0.00 | 0 |
| alg-noquant | 0.01 | 1e-16 | 11 | 3 | 3.01 | 2.29 (tropics-20min) | ISOP 1.50 (polluted-48h-30min) | 2.53 | 0.30 | 2.53 | 7.59 | 4.77 | 0.00 | 0 |
| alg-noquant | 0.002 | 2e-17 | 11 | 11 | 3.72 | 3.38 (clean-10min) | DMS 2.32 (winter60-10min) | 6.80 | 1.22 | 6.80 | 20.39 | 12.38 | 0.00 | 0 |
| alg-noquant | 0.0005 | 5e-18 | 11 | 11 | 4.34 | 3.85 (diurnal-30min) | ISOP 2.74 (diurnal-30min) | 13.92 | 2.62 | 13.92 | 41.76 | 25.22 | 0.00 | 0 |
| alg-reset | 0.01 | 1e-16 | 11 | 9 | 4.01 | 3.15 (diurnal-30min) | ISOP 1.91 (polluted-48h-30min) | 13.17 | 0.02 | 13.17 | 39.51 | 26.32 | 0.00 | 0 |
| alg-reset | 0.002 | 2e-17 | 11 | 11 | 4.33 | 3.56 (diurnal-30min) | ISOP 2.33 (diurnal-30min) | 15.57 | 0.20 | 15.57 | 46.72 | 30.95 | 0.00 | 0 |
| alg-reset | 0.0005 | 5e-18 | 11 | 11 | 4.65 | 3.98 (diurnal-30min) | ISOP 2.73 (diurnal-30min) | 19.20 | 0.72 | 19.20 | 57.60 | 37.68 | 0.00 | 0 |
| be | 1e-06 | 0.001 | 0 | 0 | 0.36 | -0.51 (diurnal-30min) | ISOP -1.92 (diurnal-30min) | 1.00 | 0.00 | 2.00 | 2.00 | 2.00 | 2.00 | 0 |
| be-tight | 0.001 | 1e-15 | 0 | 0 | 0.38 | -0.49 (polluted-48h-30min) | ISOP -1.84 (diurnal-30min) | 1.00 | 0.00 | 2.72 | 2.72 | 2.72 | 2.72 | 0 |
| be-tight | 1e-06 | 1e-15 | 0 | 0 | 0.38 | -0.49 (polluted-48h-30min) | ISOP -1.84 (diurnal-30min) | 1.00 | 0.00 | 3.29 | 3.29 | 3.29 | 3.29 | 0 |
| restart | 0.03 | 1e-15 | 11 | 6 | 3.24 | 2.32 (diurnal-30min) | ISOP 1.01 (polluted-48h-30min) | 9.40 | 0.00 | 9.40 | 28.19 | 18.79 | 0.00 | 0 |
| restart | 0.01 | 1e-15 | 11 | 6 | 3.38 | 2.52 (diurnal-30min) | ISOP 1.25 (polluted-48h-30min) | 9.69 | 0.00 | 9.69 | 29.06 | 19.37 | 0.00 | 0 |
| restart | 0.003 | 1e-15 | 11 | 8 | 3.61 | 2.81 (diurnal-30min) | ISOP 1.48 (polluted-48h-30min) | 10.14 | 0.00 | 10.14 | 30.41 | 20.27 | 0.00 | 0 |
| restart | 0.001 | 1e-15 | 11 | 9 | 3.81 | 3.02 (diurnal-30min) | ISOP 1.78 (polluted-48h-30min) | 10.71 | 0.00 | 10.71 | 32.12 | 21.41 | 0.00 | 0 |
| restart | 0.0003 | 1e-15 | 11 | 11 | 4.07 | 3.34 (diurnal-30min) | ISOP 2.12 (diurnal-30min) | 11.53 | 0.00 | 11.53 | 34.60 | 23.06 | 0.00 | 0 |
| restart | 0.0001 | 1e-15 | 11 | 11 | 4.29 | 3.66 (diurnal-30min) | ISOP 2.45 (diurnal-30min) | 12.49 | 0.00 | 12.49 | 37.46 | 24.97 | 0.00 | 0 |
| restart-a20 | 0.01 | 1e-20 | 11 | 8 | 3.53 | 2.77 (diurnal-30min) | ISOP 1.49 (polluted-48h-30min) | 10.60 | 0.13 | 10.60 | 31.81 | 21.08 | 0.00 | 0 |
| restart-a20 | 0.001 | 1e-20 | 11 | 11 | 4.10 | 3.50 (diurnal-30min) | ISOP 2.24 (diurnal-30min) | 13.99 | 0.73 | 13.99 | 41.96 | 27.25 | 0.00 | 0 |
| restart-a20 | 0.0003 | 1e-20 | 11 | 11 | 4.48 | 3.97 (diurnal-30min) | ISOP 2.79 (diurnal-30min) | 16.93 | 1.14 | 16.93 | 50.79 | 32.72 | 0.00 | 0 |
| restart-p | 0.03 | 3e-16 | 11 | 6 | 3.23 | 2.43 (diurnal-30min) | ISOP 1.07 (polluted-48h-30min) | 9.48 | 0.00 | 9.48 | 28.45 | 18.97 | 0.00 | 0 |
| restart-p | 0.01 | 1e-16 | 11 | 7 | 3.40 | 2.54 (diurnal-30min) | ISOP 1.27 (diurnal-30min) | 9.98 | 0.02 | 9.98 | 29.95 | 19.95 | 0.00 | 0 |
| restart-p | 0.003 | 3e-17 | 11 | 8 | 3.70 | 2.92 (diurnal-30min) | ISOP 1.66 (polluted-48h-30min) | 10.93 | 0.12 | 10.93 | 32.78 | 21.74 | 0.00 | 0 |
| restart-p | 0.001 | 1e-17 | 11 | 11 | 3.97 | 3.24 (diurnal-30min) | ISOP 2.02 (diurnal-30min) | 12.45 | 0.40 | 12.45 | 37.34 | 24.49 | 0.00 | 0 |
| restart-p | 0.0003 | 3e-18 | 11 | 11 | 4.33 | 3.75 (diurnal-30min) | ISOP 2.56 (diurnal-30min) | 15.13 | 0.95 | 15.13 | 45.39 | 29.31 | 0.00 | 0 |
| restart-p | 0.0001 | 1e-18 | 11 | 11 | 4.70 | 4.22 (diurnal-30min) | ISOP 3.10 (diurnal-30min) | 18.22 | 1.27 | 18.22 | 54.66 | 35.18 | 0.00 | 0 |

## Loosest relative tolerance reaching each level on all eleven scenarios

| level | config | rtol | mean D | lowest D | lowest species | steps | rej | LU | solves | rhs |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | alg | 0.01 | 3.41 | 2.70 | ISOP 1.80 | 2.67 | 0.25 | 2.67 | 8.00 | 5.08 | work 75.1 |
| 2 | alg-a20 | 0.03 | 3.11 | 2.44 | PHENO 1.67 | 2.46 | 0.19 | 2.46 | 7.39 | 4.73 | work 69.7 |
| 2 | alg-a15 | 0.01 | 3.24 | 2.58 | ISOP 1.76 | 2.28 | 0.17 | 2.28 | 6.85 | 4.40 | work 64.7 |
| 2 | restart-p | 0.03 | 3.23 | 2.43 | ISOP 1.07 | 9.48 | 0.00 | 9.48 | 28.45 | 18.97 | work 274.2 |
| 2 | restart | 0.03 | 3.24 | 2.32 | ISOP 1.01 | 9.40 | 0.00 | 9.40 | 28.19 | 18.79 | work 271.7 |
| 2 | restart-a20 | 0.01 | 3.53 | 2.77 | ISOP 1.49 | 10.60 | 0.13 | 10.60 | 31.81 | 21.08 | work 305.5 |
| 3 | alg | 0.005 | 3.65 | 3.21 | ISOP 2.12 | 3.96 | 0.53 | 3.96 | 11.88 | 7.39 | work 110.3 |
| 3 | alg-a20 | 0.01 | 3.55 | 3.04 | ISOP 2.14 | 4.28 | 0.56 | 4.28 | 12.85 | 8.00 | work 119.3 |
| 3 | alg-a15 | 0.002 | 3.81 | 3.49 | ISOP 2.20 | 4.34 | 0.58 | 4.34 | 13.03 | 8.10 | work 120.9 |
| 3 | restart-p | 0.001 | 3.97 | 3.24 | ISOP 2.02 | 12.45 | 0.40 | 12.45 | 37.34 | 24.49 | work 356.6 |
| 3 | restart | 0.0003 | 4.07 | 3.34 | ISOP 2.12 | 11.53 | 0.00 | 11.53 | 34.60 | 23.06 | work 333.4 |
| 3 | restart-a20 | 0.001 | 4.10 | 3.50 | ISOP 2.24 | 13.99 | 0.73 | 13.99 | 41.96 | 27.25 | work 398.6 |

## Steps per interval at equal mean digits (interpolated)

| config | D = 3.0 | D = 3.5 | D = 4.0 | D = 4.5 |
|---|---|---|---|---|
| alg | 1.98 | 3.09 | 6.00 | 11.83 |
| alg-a15 | 2.04 | 2.97 | 5.14 | -- |
| alg-a20 | -- | 4.00 | 7.58 | 12.70 |
| restart | -- | 9.92 | 11.31 | -- |
| restart-a20 | -- | -- | 13.31 | -- |
| restart-p | -- | 10.29 | 12.65 | 16.49 |
| alg-elementary | 2.19 | 4.30 | 8.24 | -- |
| alg-nomemory | -- | 4.40 | 8.24 | -- |
| alg-noquant | -- | 4.99 | 9.37 | -- |
| alg-reset | -- | -- | -- | 17.44 |
| be | -- | -- | -- | -- |
| be-tight | -- | -- | -- | -- |

## Rejections per interval at equal mean digits

| config | D = 3.0 | D = 3.5 | D = 4.0 | D = 4.5 |
|---|---|---|---|---|
| alg | 0.09 | 0.33 | 0.94 | 2.09 |
| alg-a15 | 0.11 | 0.31 | 0.73 | -- |
| alg-a20 | -- | 0.49 | 1.20 | 2.16 |
| alg-elementary | 0.31 | 0.90 | 2.03 | -- |
| alg-nomemory | -- | 0.58 | 1.40 | -- |
| alg-noquant | -- | 0.78 | 1.71 | -- |
| alg-reset | -- | -- | -- | 0.40 |
| be | -- | -- | -- | -- |
| be-tight | -- | -- | -- | -- |
| restart | -- | 0.00 | 0.00 | -- |
| restart-a20 | -- | -- | 0.53 | -- |
| restart-p | -- | 0.03 | 0.43 | 1.09 |

## Absolute-tolerance sweep (four scenarios)

| rtol | atol | mean D | lowest species | steps | rej |
|---|---|---|---|---|---|
| 0.01 | 1e-12 | 2.18 | BRO 1.03 | 1.59 | 0.04 |
| 0.01 | 1e-13 | 2.73 | MVK 1.95 | 1.90 | 0.08 |
| 0.01 | 1e-14 | 2.93 | MACR 2.05 | 2.33 | 0.15 |
| 0.01 | 1e-15 | 2.97 | MACR 2.07 | 2.74 | 0.23 |
| 0.01 | 1e-16 | 3.17 | ISOP 2.10 | 3.32 | 0.36 |
| 0.01 | 1e-17 | 3.21 | ISOP 2.15 | 4.35 | 0.59 |
| 0.01 | 1e-18 | 3.23 | HYDRALD 2.22 | 4.80 | 0.67 |
| 0.01 | 1e-20 | 3.27 | HYDRALD 2.30 | 5.30 | 0.72 |
| 0.002 | 1e-12 | 2.90 | MVK 2.08 | 2.34 | 0.21 |
| 0.002 | 1e-13 | 3.16 | ISOP 2.08 | 3.10 | 0.34 |
| 0.002 | 1e-14 | 3.43 | ISOP 2.27 | 4.35 | 0.58 |
| 0.002 | 1e-15 | 3.60 | ISOP 2.52 | 5.54 | 0.78 |
| 0.002 | 1e-16 | 3.80 | ISOP 2.56 | 6.75 | 1.01 |
| 0.002 | 2e-17 | 3.95 | ISOP 2.55 | 9.29 | 1.60 |
| 0.002 | 1e-17 | 3.97 | ISOP 2.55 | 10.09 | 1.78 |
| 0.002 | 1e-18 | 4.03 | ISOP 2.66 | 11.33 | 1.97 |
| 0.002 | 1e-20 | 4.14 | ISOP 2.83 | 12.60 | 2.07 |

## Restarting loop on the same four scenarios

| config | rtol | atol | mean D | lowest | steps |
|---|---|---|---|---|---|
| restart | 0.03 | 1e-15 | 3.06 | ISOP 1.46 | 9.58 |
| restart | 0.01 | 1e-15 | 3.15 | ISOP 1.55 | 9.92 |
| restart | 0.003 | 1e-15 | 3.34 | ISOP 1.79 | 10.47 |
| restart | 0.001 | 1e-15 | 3.53 | ISOP 2.02 | 11.16 |
| restart | 0.0003 | 1e-15 | 3.76 | ISOP 2.28 | 12.14 |
| restart | 0.0001 | 1e-15 | 3.99 | ISOP 2.53 | 13.26 |
| restart-a13 | 0.01 | 1e-13 | 3.11 | ISOP 1.46 | 9.45 |
| restart-a13 | 0.001 | 1e-13 | 3.33 | ISOP 1.74 | 10.04 |
| restart-a13 | 0.0001 | 1e-13 | 3.54 | ISOP 1.94 | 10.88 |
| restart-p | 0.03 | 3e-16 | 3.06 | ISOP 1.45 | 9.67 |
| restart-p | 0.01 | 1e-16 | 3.21 | ISOP 1.65 | 10.26 |
| restart-p | 0.003 | 3e-17 | 3.43 | ISOP 1.87 | 11.43 |
| restart-p | 0.001 | 1e-17 | 3.68 | ISOP 2.16 | 13.35 |
| restart-p | 0.0003 | 3e-18 | 4.06 | ISOP 2.65 | 16.18 |
| restart-p | 0.0001 | 1e-18 | 4.49 | ISOP 3.13 | 19.58 |

## Backward Euler with sub-steps, against the algorithm

| scenario | solver | substeps | D | worst | newton | work |
|---|---|---|---|---|---|---|
| tier1-30s | BE | 1 | 1.91 | ISOP 1.21 | 2.0 | 41.4 |
| tier1-30s | BE | 4 | 2.52 | ISOP 1.82 | 8.0 | 165.4 |
| tier1-30s | BE | 16 | 3.13 | ISOP 2.43 | 32.0 | 661.8 |
| tier1-30s | BE | 64 | 3.73 | ISOP 3.03 | 128.0 | 2647.2 |
| tier1-30s | alg 0.03 | -- | 5.77 | ISOP 5.00 | -- | 36.4 |
| diurnal-10min | BE | 1 | 0.42 | ISOP -0.73 | 2.0 | 41.4 |
| diurnal-10min | BE | 4 | 1.19 | ISOP 0.17 | 8.0 | 165.4 |
| diurnal-10min | BE | 16 | 1.83 | ISOP 0.85 | 32.0 | 661.8 |
| diurnal-10min | BE | 64 | 2.43 | ISOP 1.47 | 128.0 | 2647.2 |
| diurnal-10min | alg 0.03 | -- | 2.63 | PHENO 2.02 | -- | 47.7 |
| diurnal-30min | BE | 1 | -0.51 | ISOP -1.92 | 2.0 | 41.4 |
| diurnal-30min | BE | 4 | 0.61 | ISOP -0.52 | 8.0 | 165.4 |
| diurnal-30min | BE | 16 | 1.33 | ISOP 0.33 | 32.0 | 661.8 |
| diurnal-30min | BE | 64 | 1.96 | ISOP 0.98 | 128.0 | 2647.2 |
| diurnal-30min | BE | 256 | 2.56 | ISOP 1.60 | 512.0 | 10588.7 |
| diurnal-30min | alg 0.03 | -- | 2.38 | ISOP 1.66 | -- | 110.0 |
