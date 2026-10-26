# Validation of the baseline prototypes against MICM (2026-10-20)

MICM at commit e15898ce. `validation/tier1_final.py`: MICM's Tier 1 run (100 cells
of ts1_cells_100_seed1.csv in one group, 101 calls of 30 s, three-stage Rosenbrock,
relative 1e-6, absolute 1e-3). `validation/scenario_agreement.py`: MICM's scenario
driver (rosenbrock3, rosenbrock4dae) and a backward-Euler build of it
(`validation/micm_be_scenario.cpp`), one group of 32 cells, against the prototypes.

## Tier 1 final state (restarting loop, ROS3)

- Steps: 909 attempted, 909 accepted, 0 rejected (9 per call), as MICM.
- State: 206 of 209 species agree to 5e-12 or better. CL2O2, S, and SO (about
  1e-22, 1e-26, 1e-21 mol m-3, nineteen orders below the absolute tolerance, so
  uncontrolled) differ by 3.5e-5, 4.6e-6, 4.9e-7. Changing the prototype's own LU
  ordering (COLAMD to MMD) moves them by the same amounts (4.8e-5, 4.6e-6, 4.9e-7);
  with MMD, S and SO agree with MICM to 2e-15 and 4e-15. Roundoff in species the
  error control does not see.

## diurnal-10min, first 32 cells, 36 calls of 600 s (06:00 to 12:00), relative 1e-3, absolute 1e-15

| Solver | Largest relative state difference (values >= 1e-14) | MICM counts | Prototype counts |
|---|---|---|---|
| ROS3 loop | 8.2e-12 | 491 steps, 982 forcings, 491 Jacobians, 491 LU, 1473 solves | identical |
| RODAS3 loop | 2.7e-12 | 428 steps, 1284 forcings, 428 Jacobians, 428 LU, 1712 solves | identical |
| backward Euler | 3.1e-14 | 119 Newton iterations, 36 steps | identical |

## diurnal-30min, first 32 cells, 24 calls of 1800 s, relative 1e-4, absolute 1e-18

| Solver | Largest relative state difference | MICM | Prototype |
|---|---|---|---|
| ROS3 loop | 4.6e-12 | 1116 steps (1025 accepted), 2141 forcings, 1116 Jacobians, 1116 LU, 3348 solves | 1116 steps (1025 accepted, 91 rejected), 2141 forcings, 1025 Jacobians, 1116 LU, 3348 solves |
| RODAS3 loop | 3.0e-12 | 741 steps (731 accepted), 2213 forcings, 741 Jacobians, 741 LU, 2964 solves | 741 (731, 10 rejected), 2213, 731 Jacobians, 741, 2964 |
| backward Euler | 2.0e-14 | 103 Newton iterations, 24 steps | identical |

The two count differences are bookkeeping: MICM's in-place LU re-evaluates the
Jacobian after a rejection (the prototypes keep it, as the paper's algorithm
does), and MICM's rejection statistic counts only rejections after a call's first
accepted step.
