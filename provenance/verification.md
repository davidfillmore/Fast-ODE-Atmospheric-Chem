# Import verification, 2026-10-26

Environment: Python 3.14.6, NumPy 2.5.4, SciPy 1.18.1, pytest 9.1.1 on macOS
arm64. Exact verification dependencies are in `requirements-tested.txt`.

- The complete numerical package, existing tests, and archived result files
  match their source SHA-256 values. No solver code was changed.
- Each of the eleven datasets contains exactly source cells 0-63, with the
  matching final references and every snapshot. Selected CSV rows preserve
  the original bytes.
- `python -m pytest -q`: **20 passed**.
- `python tools/check_tableaus.py`: ROS3 and RODAS3 orders, embedded orders,
  stability limits, and index-1 values passed.
- `python examples/ts1.py`: algorithm and restarting loop both converged on
  all 101 intervals of the 64-cell example and passed the state gate.
- `python validation/tier1_final.py`: **909 attempted and accepted steps,
  zero rejections**, matching MICM. The largest relative state difference
  was `4.831e-5` in CL2O2 (an uncontrolled trace species in the original
  validation report); nine of 20,900 entries differed by more than `1e-10`.
  This diagnostic prints differences rather than asserting a tolerance.
- The command below reran `Algorithm(ROS3)` on all eleven 64-cell scenarios:

  ```sh
  python experiments/ts1_experiments.py --configs alg --rtols 1e-3 \
    --workers 3 --out results/smoke/ts1.csv
  ```

  All calls converged and every scenario passed the finite/nonnegative-state
  gate. Attempted, accepted, and rejected steps; right-hand-side and Jacobian
  evaluations; factorizations; and solves matched the corresponding archived
  `alg`, `rtol = 1e-3` records **exactly**. Overall and worst-species digits
  matched at the recorded four decimal places, as did the worst species.
- Both table generators ran against `reference-results/`. All seven TS1
  tables matched `paper_fast_ode_ts1_v2/tables/` byte for byte. The additional
  small-system, DAE, projection, and TS1-family tables matched the longer
  paper's tables as well.
- A wheel was built and installed. An isolated import loaded `fastode` from
  `site-packages` and instantiated TS1 with **209 species and 547 reactions**.
  The wheel includes the TS1 JSON and the GPL and MICM license/notice files.

The complete 588-run tolerance/ablation suite and the optional live MICM
comparisons were not rerun during this import. The original result records
are retained separately from fresh outputs.
