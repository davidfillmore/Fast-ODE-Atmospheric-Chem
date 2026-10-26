# Baseline validation

`python validation/tier1_final.py` compares the Python restarting ROS3 loop
against the included MICM Tier 1 final states: 100 cells, 101 calls of 30 s,
`rtol = 1e-6`, and `atol = 1e-3 mol m-3`. It prints differences and counts;
it is a diagnostic script, not an assertion-based test. The original report
is in `../reference-results/validation.md`.

To repeat the additional comparisons against MICM itself, obtain MICM at
`e15898ce482a979594397144c76258e1ada6b79d` and build the two supplied C++
drivers. For example, with a C++20 compiler and that checkout at `$MICM_SRC`:

```sh
c++ -std=c++20 -O3 -DNDEBUG -pthread \
  -I "$MICM_SRC/include" -I "$MICM_SRC/test/benchmark" \
  validation/micm_scenario.cpp -o /tmp/micm_scenario
c++ -std=c++20 -O3 -DNDEBUG -pthread \
  -I "$MICM_SRC/include" -I "$MICM_SRC/test/benchmark" \
  validation/micm_be_scenario.cpp -o /tmp/micm_be_scenario
python validation/scenario_agreement.py data/scenarios/diurnal-10min \
  --micm /tmp/micm_scenario --micm-be /tmp/micm_be_scenario \
  --cells 32 --steps 36 --rtol 1e-3 --atol 1e-15
```

On macOS a compiler may also need `SDKROOT` set to the installed macOS SDK.
This optional comparison requires an external MICM checkout. The regular
Python tests, examples, and experiment drivers do not.
