# MICM attribution

The TS1 mechanism in `fastode/data/ts1.json` was translated from NCAR's
Model Independent Chemistry Module (MICM), specifically
`test/benchmark/ts1_mechanism.hpp` at commit
`e15898ce482a979594397144c76258e1ada6b79d`. The rate-constant implementation
in `fastode/ts1.py` follows that version's formulas. The baseline prototypes
were checked against the same MICM version.

MICM -- Model Independent Chemistry Module
Copyright (C) 2023-2026 University Corporation for Atmospheric Research

The mechanism has been translated to JSON and array-based Python operations;
it is not an unmodified MICM source file. Scenario data and reference results
were produced for the paper's benchmark experiments.

The original [Apache 2.0 license](MICM-LICENSE) and [MICM notice](MICM-NOTICE)
are copied from that checkout. No MICM build or bundled dependency is needed
to use the Python solver.
