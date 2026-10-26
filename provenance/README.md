# Source provenance

Imported on 2026-10-26 from `davidfillmore/ODE`, commit
`bbb933a4447dd730c24e78e6bbaed62c3f725c96`.
The paper is `paper_fast_ode_ts1_v2`; its Appendix E identifies this Python
implementation. The longer `paper_fast_ode` also describes the included DAE
extension and small-system experiments.

All files under `fastode/`, including `data/ts1.json`, are byte-for-byte copies
of `fastode/fastode/` in that source commit. The existing unit tests are also
unchanged. No numerical algorithm, coefficients, controller constants, or
scoring formulas were changed for this import.

The experiment and validation drivers were adapted only to use repository
paths, create missing output directories, check the bundled cell limit in
the main sweep, and allow table generation from the archived result records.
The `disclosed` selector still names the original six scenarios even though
all eleven now share one data directory. The example, package metadata, and
documentation were added here.

`import-manifest.json` records each original path, its source SHA-256, the
SHA-256 after copying or selecting cells, and its final SHA-256 after the
portability edits. The files marked as selected contain cells 0 through 63;
CSV headers and selected numeric rows retain their original bytes. Scenario
files and reference-generation logs retain the original 1,000-cell metadata;
the supplied subset is documented in `data/scenarios/README.md`.

The original eleven-scenario inputs came from `benchmarks/tier3/scenarios/`
(six disclosed scenarios) and the corresponding `tier3-heldout/heldout/`
datasets (five scenarios). All eleven are included with the author's
authorization. The 100-cell baseline validation inputs are separate files
from `benchmarks/reference/`.

`reference-results/` contains the current paper records, copied unchanged.
Superseded runs were excluded. These are historical results, not claims that
every experiment was rerun during the import. Newly generated results belong
under the ignored `results/` directory.
