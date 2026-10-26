"""Summaries and LaTeX tables from the TS1 results, for the paper.

    python3 experiments/make_tables.py [--paper ../paper_fast_ode/tables]

Counts are per interval and per 32-cell group (so per cell), summed over the
scenarios and divided by their total number of intervals. The absolute-tolerance
sweep (results/ts1/atol.csv) runs on four scenarios; its table also takes the
main runs at the same tolerances on those four.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCEN_ORDER = ["tier1-30s", "diurnal-30min", "diurnal-10min", "diurnal-5min", "polluted-10min", "clean-10min",
              "summer-15min", "tropics-20min", "winter60-10min", "polluted-48h-30min", "clean-48h-5min"]
ATOL_SCENARIOS = ["diurnal-10min", "polluted-10min", "tropics-20min", "winter60-10min"]
LABEL = {"alg": "Algorithm", "alg-a20": r"Algorithm, atol $10^{-20}$", "alg-a15": r"Algorithm, atol $10^{-15}$",
         "restart-p": "Restarting loop", "restart": r"Restarting loop, atol $10^{-15}$",
         "restart-a20": r"Restarting loop, atol $10^{-20}$", "be": "Backward Euler",
         "be-tight": r"Backward Euler, atol $10^{-15}$", "alg-elementary": "elementary rule",
         "alg-reset": "restarted, not rounded", "alg-noquant": "carried, not rounded",
         "alg-nomemory": "filter memory reset", "restart-a13": r"Restarting loop, atol $10^{-13}$"}
LEVEL_CONFIGS = ("alg", "alg-a20", "alg-a15", "restart-p", "restart", "restart-a20")
# multiply-adds per cell on TS1 (Section 3): factorization, solve, Jacobian, forcing
WEIGHTS = {"lu": 10.187, "solves": 2.520, "jac": 4.785, "rhs": 3.189}


def work(a):
    """Thousands of multiply-adds per interval and cell."""
    return sum(w * a[k] for k, w in WEIGHTS.items())


def load(path):
    rows = []
    if not Path(path).exists():
        return rows
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            for k in ("rtol", "atol", "digits", "worst_digits"):
                r[k] = float(r[k])
            for k in ("cells", "intervals", "failed_calls", "steps", "accepted", "rejected", "rhs", "jac", "lu",
                      "solves", "newton", "failures", "unconverged"):
                r[k] = int(r[k])
            for k in ("gate", "level2", "level3"):
                r[k] = r[k] == "True"
            rows.append(r)
    return rows


def summarize(rs):
    """A summary of runs of one configuration over several scenarios."""
    groups = math.ceil(rs[0]["cells"] / 32)
    ints = sum(r["intervals"] for r in rs)
    per = lambda k: sum(r[k] for r in rs) / ints / groups
    worst = min(rs, key=lambda r: r["worst_digits"])
    lowest = min(rs, key=lambda r: r["digits"])
    return {
        "scenarios": len(rs), "atol": rs[0]["atol"],
        "all_level2": all(r["level2"] for r in rs), "all_level3": all(r["level3"] for r in rs),
        "n_level2": sum(r["level2"] for r in rs), "n_level3": sum(r["level3"] for r in rs),
        "mean_digits": sum(r["digits"] for r in rs) / len(rs),
        "min_digits": lowest["digits"], "min_digits_scenario": lowest["scenario"],
        "min_worst": worst["worst_digits"], "min_worst_species": worst["worst_species"],
        "min_worst_scenario": worst["scenario"],
        "steps": per("steps"), "rejected": per("rejected"), "lu": per("lu"), "solves": per("solves"),
        "rhs": per("rhs"), "jac": per("jac"), "newton": per("newton"), "failures": per("failures"),
        "failed_calls": sum(r["failed_calls"] for r in rs), "unconverged": sum(r["unconverged"] for r in rs),
        "rows": {r["scenario"]: r for r in rs},
    }


def aggregate(rows, scenarios=None):
    """(config, rtol) -> summary over the scenarios present (or the given ones)."""
    by = defaultdict(list)
    for r in rows:
        if scenarios is None or r["scenario"] in scenarios:
            by[(r["config"], r["rtol"])].append(r)
    out = {}
    for key, rs in by.items():
        out[key] = {"config": key[0], "rtol": key[1], **summarize(rs)}
    return out


def fmt(x, d=2):
    return f"{x:.{d}f}"


def tex_pow(x):
    m, e = f"{x:.0e}".split("e")
    e = int(e)
    return f"$10^{{{e}}}$" if m == "1" else f"${m}\\times10^{{{e}}}$"


def species(s):
    return s.replace("_", r"\_")


def at_digits(curve, target):
    """Steps at a mean digit count, interpolated linearly in log(steps) between the
    bracketing tolerances of one configuration; None outside the curve."""
    pts = sorted((a["mean_digits"], a["steps"]) for a in curve)
    for (d0, s0), (d1, s1) in zip(pts, pts[1:]):
        if d0 <= target <= d1 and d1 > d0:
            w = (target - d0) / (d1 - d0)
            return math.exp((1 - w) * math.log(s0) + w * math.log(s1))
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", default=str(ROOT / "results/ts1/runs.csv"))
    ap.add_argument("--atol", default=str(ROOT / "results/ts1/atol.csv"))
    ap.add_argument("--be-substeps", default=str(ROOT / "results/ts1/be_substeps.csv"))
    ap.add_argument("--paper", default=str(ROOT / "tables"))
    args = ap.parse_args()
    rows = load(args.runs)
    agg = {k: a for k, a in aggregate(rows).items() if a["scenarios"] == 11}
    paper = Path(args.paper)
    paper.mkdir(parents=True, exist_ok=True)
    md = ["# TS1 summary", "",
          f"{len(rows)} runs on 64 cells; counts per interval and per cell over the eleven scenarios.", ""]

    # ---- every configuration ---------------------------------------------------------
    md += ["| config | rtol | atol | L2 | L3 | mean D | lowest D (scenario) | lowest species | steps | rej | LU | "
           "solves | rhs | newton | failed |", "|" + "---|" * 15]
    for key in sorted(agg, key=lambda k: (k[0], -k[1])):
        a = agg[key]
        md.append(f"| {a['config']} | {a['rtol']:g} | {a['atol']:g} | {a['n_level2']} | {a['n_level3']} | "
                  f"{a['mean_digits']:.2f} | {a['min_digits']:.2f} ({a['min_digits_scenario']}) | "
                  f"{a['min_worst_species']} {a['min_worst']:.2f} ({a['min_worst_scenario']}) | {a['steps']:.2f} | "
                  f"{a['rejected']:.2f} | {a['lu']:.2f} | {a['solves']:.2f} | {a['rhs']:.2f} | {a['newton']:.2f} | "
                  f"{a['failed_calls']} |")

    # ---- level points: the loosest rtol reaching each level on every scenario ----------
    level_rows = []
    for L in (2, 3):
        for config in LEVEL_CONFIGS:
            cands = [a for k, a in agg.items() if k[0] == config and a[f"all_level{L}"]]
            if cands:
                level_rows.append((L, config, max(cands, key=lambda a: a["rtol"])))
    md += ["", "## Loosest relative tolerance reaching each level on all eleven scenarios", "",
           "| level | config | rtol | mean D | lowest D | lowest species | steps | rej | LU | solves | rhs |",
           "|" + "---|" * 11]
    for L, config, a in level_rows:
        md.append(f"| {L} | {config} | {a['rtol']:g} | {a['mean_digits']:.2f} | {a['min_digits']:.2f} | "
                  f"{a['min_worst_species']} {a['min_worst']:.2f} | {a['steps']:.2f} | {a['rejected']:.2f} | "
                  f"{a['lu']:.2f} | {a['solves']:.2f} | {a['rhs']:.2f} | work {work(a):.1f} |")

    # ---- steps at equal mean digits ------------------------------------------------------
    curves = defaultdict(list)
    for (config, _), a in agg.items():
        curves[config].append(a)
    targets = (3.0, 3.5, 4.0, 4.5)
    md += ["", "## Steps per interval at equal mean digits (interpolated)", "",
           "| config | " + " | ".join(f"D = {t}" for t in targets) + " |", "|" + "---|" * (len(targets) + 1)]
    equal = {}
    for config in sorted(curves, key=lambda c: (c not in LEVEL_CONFIGS, c)):
        equal[config] = [at_digits(curves[config], t) for t in targets]
        md.append(f"| {config} | " + " | ".join("--" if v is None else f"{v:.2f}" for v in equal[config]) + " |")
    # rejections at equal mean digits, for the controllers
    md += ["", "## Rejections per interval at equal mean digits", "",
           "| config | " + " | ".join(f"D = {t}" for t in targets) + " |", "|" + "---|" * (len(targets) + 1)]
    for config in sorted(curves):
        rej = [at_digits([{**a, "steps": max(a["rejected"], 1e-6)} for a in curves[config]], t) for t in targets]
        md.append(f"| {config} | " + " | ".join("--" if v is None else f"{v:.2f}" for v in rej) + " |")

    # ---- LaTeX: the level table ---------------------------------------------------
    lines = [r"\begin{tabular}{llllrrr}", r"\toprule",
             r"Level & Solver & rtol & Mean $D$; lowest $D$; worst species & Steps & Rej. & Work \\",
             r"\midrule"]
    last = None
    for L, config, a in level_rows:
        if last is not None and L != last:
            lines.append(r"\midrule")
        last = L
        lines.append(f"{L} & {LABEL[config]} & {tex_pow(a['rtol'])} & {a['mean_digits']:.2f}; "
                     f"{a['min_digits']:.2f}; {species(a['min_worst_species'])} {a['min_worst']:.2f} & "
                     f"{fmt(a['steps'])} & {fmt(a['rejected'])} & {work(a):.0f} \\\\")
    be = agg.get(("be", 1e-6))
    if be:
        lines += [r"\midrule",
                  f"none & {LABEL['be']} & {tex_pow(1e-6)} & {be['mean_digits']:.2f}; {be['min_digits']:.2f}; "
                  f"{species(be['min_worst_species'])} {be['min_worst']:.2f} & {fmt(be['steps'])} & "
                  f"{fmt(be['rejected'])} & {work(be):.0f} \\\\"]
    lines += [r"\bottomrule", r"\end{tabular}"]
    (paper / "ts1_levels.tex").write_text("\n".join(lines) + "\n")

    # ---- LaTeX: steps at equal digits ---------------------------------------------------
    lines = [r"\begin{tabular}{l" + "r" * len(targets) + "}", r"\toprule",
             r"Solver & " + " & ".join(f"$D = {t}$" for t in targets) + r" \\", r"\midrule"]
    for config in LEVEL_CONFIGS:
        lines.append(f"{LABEL[config]} & " + " & ".join("--" if v is None else fmt(v) for v in equal[config]) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (paper / "ts1_equal_digits.tex").write_text("\n".join(lines) + "\n")

    # ---- LaTeX: the sweep for the algorithm and the restarting loop -------------------
    lines = [r"\begin{tabular}{llrrllrrr}", r"\toprule",
             r"Solver & rtol & L2 & L3 & Mean $D$ & Lowest $D$; worst species & Steps & Rej. & Work \\",
             r"\midrule"]
    for config in ("alg", "restart-p"):
        for key in sorted((k for k in agg if k[0] == config), key=lambda k: -k[1]):
            a = agg[key]
            lines.append(f"{LABEL[config]} & {tex_pow(a['rtol'])} & {a['n_level2']} & {a['n_level3']} & "
                         f"{a['mean_digits']:.2f} & {a['min_digits']:.2f}; {species(a['min_worst_species'])} "
                         f"{a['min_worst']:.2f} & {fmt(a['steps'])} & {fmt(a['rejected'])} & {work(a):.0f} \\\\")
        lines.append(r"\midrule")
    lines[-1] = r"\bottomrule"
    lines.append(r"\end{tabular}")
    (paper / "ts1_sweep.tex").write_text("\n".join(lines) + "\n")

    # ---- LaTeX: the changes taken apart, at matched tolerances ---------------------
    variants = ("alg", "alg-nomemory", "alg-elementary", "alg-noquant", "alg-reset", "alg-a20", "alg-a15")
    names = {"alg": "the algorithm", "alg-a20": r"atol $10^{-20}$", "alg-a15": r"atol $10^{-15}$"}
    lines = [r"\begin{tabular}{llrrrll}", r"\toprule",
             r"rtol & Variant & Steps & Rej. & Work & Mean $D$ & Lowest $D$; worst species \\", r"\midrule"]
    for rtol in (1e-2, 2e-3, 5e-4):
        for config in variants:
            a = agg.get((config, rtol))
            if a is None:
                continue
            lines.append(f"{tex_pow(rtol)} & {names.get(config, LABEL[config])} & {fmt(a['steps'])} & "
                         f"{fmt(a['rejected'])} & {work(a):.0f} & {a['mean_digits']:.2f} & {a['min_digits']:.2f}; "
                         f"{species(a['min_worst_species'])} {a['min_worst']:.2f} \\\\")
        lines.append(r"\midrule")
    lines[-1] = r"\bottomrule"
    lines.append(r"\end{tabular}")
    (paper / "ts1_ablation.tex").write_text("\n".join(lines) + "\n")

    # ---- LaTeX: per-scenario digits and steps at the level-3 points, and backward Euler ----
    pick = {c: a for L, c, a in level_rows if L == 3}
    base = min((c for c in ("restart-p", "restart", "restart-a20") if c in pick), key=lambda c: pick[c]["steps"],
               default=None)
    cols = [("alg", pick.get("alg")), (base, pick.get(base) if base else None)]
    be = agg.get(("be", 1e-6))
    head = lambda c, a: {"alg": "Algorithm"}.get(c, "Restarting loop" if c else "--")
    lines = [r"\begin{tabular}{lrrrrrr}", r"\toprule",
             r"& \multicolumn{2}{c}{" + head("alg", cols[0][1]) + r"} & \multicolumn{2}{c}{"
             + (head(base, cols[1][1]) if base else "--") + r"} & \multicolumn{2}{c}{Backward Euler} \\",
             r"Scenario & $D$ (worst) & Steps & $D$ (worst) & Steps & $D$ (worst) & Newton \\", r"\midrule"]
    for s in SCEN_ORDER:
        cells = []
        for _, a in cols:
            r = a["rows"].get(s) if a else None
            groups = math.ceil(r["cells"] / 32) if r else 1
            cells += ([f"{r['digits']:.2f} ({r['worst_digits']:.2f})", fmt(r['steps'] / r['intervals'] / groups)]
                      if r else ["--", "--"])
        r = be["rows"].get(s) if be else None
        groups = math.ceil(r["cells"] / 32) if r else 1
        cells += ([f"{r['digits']:.2f} ({r['worst_digits']:.2f})", fmt(r['newton'] / r['intervals'] / groups)]
                  if r else ["--", "--"])
        lines.append(f"\\texttt{{{s}}} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}"]
    (paper / "ts1_scenarios.tex").write_text("\n".join(lines) + "\n")

    # ---- LaTeX: the absolute-tolerance sweep on four scenarios -------------------------
    sweep = load(args.atol) + [r for r in rows if r["config"] in ("alg", "alg-a15", "alg-a20")]
    four = aggregate(sweep, ATOL_SCENARIOS)
    lines = [r"\begin{tabular}{llrlr}", r"\toprule",
             r"rtol & atol & Mean $D$ & Lowest $D_i$, species & Steps \\", r"\midrule"]
    md += ["", "## Absolute-tolerance sweep (four scenarios)", "",
           "| rtol | atol | mean D | lowest species | steps | rej |", "|---|---|---|---|---|---|"]
    for rtol in (1e-2, 2e-3):
        entries = [a for (c, r), a in four.items() if r == rtol and c.startswith("alg") and a["scenarios"] == 4]
        rule = [a["atol"] for a in entries if a["config"] == "alg"]
        entries = [a for a in entries if a["config"] == "alg" or not any(math.isclose(a["atol"], x) for x in rule)]
        for a in sorted(entries, key=lambda a: -a["atol"]):
            tag = r" $= \rtol \cdot a$" if a["config"] == "alg" else ""
            lines.append(f"{tex_pow(rtol)} & {tex_pow(a['atol'])}{tag} & {a['mean_digits']:.2f} & "
                         f"{a['min_worst']:.2f}, {species(a['min_worst_species'])} & {fmt(a['steps'])} \\\\")
            md.append(f"| {rtol:g} | {a['atol']:g} | {a['mean_digits']:.2f} | {a['min_worst_species']} "
                      f"{a['min_worst']:.2f} | {a['steps']:.2f} | {a['rejected']:.2f} |")
        lines.append(r"\midrule")
    lines[-1] = r"\bottomrule"
    lines.append(r"\end{tabular}")
    (paper / "ts1_atol.tex").write_text("\n".join(lines) + "\n")
    rs13 = aggregate([r for r in load(args.atol) if r["config"] == "restart-a13"] +
                     [r for r in rows if r["config"] in ("restart", "restart-p")], ATOL_SCENARIOS)
    md += ["", "## Restarting loop on the same four scenarios", "", "| config | rtol | atol | mean D | lowest | steps |",
           "|---|---|---|---|---|---|"]
    for (c, r), a in sorted(rs13.items(), key=lambda kv: (kv[0][0], -kv[0][1])):
        if a["scenarios"] == 4:
            md.append(f"| {c} | {r:g} | {a['atol']:g} | {a['mean_digits']:.2f} | {a['min_worst_species']} "
                      f"{a['min_worst']:.2f} | {a['steps']:.2f} |")

    # ---- LaTeX: backward Euler sub-stepped, against the algorithm, on three scenarios ----
    be_path = Path(args.be_substeps)
    if be_path.exists():
        with open(be_path, newline="") as f:
            be_rows = list(csv.DictReader(f))
        lines = [r"\begin{tabular}{llrlrr}", r"\toprule",
                 r"Scenario & Solver & Sub-steps & $D$ (worst species) & Newton & Work \\", r"\midrule"]
        md += ["", "## Backward Euler with sub-steps, against the algorithm", "",
               "| scenario | solver | substeps | D | worst | newton | work |", "|---|---|---|---|---|---|---|"]
        for s_name in ("tier1-30s", "diurnal-10min", "diurnal-30min"):
            rs = sorted((r for r in be_rows if r["scenario"] == s_name), key=lambda r: int(r["substeps"]))
            for r in rs:
                lines.append(f"\\texttt{{{s_name}}} & Backward Euler & {r['substeps']} & {float(r['digits']):.2f} "
                             f"({float(r['worst_digits']):.2f}) & {float(r['newton_per_interval']):.0f} & "
                             f"{float(r['work_kmadd']):.0f} \\\\")
                md.append(f"| {s_name} | BE | {r['substeps']} | {float(r['digits']):.2f} | "
                          f"{r['worst_species']} {float(r['worst_digits']):.2f} | {r['newton_per_interval']} | "
                          f"{r['work_kmadd']} |")
            alg = [r for r in rows if r["config"] == "alg" and r["scenario"] == s_name and r["level2"]]
            if alg:
                a = summarize([max(alg, key=lambda r: r["rtol"])])
                lines.append(f"& Algorithm, {tex_pow(a['rows'][s_name]['rtol'])} & -- & {a['mean_digits']:.2f} "
                             f"({a['min_worst']:.2f}) & -- & {work(a):.0f} \\\\")
                md.append(f"| {s_name} | alg {a['rows'][s_name]['rtol']:g} | -- | {a['mean_digits']:.2f} | "
                          f"{a['min_worst_species']} {a['min_worst']:.2f} | -- | {work(a):.1f} |")
            lines.append(r"\midrule")
        lines[-1] = r"\bottomrule"
        lines.append(r"\end{tabular}")
        (paper / "ts1_be_substeps.tex").write_text("\n".join(lines) + "\n")

    # ---- worst species ----------------------------------------------------------------
    counts = defaultdict(int)
    for r in rows:
        if r["config"] in ("alg", "restart-p"):
            counts[r["worst_species"]] += 1
    (ROOT / "results/ts1").mkdir(parents=True, exist_ok=True)
    (ROOT / "results/ts1/worst_species.json").write_text(
        json.dumps(dict(sorted(counts.items(), key=lambda kv: -kv[1])), indent=1))
    (ROOT / "results/ts1/summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"... {len(agg)} configurations; tables written to {paper}")


if __name__ == "__main__":
    main()
