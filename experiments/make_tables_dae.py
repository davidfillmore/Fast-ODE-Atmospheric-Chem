"""LaTeX tables for the small problems and the constrained species.

    python3 experiments/make_tables_dae.py [--paper tables] [--results results]
"""
from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"


def read(name, sub="dae"):
    with open(RESULTS / sub / f"{name}.csv", newline="") as f:
        return list(csv.DictReader(f))


def slope(rows, key, n1, n2, field):
    e = {int(r["N"]): float(r[field]) for r in rows if all(r[k] == v for k, v in key.items())}
    return math.log(e[n1] / e[n2]) / math.log(n2 / n1)


def table(lines):
    return "\n".join(lines) + "\n"


def pow10(x):
    """$10^{-3}$ for 1e-3, $3\\times10^{-4}$ for 3e-4."""
    m, e = f"{float(x):.0e}".split("e")
    return f"$10^{{{int(e)}}}$" if m == "1" else f"${m}\\times10^{{{int(e)}}}$"


def species(name):
    return name.replace("_", "\\_")


def main():
    global RESULTS
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper", default=str(ROOT / "tables"))
    ap.add_argument("--results", default=str(RESULTS))
    args = ap.parse_args()
    RESULTS = Path(args.results)
    out = Path(args.paper)
    out.mkdir(parents=True, exist_ok=True)

    # ---- small problems: the three solvers at rtol 1e-3 -------------------------------
    small = read("runs", "small")
    pick = defaultdict(dict)
    for r in small:
        if r["solver"] == "backward-euler" or float(r["rtol"]) == 1e-3:
            pick[r["problem"]][r["solver"]] = r
    names = {"robertson": "Robertson", "chapman": "Chapman", "tropospheric": "Tropospheric",
             "tropospheric-diurnal": "Tropospheric, diurnal", "triples": "Triples ($N = 10$)", "hires": "HIRES"}
    L = [r"\begin{tabular}{lrrrrrrrr}", r"\toprule",
         r"& \multicolumn{4}{c}{Algorithm} & \multicolumn{2}{c}{Restarting loop} & \multicolumn{2}{c}{Backward Euler} \\",
         r"\cmidrule(lr){2-5} \cmidrule(lr){6-7} \cmidrule(lr){8-9}",
         r"& & \multicolumn{3}{c}{Steps at atol} & & & & \\",
         r"Problem & $D$ & $\rtol \cdot a$ & $a/10$ & $10^{-6} a$ & $D$ & Steps & $D$ & Newton \\", r"\midrule"]
    for p, label in names.items():
        c = pick[p]
        per = lambda sv, k: float(c[sv][k]) / int(c[sv]["intervals"])
        L.append(f"{label} & {float(c['algorithm']['digits']):.2f} & {per('algorithm', 'steps'):.2f} & "
                 f"{per('algorithm-a/10', 'steps'):.2f} & {per('algorithm-tiny', 'steps'):.2f} & "
                 f"{float(c['restarting']['digits']):.2f} & {per('restarting', 'steps'):.2f} & "
                 f"{float(c['backward-euler']['digits']):.2f} & {per('backward-euler', 'newton'):.2f} \\\\")
    L += [r"\bottomrule", r"\end{tabular}"]
    (out / "small.tex").write_text(table(L))

    # ---- orders ------------------------------------------------------------------------
    orders = read("orders")
    pr = read("prothero_robinson")
    qde = read("qde")
    L = [r"\begin{tabular}{llcc}", r"\toprule", r"Problem & Tableau & Differential & Algebraic \\", r"\midrule"]
    for prob, label in (("circle", r"$\dot y = z$, $0 = y^2 + z^2 - 1$"), ("coupled", "the coupled problem")):
        for tab in ("ros3", "rodas3"):
            d = slope(orders, {"problem": prob, "tableau": tab}, 128, 256, "err_differential")
            a = slope(orders, {"problem": prob, "tableau": tab}, 128, 256, "err_algebraic")
            L.append(f"{label} & \\method{{{tab}}} & {d:.2f} & {a:.2f} \\\\")
    L.append(r"\midrule")
    for tab in ("ros3", "rodas3"):
        s = slope(pr, {"tableau": tab, "lambda": "-1000000.0"}, 512, 1024, "error")
        L.append(f"Prothero--Robinson, $\\lambda = -10^6$ & \\method{{{tab}}} & {s:.2f} & -- \\\\")
    for prob in ("qde2", "qde3"):
        for tab in ("ros3", "rodas3"):
            s = slope(qde, {"problem": prob, "tableau": tab}, 64, 128, "error")
            L.append(f"{prob.upper().replace('QDE', 'QDE-')} & \\method{{{tab}}} & {s:.2f} & -- \\\\")
    L += [r"\bottomrule", r"\end{tabular}"]
    (out / "dae_orders.tex").write_text(table(L))

    # ---- Robertson and the tropospheric mechanism: forms and loops -----------------------
    rob = read("robertson")
    tro = read("tropospheric")
    ctl = read("rows_control")
    L = [r"\begin{tabular}{llrrr}", r"\toprule", r"Problem & Formulation and solver & $D$ & Steps & Failed calls \\",
         r"\midrule"]
    for r in rob:
        if float(r["k2"]) == 3e7 and float(r["rtol"]) == 1e-4 and r["solver"] == "algorithm":
            L.append(f"Robertson & {r['form'].upper()}, \\method{{{r['tableau']}}}, algorithm & {float(r['digits_AC']):.2f} & "
                     f"{r['steps']} & {r['failed_calls']} \\\\")
    for r in ctl:
        if r["problem"] == "robertson" and float(r["rtol"]) == 1e-4 and r["solver"] != "algorithm":
            L.append(f"Robertson & DAE, \\method{{rodas3}}, {r['solver']} & {float(r['digits']):.2f} & {r['steps']} & "
                     f"{r['failed_calls']} \\\\")
    L.append(r"\midrule")
    for r in tro:
        if float(r["j"]) == 1.0 and float(r["rtol"]) == 1e-4:
            L.append(f"Tropospheric & {r['form'].upper()}, \\method{{rodas3}}, {r['solver']} & {float(r['digits']):.2f} & "
                     f"{r['steps']} & {r['failed_calls']} \\\\")
    for r in ctl:
        if r["problem"] == "tropospheric" and float(r["rtol"]) == 1e-4 and r["solver"] == "restarting, embedded rows":
            L.append(f"Tropospheric & DAE, \\method{{rodas3}}, {r['solver']} & {float(r['digits']):.2f} & {r['steps']} & "
                     f"{r['failed_calls']} \\\\")
    L += [r"\bottomrule", r"\end{tabular}"]
    (out / "dae_forms.tex").write_text(table(L))

    # ---- sweeps in the stiffness and photolysis --------------------------------------------
    L = [r"\begin{tabular}{lrrrrr}", r"\toprule", r"Problem & Parameter & Steps, ODE & Steps, DAE & $D$, ODE & $D$, DAE \\",
         r"\midrule"]
    by = defaultdict(dict)
    for r in rob:
        if float(r["rtol"]) == 1e-4 and r["solver"] == "algorithm" and r["tableau"] == "rodas3":
            by[("Robertson", "$k_2 = " + f"{float(r['k2']):.0e}".replace("e+0", r"\times10^{") + "}$")][r["form"]] = (
                int(r["steps"]), float(r["digits_AC"]))
    for r in read("chapman"):
        if float(r["rtol"]) == 1e-4:
            by[("Chapman", f"scale ${float(r['scale']):g}$")][r["form"]] = (int(r["steps"]), float(r["digits"]))
    for r in tro:
        if float(r["rtol"]) == 1e-4 and r["solver"] == "algorithm":
            by[("Tropospheric", f"$j = {float(r['j']):g}$")][r["form"]] = (int(r["steps"]), float(r["digits"]))
    trip = defaultdict(dict)
    for r in read("triples"):
        trip[r["form"]][r["N"]] = (int(r["steps"]), float(r["digits"]))
    for form, v in trip.items():
        if len(set(v.values())) == 1:      # the same at every N: one row
            by[("Triples", "$N = " + ", ".join(v) + "$")][form] = next(iter(v.values()))
        else:
            for N, x in v.items():
                by[("Triples", f"$N = {N}$")][form] = x
    for (p, par), v in by.items():
        if "ode" in v and "dae" in v:
            L.append(f"{p} & {par} & {v['ode'][0]} & {v['dae'][0]} & {v['ode'][1]:.2f} & {v['dae'][1]:.2f} \\\\")
    L += [r"\bottomrule", r"\end{tabular}"]
    (out / "dae_sweeps.tex").write_text(table(L))

    # ---- the projection sweep ---------------------------------------------------------
    proj = read("projection")
    agg = defaultdict(lambda: [0, 0, 0, 0])
    for r in proj:
        a = agg[r["problem"]]
        a[0] += 1
        a[1] += (r["weighted_status"] != "converged") or (r["weighted_false_conv"] == "True")
        a[2] += r["residual_false_conv"] == "True"
        a[3] += r["residual_false_fail"] == "True"
    L = [r"\begin{tabular}{lrrrr}", r"\toprule",
         r"Problem & Cases & Weighted test: wrong & Residual test: false convergence & false failure \\", r"\midrule"]
    tot = [0, 0, 0, 0]
    for p, a in agg.items():
        label = {"cold-start": "Cold start"}.get(p, p.capitalize())
        L.append(f"{label} & {a[0]} & {a[1]} & {a[2]} & {a[3]} \\\\")
        tot = [x + y for x, y in zip(tot, a)]
    L += [r"\midrule", f"All & {tot[0]} & {tot[1]} & {tot[2]} & {tot[3]} \\\\", r"\bottomrule", r"\end{tabular}"]
    (out / "projection.tex").write_text(table(L))

    # ---- TS1 with constrained families ----------------------------------------------
    fam = read("ts1_families")
    L = [r"\begin{tabular}{lllrrrrr}", r"\toprule",
         r"rtol & Formulation & Solver & $D$ & Worst species & Steps & LU & Projections \\", r"\midrule"]
    for r in fam:
        form = {"ode": "ODE"}.get(r["form"], r["form"].replace("dae-", "DAE, ") + " constrained")
        L.append(f"{pow10(r['rtol'])} & {form}, \\method{{{r['tableau']}}} & {r['solver']} & {float(r['digits']):.2f} & "
                 f"{species(r['worst'])} {float(r['worst_digits']):.2f} & {float(r['steps_per_interval']):.2f} & "
                 f"{float(r['lu_per_interval']):.2f} & {float(r['proj_passes_per_interval']):.2f} \\\\")
    L += [r"\bottomrule", r"\end{tabular}"]
    (out / "ts1_families.tex").write_text(table(L))
    print("written:", sorted(p.name for p in out.glob("*.tex")))


if __name__ == "__main__":
    main()
